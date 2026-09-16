#!/usr/bin/env python

"""
Compare two pipeline runs that differ only in the aligner (e.g. STAR vs
rustar-aligner) and report whether outputs and runtimes agree.

Inputs are the two --outdir directories (typically results-star/ and
results-rustar/). The script reads:

- {outdir}/star_salmon/log/*.Log.final.out         per-sample STAR-format
                                                   alignment summary
- {outdir}/star_salmon/salmon.merged.gene_tpm.tsv  merged TPM matrix
- {outdir}/star_salmon/salmon.merged.gene_counts.tsv  merged counts matrix
- {outdir}/pipeline_info/execution_trace_*.txt     Nextflow trace report

Outputs a JSON summary (machine-readable) and a Markdown table
(reviewer-readable, for pasting into a PR body).
"""

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

PASS_TPM_PEARSON = 0.999
PASS_PERCENT_MAPPED_PP = 0.5  # percentage-points tolerance

ALIGN_PROCESSES = {
    "STAR_ALIGN",
    "RUSTAR_ALIGN",
    "SENTIEON_STAR_ALIGN",
    "PARABRICKS_RNA_FQ2BAM",
}
GENOMEGEN_PROCESSES = {"STAR_GENOMEGENERATE", "RUSTAR_GENOMEGENERATE"}


def parse_log_final(path: Path) -> dict:
    """Parse a STAR-format Log.final.out into a flat dict keyed by metric name."""
    out = {}
    with path.open() as fh:
        for line in fh:
            if "|" not in line:
                continue
            key, _, value = line.partition("|")
            key = key.strip()
            value = value.strip().rstrip("%")
            try:
                out[key] = float(value)
            except ValueError:
                out[key] = value
    return out


def collect_log_finals(outdir: Path) -> dict:
    """Map sample id -> parsed Log.final.out for the star_salmon path."""
    results = {}
    for log_path in (outdir / "star_salmon" / "log").glob("*.Log.final.out"):
        sample = log_path.name.split(".")[0]
        results[sample] = parse_log_final(log_path)
    return results


def read_tsv_matrix(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    """Read a salmon merged TSV. First column is the feature id; any extra
    leading non-numeric columns (e.g. gene_name) are discarded; the remaining
    columns are sample values."""
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        first_data = next(reader)
        sample_start = 1
        for i in range(1, len(first_data)):
            try:
                float(first_data[i])
                sample_start = i
                break
            except ValueError:
                continue
        samples = header[sample_start:]
        ids = [first_data[0]]
        rows = [[float(v) for v in first_data[sample_start:]]]
        for row in reader:
            if not row:
                continue
            ids.append(row[0])
            rows.append([float(v) for v in row[sample_start:]])
    return samples, ids, rows


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n == 0:
        return float("nan")
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    cov = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    denom = math.sqrt(var_x * var_y)
    return cov / denom if denom else float("nan")


def spearman(xs: list[float], ys: list[float]) -> float:
    # Spearman = Pearson on ranks; ties get the average rank.
    return pearson(rank(xs), rank(ys))


def rank(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda iv: iv[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j + 1 < len(indexed) and indexed[j + 1][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[indexed[k][0]] = avg_rank
        i = j + 1
    return ranks


def compare_matrix(path_a: Path, path_b: Path) -> dict:
    samples_a, ids_a, rows_a = read_tsv_matrix(path_a)
    samples_b, ids_b, rows_b = read_tsv_matrix(path_b)
    if samples_a != samples_b:
        return {"error": f"sample columns differ: {samples_a} vs {samples_b}"}
    if ids_a != ids_b:
        return {"error": "feature ids differ between matrices"}

    per_sample = {}
    for col_idx, sample in enumerate(samples_a):
        xs = [row[col_idx] for row in rows_a]
        ys = [row[col_idx] for row in rows_b]
        per_sample[sample] = {
            "pearson": pearson(xs, ys),
            "spearman": spearman(xs, ys),
            "n_features": len(xs),
        }
    return per_sample


def parse_trace(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def trace_summary(trace_path: Path, process_filter: set[str]) -> dict:
    rows = parse_trace(trace_path)
    by_process = {}
    for row in rows:
        # Nextflow trace process names look like "NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (sample)".
        # The trace file column is "name" in modern Nextflow; older versions used "process".
        name = row.get("name") or row.get("process", "")
        base = name.split(":")[-1].split(" ")[0]
        if base not in process_filter:
            continue
        by_process.setdefault(base, []).append(row)

    out = {}
    for proc, proc_rows in by_process.items():
        durations_ms = [parse_duration_ms(r.get("realtime", "")) for r in proc_rows]
        peak_bytes = [parse_size_bytes(r.get("peak_rss", "")) for r in proc_rows]
        out[proc] = {
            "n_tasks": len(proc_rows),
            "wall_seconds_median": median([d / 1000 for d in durations_ms if d]),
            "wall_seconds_total": sum(d / 1000 for d in durations_ms if d),
            "peak_rss_gb_max": max((b for b in peak_bytes if b), default=0) / 1024**3,
        }
    return out


_DURATION_RE = re.compile(r"(\d+(?:\.\d+)?)\s*(ms|s|m|h)")


def parse_duration_ms(text: str) -> float:
    text = (text or "").strip()
    if not text or text == "-":
        return 0.0
    total = 0.0
    for value, unit in _DURATION_RE.findall(text):
        v = float(value)
        if unit == "ms":
            total += v
        elif unit == "s":
            total += v * 1000
        elif unit == "m":
            total += v * 60 * 1000
        elif unit == "h":
            total += v * 3600 * 1000
    return total


_SIZE_RE = re.compile(r"([\d.]+)\s*([KMGT]?)B?")


def parse_size_bytes(text: str) -> float:
    text = (text or "").strip()
    if not text or text == "-" or text == "0":
        return 0.0
    m = _SIZE_RE.match(text)
    if not m:
        return 0.0
    value = float(m.group(1))
    unit = m.group(2)
    factor = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}[unit]
    return value * factor


def median(xs: list[float]) -> float:
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return 0.0
    mid = len(xs) // 2
    if len(xs) % 2:
        return xs[mid]
    return (xs[mid - 1] + xs[mid]) / 2


def latest_trace(outdir: Path) -> Path | None:
    candidates = sorted((outdir / "pipeline_info").glob("execution_trace_*.txt"))
    return candidates[-1] if candidates else None


def render_markdown(report: dict) -> str:
    lines = []
    lines.append("# Aligner comparison")
    lines.append("")
    lines.append(f"- STAR run: `{report['paths']['star']}`")
    lines.append(f"- rustar run: `{report['paths']['rustar']}`")
    lines.append("")

    lines.append("## % Uniquely mapped reads")
    lines.append("")
    lines.append("| Sample | STAR | rustar | Δ (pp) | Pass |")
    lines.append("|---|---|---|---|---|")
    for sample, m in sorted(report["mapping"].items()):
        delta = m["delta_pp"]
        passed = "✅" if m["pass"] else "❌"
        lines.append(
            f"| {sample} | {m['star']:.2f} | {m['rustar']:.2f} | {delta:+.2f} | {passed} |"
        )
    lines.append("")

    lines.append("## Salmon merged matrices (per-sample Pearson)")
    lines.append("")
    lines.append("| Sample | gene_tpm | gene_counts | Pass |")
    lines.append("|---|---|---|---|")
    for sample in sorted(report["matrices"].get("gene_tpm", {})):
        tpm = report["matrices"]["gene_tpm"][sample]["pearson"]
        cnt = report["matrices"].get("gene_counts", {}).get(sample, {}).get("pearson", float("nan"))
        passed = "✅" if tpm >= PASS_TPM_PEARSON else "❌"
        lines.append(f"| {sample} | {tpm:.6f} | {cnt:.6f} | {passed} |")
    lines.append("")

    lines.append("## Trace timings")
    lines.append("")
    lines.append("| Process | n | Wall median (s) STAR → rustar | Peak RSS (GB) STAR → rustar |")
    lines.append("|---|---|---|---|")
    for proc in sorted(set(report["trace"]["star"]) | set(report["trace"]["rustar"])):
        a = report["trace"]["star"].get(proc, {})
        b = report["trace"]["rustar"].get(proc, {})
        n = b.get("n_tasks", a.get("n_tasks", 0))
        wa = a.get("wall_seconds_median", 0)
        wb = b.get("wall_seconds_median", 0)
        ra = a.get("peak_rss_gb_max", 0)
        rb = b.get("peak_rss_gb_max", 0)
        lines.append(f"| {proc} | {n} | {wa:.1f} → {wb:.1f} | {ra:.2f} → {rb:.2f} |")
    lines.append("")

    summary = report["summary"]
    verdict = "✅ PASS" if summary["pass"] else "❌ FAIL"
    lines.append(f"## Overall: {verdict}")
    if not summary["pass"]:
        for reason in summary["failures"]:
            lines.append(f"- {reason}")
    lines.append("")
    return "\n".join(lines)


def build_report(star_dir: Path, rustar_dir: Path) -> dict:
    star_logs = collect_log_finals(star_dir)
    rustar_logs = collect_log_finals(rustar_dir)
    mapping = {}
    for sample in sorted(set(star_logs) & set(rustar_logs)):
        a = star_logs[sample].get("Uniquely mapped reads %", 0.0)
        b = rustar_logs[sample].get("Uniquely mapped reads %", 0.0)
        delta = b - a
        mapping[sample] = {
            "star": a,
            "rustar": b,
            "delta_pp": delta,
            "pass": abs(delta) <= PASS_PERCENT_MAPPED_PP,
        }

    matrices = {}
    for name in ("gene_tpm", "gene_counts"):
        a = star_dir / "star_salmon" / f"salmon.merged.{name}.tsv"
        b = rustar_dir / "star_salmon" / f"salmon.merged.{name}.tsv"
        if a.exists() and b.exists():
            matrices[name] = compare_matrix(a, b)

    star_trace = latest_trace(star_dir)
    rustar_trace = latest_trace(rustar_dir)
    process_filter = ALIGN_PROCESSES | GENOMEGEN_PROCESSES
    trace = {
        "star": trace_summary(star_trace, process_filter) if star_trace else {},
        "rustar": trace_summary(rustar_trace, process_filter) if rustar_trace else {},
    }

    failures = []
    for sample, m in mapping.items():
        if not m["pass"]:
            failures.append(
                f"% mapped delta {m['delta_pp']:+.2f} pp for {sample} exceeds ±{PASS_PERCENT_MAPPED_PP} pp"
            )
    for sample, m in matrices.get("gene_tpm", {}).items():
        if isinstance(m, dict) and "pearson" in m and m["pearson"] < PASS_TPM_PEARSON:
            failures.append(
                f"gene_tpm Pearson {m['pearson']:.6f} for {sample} below {PASS_TPM_PEARSON}"
            )

    return {
        "paths": {"star": str(star_dir), "rustar": str(rustar_dir)},
        "mapping": mapping,
        "matrices": matrices,
        "trace": trace,
        "summary": {"pass": not failures, "failures": failures},
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("star_outdir", type=Path, help="Pipeline outdir of the STAR run")
    parser.add_argument("rustar_outdir", type=Path, help="Pipeline outdir of the rustar run")
    parser.add_argument("--json", type=Path, help="Write JSON report here")
    parser.add_argument("--md", type=Path, help="Write Markdown report here")
    args = parser.parse_args(argv)

    if not args.star_outdir.is_dir() or not args.rustar_outdir.is_dir():
        parser.error("both outdir paths must exist")

    report = build_report(args.star_outdir, args.rustar_outdir)
    if args.json:
        args.json.write_text(json.dumps(report, indent=2))
    if args.md:
        args.md.write_text(render_markdown(report))
    if not args.json and not args.md:
        print(render_markdown(report))
    return 0 if report["summary"]["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
