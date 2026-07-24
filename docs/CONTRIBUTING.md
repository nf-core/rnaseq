---
title: Contributing
markdownPlugin: checklist
---

# `nf-core/rnaseq`: Contributing guidelines

Hi there!
Thanks for taking an interest in improving nf-core/rnaseq.

This page describes the recommended nf-core way to contribute to both nf-core/rnaseq and nf-core pipelines in general, including:

- [General contribution guidelines](#general-contribution-guidelines): common procedures or guides across all nf-core pipelines.
- [Pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines): procedures or guides specific to the development conventions of nf-core/rnaseq.

> [!NOTE]
> If you need help using or modifying nf-core/rnaseq, ask on the nf-core Slack [#rnaseq](https://nfcore.slack.com/channels/rnaseq) channel ([join our Slack here](https://nf-co.re/join/slack)).

## General contribution guidelines

### Contribution quick start

To contribute code to any nf-core pipeline:

- [ ] Ensure you have Nextflow, nf-core tools, and nf-test installed. See the [nf-core/tools repository](https://github.com/nf-core/tools) for instructions.
- [ ] Check whether a GitHub [issue](https://github.com/nf-core/rnaseq/issues) about your idea already exists. If an issue does not exist, create one so that others are aware you are working on it.
- [ ] [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/rnaseq repository](https://github.com/nf-core/rnaseq) to your GitHub account.
- [ ] Create a branch on your forked repository and make your changes following [pipeline conventions](#pipeline-contribution-conventions) (if applicable).
- [ ] To fix major bugs, name your branch `patch` and follow the [patch release](#patch-release) process.
- [ ] Update relevant documentation within the `docs/` folder, use nf-core/tools to update `nextflow_schema.json`, and update `CITATIONS.md`.
- [ ] Run and/or update tests. See [Testing](#testing) for more information.
- [ ] [Lint](#lint-tests) your code with nf-core/tools.
- [ ] Submit a pull request (PR) against the `dev` branch and request a review.

If you are not used to this workflow with Git, see the [GitHub documentation](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or [Git resources](https://try.github.io/) for more information.

## Use of AI and LLMs

The nf-core stance on the use of AI and LLMs is that humans are still ultimately responsible for their submitted code, regardless of the tools they use.

If you’re using AI tools, try to stick by these guidelines:

- Keep PRs as small and focused as possible
- Avoid any unnecessary changes, such as moving or refactoring code (unless that is the explicit intention of the PR)
- Review all generated code yourself before opening a PR, and ensure that you understand it
- Engage with the community review process and expect to make revisions

For more detail, see the the [blog post](https://nf-co.re/blog/2026/statement-on-ai) for a statement from the nf-core/core team.

### Getting help

For further information and help, see the [nf-core/rnaseq documentation](https://nf-co.re/rnaseq/usage) or ask on the nf-core [#rnaseq](https://nfcore.slack.com/channels/rnaseq) Slack channel ([join our Slack here](https://nf-co.re/join/slack)).

### GitHub Codespaces

You can contribute to nf-core/rnaseq without installing a local development environment on your machine by using [GitHub Codespaces](https://github.com/codespaces).

[GitHub Codespaces](https://github.com/codespaces) is an online developer environment that runs in your browser, complete with VS Code and a terminal.
Most nf-core repositories include a devcontainer configuration, which creates a GitHub Codespaces environment specifically for Nextflow development.
The environment includes pre-installed nf-core tools, Nextflow, and a few other helpful utilities via a Docker container.

To get started, open the repository in [Codespaces](https://github.com/nf-core/rnaseq/codespaces).

### Testing

Once you have made your changes, run the pipeline with nf-test to test them locally.
For additional information, use the `--verbose` flag to view the Nextflow console log output.

> [!IMPORTANT]
> Unlike most nf-core pipelines, this pipeline does **not** set a default `profile "test"` in `nf-test.config`. This is because the pipeline supports both CPU and GPU test profiles (`test` and `test_gpu`) with different resource limits, and hardcoding one would prevent the other from being used in CI. You must always include the `test` profile explicitly when running tests locally (e.g. `--profile=+test,docker`).

```bash
nf-test test --tag test --profile +docker --verbose
```

If you have added new functionality, ensure you update the test assertions in the `.nf.test` files in the `tests/` directory.
Update the snapshots with the following command:

```bash
nf-test test --tag test --profile +docker --verbose --update-snapshots
```

When you create a pull request with changes, GitHub Actions will run automatic tests.
Pull requests are typically reviewed when these tests are passing.

Two types of tests are typically run:

#### Lint tests

nf-core has a [set of guidelines](https://nf-co.re/docs/specifications/overview) which all pipelines must follow.
To enforce these, run linting with nf-core/tools:

```bash
nf-core pipelines lint <pipeline_directory>
```

If you encounter failures or warnings, follow the linked documentation printed to screen.
For more information about linting tests, see [nf-core/tools API documentation](https://nf-co.re/docs/nf-core-tools/api_reference/latest/pipeline_lint_tests/actions_awsfulltest).

#### Pipeline tests

Each nf-core pipeline should be set up with a minimal set of test data.
GitHub Actions runs the pipeline on this data to ensure it runs through and exits successfully.
If there are any failures then the automated tests fail.
These tests are run with the latest available version of Nextflow and the minimum required version specified in the pipeline code.

### Patch release

> [!WARNING]
> Only in the unlikely event of a release that contains a critical bug.

- [ ] Create a new branch `patch` on your fork based on `upstream/main` or `upstream/master`.
- [ ] Fix the bug and use nf-core/tools to bump the version to the next semantic version, for example, `1.2.3` → `1.2.4`.
- [ ] Open a Pull Request from `patch` directly to `main`/`master` with the changes.

### Pipeline contribution conventions

nf-core semi-standardises how you write code and other contributions to make the nf-core/rnaseq code and processing logic more understandable for new contributors and to ensure quality.

#### Add a new pipeline step

To contribute a new step to the pipeline, follow the general nf-core coding procedure.
Please also refer to the [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines):

- [ ] Define the corresponding [input channel](#channel-naming-schemes) into your new process from the expected previous process channel.
- [ ] Install a module with nf-core/tools, or write a local module (see [default processes resource requirements](#default-processes-resource-requirements)), and add it to the target `<workflow>.nf`.
- [ ] Define the output channel if needed. Mix the version output channel into `ch_versions` and relevant files into `ch_multiqc`.
- [ ] Add new or updated parameters to `nextflow.config` with a [default value](#default-parameter-values).
- [ ] Add new or updated parameters and relevant help text to `nextflow_schema.json` with [nf-core/tools](#default-parameter-values).
- [ ] Add validation for relevant parameters to the pipeline utilisation section of `utils_nfcore_\_pipeline/main.nf` subworkflow.
- [ ] Perform local tests to validate that the new code works as expected.
  - [ ] If applicable, add a new test in the `tests` directory.
- [ ] Update `usage.md`, `output.md`, and `citation.md` as appropriate.
- [ ] [Lint](lint) the code with nf-core/tools.
- [ ] Update any diagrams or pipeline images as necessary.
- [ ] Update MultiQC config `assets/multiqc_config.yml` so relevant suffixes, file name cleanup, and module plots are in the appropriate order.
- [ ] If applicable, create a [MultiQC](https://seqera.io/multiqc/) module.
- [ ] Add a description of the output files and, if relevant, images from the MultiQC report to `docs/output.md`.

To update the minimum required Nextflow version, see the [Nextflow version bumping](#nextflow-version-bumping) section below. For more information about pipeline contributions, see [pipeline-specific contribution guidelines](#pipeline-specific-contribution-guidelines).

#### Channel naming schemes

Use the following naming schemes for channels to make the channel flow easier to understand:

- Initial process channel: `ch_output_from_<process>`
- Intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

#### Default parameter values

Parameters should be initialised and defined with default values within the `params` scope in `nextflow.config`.
They should also be documented in the pipeline JSON schema.

To update `nextflow_schema.json`, run:

```bash
nf-core pipelines schema build
```

The schema builder interface that loads in your browser should automatically update the defaults in the parameter documentation.

#### Default processes resource requirements

If you write a local module, specify a default set of resource requirements for the process.

Sensible defaults for process resource requirements (CPUs, memory, time) should be defined in `conf/base.config`.
Specify these with generic `withLabel:` selectors, so they can be shared across multiple processes and steps of the pipeline.

nf-core provides a set of standard labels that you should follow where possible, as seen in the [nf-core pipeline template](https://github.com/nf-core/tools/blob/main/nf_core/pipeline-template/conf/base.config).
These labels define resource defaults for single-core processes, modules that require a GPU, and different levels of multi-core configurations with increasing memory requirements.

Values assigned within these labels can be dynamically passed to a tool using the the `${task.cpus}` and `${task.memory}` Nextflow variables in the `script:` block of a module (see an example in the [modules repository](https://github.com/nf-core/modules/blob/bd1b6a40f55933d94b8c9ca94ec8c1ea0eaf4b82/modules/nf-core/samtools/bam2fq/main.nf#L30)).

#### Nextflow version bumping

If you use a new feature from core Nextflow, bump the minimum required Nextflow version in the pipeline with:

```bash
nf-core pipelines bump-version --nextflow . <min_nf_version>
```

#### Images and figures guidelines

If you update images or graphics, follow the nf-core [style guidelines](https://nf-co.re/docs/community/brand/workflow-schematics).

## Pipeline specific contribution guidelines

A few conventions that are specific to nf-core/rnaseq and tend to surprise new contributors:

#### Test profiles

The pipeline ships three test profile families: `test` (CPU smoke test, ~15 GB resourceLimits), `test_prokaryotic` (composes on top of `test`, swaps in the bacterial/archaeal samplesheet and `bowtie2_salmon` aligner), and `test_gpu` (~30 GB resourceLimits for GPU CI). `nf-test.config` deliberately does **not** set a default profile - you must always pass one explicitly (e.g. `--profile +test,docker`).

GPU/license-server-bound CI cases are gated on env vars so they don't run by default on contributor laptops:

- `SKIP_GPU=1` skips Parabricks and GPU ribodetector tests
- `SKIP_SENTIEON=1` skips Sentieon STAR tests
- `SKIP_PARABRICKS=1` is a finer-grained subset of `SKIP_GPU`

#### Test data

Eukaryotic test data is pulled from the iGenomes S3 mirror (`pipelines_testdata_base_path = s3://ngi-igenomes/testdata/nf-core/pipelines/rnaseq/3.15/`); prokaryotic test data is pulled from `nf-core/test-datasets` on GitHub (Salmonella SL1344 subset). Both are pinned in `conf/test.config` and `conf/test_prokaryotic.config`; if you need to add a new fixture, prefer extending one of those rather than introducing a new bucket.

#### nf-core modules and subworkflows

Modules and subworkflows under `modules/nf-core/` and `subworkflows/nf-core/` are managed by the nf-core tooling - install or update them with `nf-core modules install` / `nf-core subworkflows update`, do **not** edit them directly. Edits should go through a PR to [nf-core/modules](https://github.com/nf-core/modules) first; the pipeline then picks up the new SHA via a `modules.json` bump. Pipeline-specific code (anything not portable to other pipelines) lives under `modules/local/` and `subworkflows/local/`.

#### Module configs

Per-tool publishDir, ext.args, and ext.prefix settings are split into one file per logical group under `conf/modules/` (e.g. `conf/modules/align_star.config`, `conf/modules/quantify_rsem.config`) and included from `nextflow.config`. When you add a new local module, add or extend the matching file rather than dropping settings into `nextflow.config` directly.

#### Version reporting

Modules emit their versions onto the `versions` channel topic so the calling workflow does not have to thread a `ch_versions` through every process (PR #1689). Modules that still also declare a `path "versions.yml", emit: versions` output do so because they are templated (the `.r`/`.py` template script writes the YAML); those modules populate the topic too and don't need migrating - leave them alone.

#### `--genome` reference catalogues

`--genome <key>` resolves a bundle of reference paths from the `params.genomes` map. The pipeline ships the iGenomes catalogue out of the box, but the same mechanism works with a user-authored catalogue - that is the recommended path for modern reference data, since the iGenomes annotations are stale. If you add new fields to genome map entries, make them optional and gate behaviour on their presence (see e.g. the `star_legacy` flag in `conf/igenomes.config`).

#### Snapshots

Non-deterministic outputs (STAR, Salmon, Kallisto, RSEM, HISAT2 indices; qualimap reports) are snapshotted by file-name-only (`getSnapshot()` filtered) rather than content. Deterministic text outputs are snapshotted by md5. Verbose JSON test output (e.g. helper-function tests) should snapshot `.md5()` of the result rather than inlining the JSON. Don't snapshot timestamps or paths that contain hash directories.

#### `.nftignore`

`tests/.nftignore` (and `tests/.nftignore_rustqc` for the RustQC variant) is a list of glob patterns that nf-test excludes from `${outputDir}` snapshots at the pipeline-test level. It is the right place to drop outputs that are content-stable but not byte-stable (e.g. files that include a timestamp, paths under multiqc/multiqc_data, log files where ordering varies), or that are already covered elsewhere. If a pipeline-level snapshot is fluttering on a file you don't actually need to assert on, add it here rather than rerunning until you get lucky.

#### CHANGELOG

One-line entry per PR under the unreleased section, focused on the _what_, not the implementation history. A `### Software dependencies` table at the end of each release section captures tool version bumps (Old → New).
