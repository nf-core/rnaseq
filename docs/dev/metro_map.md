# Metro map

The pipeline overview metro map is generated from `assets/metro_map.mmd` using [nf-metro](https://github.com/pinin4fjords/nf-metro). If you add or rename pipeline steps, update the `.mmd` source and regenerate the images:

```bash
pip install 'nf-metro>=1.1.0' cairosvg

# Theme, diamond style, and logo are set via %%metro directives in the .mmd
# source, so the render command only needs output-serialization flags.

# Static SVG + PNG
# --no-chrome-css bakes concrete colors so cairosvg can rasterize the file;
# without it cairosvg aborts on the CSS custom properties used for theming.
# --no-chrome-css also drops the light/dark pair entirely (bakes to one
# resolved mode), so pin one explicitly rather than relying on the default.
nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-rnaseq_metro_map_grey.svg \
  --no-chrome-css --mode light

python -c "import cairosvg; cairosvg.svg2png(
    url='docs/images/nf-core-rnaseq_metro_map_grey.svg',
    write_to='docs/images/nf-core-rnaseq_metro_map_grey.png', output_width=2265)"

# Animated SVG (used in README) - no --no-chrome-css here: this file is
# viewed directly in a browser, never rasterized, so it keeps light-dark()
# and follows the viewer's (or GitHub's) light/dark setting.
nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-rnaseq_metro_map_grey_animated.svg \
  --animate

# Copy static PNG to docs subdir
cp docs/images/nf-core-rnaseq_metro_map_grey.png \
  docs/usage/differential_expression_analysis/img/

# Ensure trailing newlines on SVGs (required by pre-commit)
for f in docs/images/nf-core-rnaseq_metro_map_grey.svg \
         docs/images/nf-core-rnaseq_metro_map_grey_animated.svg; do
  sed -i '' -e '$a\' "$f"
done
```

## Live progress overlay

The `.mmd` file includes `%%metro process:` directives that tie each station to its
Nextflow fully-qualified process name. These are embedded in the SVG manifest at
render time and enable `nf-metro serve` to light up stations in real time as the
pipeline runs:

```bash
pip install 'nf-metro>=1.1.0'

# Serve the map and start the pipeline (one-liner)
nf-metro serve assets/metro_map.mmd --open --shutdown-after-complete -- \
    nextflow run nf-core/rnaseq -profile test,docker --outdir results

# Or serve from the committed SVG (no source needed)
nf-metro serve docs/images/nf-core-rnaseq_metro_map_grey.svg --open --shutdown-after-complete -- \
    nextflow run nf-core/rnaseq -profile test,docker --outdir results
```

To verify all stations are correctly wired after editing the map:

```bash
nf-metro check-mapping assets/metro_map.mmd \
    --run-log path/to/nextflow/.nextflow.log
```
