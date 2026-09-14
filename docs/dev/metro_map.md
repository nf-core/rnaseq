# Metro map

The pipeline overview metro map is generated from `assets/metro_map.mmd` using [nf-metro](https://github.com/pinin4fjords/nf-metro). If you add or rename pipeline steps, update the `.mmd` source and regenerate the images:

```bash
pip install 'nf-metro>=2.0.0'

# Static SVG + dark-mode PNG
nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-rnaseq_metro_map.svg \
  -o docs/images/nf-core-rnaseq_metro_map_dark.png \
  -o docs/usage/differential_expression_analysis/img/nf-core-rnaseq_metro_map.svg

# Static light-mode PNG
nf-metro render assets/metro_map.mmd --mode light \
  -o docs/images/nf-core-rnaseq_metro_map_light.png

# Animated SVG (used in manifest + README)
nf-metro render assets/metro_map.mmd --animate \
  -o docs/images/nf-core-rnaseq_metro_map_animated.svg
```

## Live progress overlay

The `.mmd` file includes `%%metro process:` directives that tie each station to its
Nextflow fully-qualified process name. These are embedded in the SVG manifest at
render time and enable `nf-metro serve` to light up stations in real time as the
pipeline runs:

```bash
pip install 'nf-metro>=2.0.0'

# Serve the map and start the pipeline (one-liner)
nf-metro serve assets/metro_map.mmd --open --shutdown-after-complete -- \
    nextflow run nf-core/rnaseq -profile test,docker --outdir results

# Or serve from the committed SVG (no source needed)
nf-metro serve docs/images/nf-core-rnaseq_metro_map.svg --open --shutdown-after-complete -- \
    nextflow run nf-core/rnaseq -profile test,docker --outdir results
```

To verify all stations are correctly wired after editing the map, export a
Nextflow DAG from a run and check the mapping against it:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir results -with-dag dag.mmd
nf-metro check-mapping assets/metro_map.mmd --dag dag.mmd
```
