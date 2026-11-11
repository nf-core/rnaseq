# Metro map

The pipeline overview metro map is generated from `assets/metro_map.mmd` using [nf-metro](https://github.com/seqeralabs/nf-metro). See the [nf-core workflow schematics guide](https://nf-co.re/docs/community/brand/workflow-schematics) for `%%metro` directive syntax and general usage.

If you add or rename pipeline steps, update the `.mmd` source and regenerate all three copies of the image:

```bash
nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-rnaseq_metro_map.svg \
  -o docs/images/nf-core-rnaseq_metro_map_dark.png \
  -o docs/usage/differential_expression_analysis/img/nf-core-rnaseq_metro_map.svg

nf-metro render assets/metro_map.mmd --mode light \
  -o docs/images/nf-core-rnaseq_metro_map_light.png

nf-metro render assets/metro_map.mmd --animate \
  -o docs/images/nf-core-rnaseq_metro_map_animated.svg
```

## Live progress overlay

`%%metro process:` directives tie each station to its Nextflow process name, so `nf-metro serve` can light stations up as a run progresses:

```bash
nf-metro serve assets/metro_map.mmd --open --shutdown-after-complete -- \
    nextflow run nf-core/rnaseq -profile test,docker --outdir results
```

After editing the map, verify the mapping against a real run:

```bash
nextflow run nf-core/rnaseq -profile test,docker --outdir results -with-dag dag.mmd
nf-metro check-mapping assets/metro_map.mmd --dag dag.mmd
```
