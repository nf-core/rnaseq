# Metro map

The pipeline overview metro map is generated from `assets/metro_map.mmd` using [nf-metro](https://github.com/pinin4fjords/nf-metro). If you add or rename pipeline steps, update the `.mmd` source and regenerate the images:

```bash
pip install 'nf-metro>=0.5.4' cairosvg

# Static SVG + PNG
nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-rnaseq_metro_map_grey.svg \
  --theme light --x-spacing 60 --y-spacing 40 \
  --no-straight-diamonds \
  --logo docs/images/nf-core-rnaseq_logo_light.png

python -c "import cairosvg; cairosvg.svg2png(
    url='docs/images/nf-core-rnaseq_metro_map_grey.svg',
    write_to='docs/images/nf-core-rnaseq_metro_map_grey.png', output_width=2265)"

# Animated SVG (used in README)
nf-metro render assets/metro_map.mmd \
  -o docs/images/nf-core-rnaseq_metro_map_grey_animated.svg \
  --theme light --x-spacing 60 --y-spacing 40 --animate \
  --no-straight-diamonds \
  --logo docs/images/nf-core-rnaseq_logo_light.png

# Copy static PNG to docs subdir
cp docs/images/nf-core-rnaseq_metro_map_grey.png \
  docs/usage/differential_expression_analysis/img/

# Ensure trailing newlines on SVGs (required by pre-commit)
for f in docs/images/nf-core-rnaseq_metro_map_grey.svg \
         docs/images/nf-core-rnaseq_metro_map_grey_animated.svg; do
  sed -i '' -e '$a\' "$f"
done
```
