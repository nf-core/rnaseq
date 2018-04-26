# nfcore/rnaseq: C3SE (Hebbe) Configuration

This pipeline has been successfully used on the [Hebbe cluster](http://www.c3se.chalmers.se/index.php/Hebbe) in Gothenburg, though it hasn't had as much testing as some of the other profiles.

To use, run the pipeline with `-profile hebbe --project [project-id]`. This will launch the [hebbe config](../../conf/hebbe.config) which has been pre-configured with a setup suitable for the Hebbe cluster. It will download a singularity image with all of the required software - see [the c3se Singularity documentation](http://www.c3se.chalmers.se/index.php/Singularity) for more details.

If running regularly, we recommend creating a config file with paths to your reference genome indices (see [reference-genomes.md](reference-genomes.md) for instructions).
