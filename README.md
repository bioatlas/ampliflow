# ![bioatlas/ampliflow](https://raw.githubusercontent.com/bioatlas/artwork/master/images/bioatlas-logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A518.10.1-brightgreen.svg)](https://www.nextflow.io/)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)

[![Docker](https://img.shields.io/docker/automated/nfcore/ampliseq.svg)](https://hub.docker.com/r/nfcore/ampliseq)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
**bioatlas/ampliflow** is a bioinformatics analysis pipeline used for several kinds of rRNA amplicon sequencing data.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)), trims primer sequences from the reads ([Cutadapt](https://journal.embnet.org/index.php/embnetjournal/article/view/200)), performs denoising and generates amplicon sequencing variants (ASV, [DADA2](https://www.nature.com/articles/nmeth.3869)), classifies features against prokaryotic and eucariotic databases upon user's demand, including [SILVA](https://www.arb-silva.de/) [v132](https://www.arb-silva.de/documentation/release-132/) [GTDB](https://gtdb.ecogenomic.org/) [r86](https://zenodo.org/record/2541239), [UNITE](https://unite.ut.ee/) [general release](https://unite.ut.ee/repository.php) [PR2](https://github.com/pr2database/pr2database) [v4.12.0](https://github.com/pr2database/pr2database/releases/tag/v4.12.0), produces relative feature/taxa count tables.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible. It also revolves on DADA2 R scripts deposited at [eemis-dada2](https://github.com/erikrikarddaniel/eemisdada2/) (Author: [Daniel Lundin](https://github.com/erikrikarddaniel))

### Documentation
The bioatlas/ampliflow pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Local installation](https://nf-co.re/usage/installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](https://nf-co.re/usage/troubleshooting)

### Credits
This pipeline has been developed by [Diego Brambilla](https://github.com/DiegoBrambilla) (developer of Nextflow code) and [Daniel Lundin](https://github.com/erikrikarddaniel) (supervisor and developer of R scripts) at the [Marine Microbiology research group](https://lnu.se/en/research/searchresearch/marin-microbiology/), part of [Linnaeus University Centre for Ecology and Evolution in Micobial model Systems](https://lnu.se/en/research/searchresearch/linnaeus-university-centre-for-ecology-and-evolution-in-microbial-model-systems/). 
These scripts were originally forked from [nf-core/ampliseq](https://github.com/nf-core/ampliseq), which has been written for use at the [Quantitative Biology Center (QBiC)](http://www.qbic.life) and [Microbial Ecology, Center for Applied Geosciences](http://www.uni-tuebingen.de/de/104325), part of Eberhard Karls Universität Tübingen (Germany) by Daniel Straub ([@d4straub](https://github.com/d4straub)) and Alexander Peltzer ([@apeltzer](https://github.com/apeltzer)).
