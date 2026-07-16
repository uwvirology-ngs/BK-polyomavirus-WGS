# BK-polyomavirus-WGS

[![UW DLMP](https://img.shields.io/badge/dev-UW_Medicine-32006e)](https://dlmp.uw.edu/)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.6-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![GitHub Actions CI Status](https://github.com/uwvirology-ngs/BK-polyomavirus-WGS/actions/workflows/ci-tests.yml/badge.svg)](https://github.com/uwvirology-ngs/BK-polyomavirus-WGS/actions/workflows/ci-tests.yml)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Description

Analyze sequencing libraries prepared with the Twist Unique Molecular Identifier (UMI) Adapter System.

Here we implement the data analysis guideline provided by [Twist Bioscience](https://www.twistbioscience.com/) for processing sequencing data prepared with their [UMI adapter system](https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system) and incorporate automated reference genome selection. 

## Workflow

<img src="./assets/img/bkv_workflow.drawio.png">

## Requirements

This pipeline is tested with [nf-test](https://www.nf-test.com/) on Ubuntu 26.04 LTS with Nextflow v26.04.6+ and Docker v29.4.2. 

Install [`Nextflow`](https://www.nextflow.io/docs/latest/install.html)

Install [`Docker`](https://docs.docker.com/engine/install/)

## Usage

### Reference Selection

By default, a reference genome is selected from a database of BK Polyomavirus (BKPyV) genomes, however any multifasta may be provided on the command line as a database. Alternatively, automated selection can be entirely circumvented by specifying a reference genome directly. 

These default settings stem from this pipeline's first application: analyzing BKPyV genomes. 

### Example Run

An example dataset is provided for users to try out the pipeline locally. To run the example, first clone this repository, then execute the following command in the project root directory: 

```bash
nextflow run main.nf \
    -params-file assets/example_run/params.json
```

### Integration with Amazon Web Services

A profile is included for pipeline execution with [awsbatch](https://docs.seqera.io/nextflow/executor#aws-batch). This profile is intended for users at the University of Washington but can be used freely, though we suggest tuning the resource limits to your HPC specs. [Fusion file system](https://docs.seqera.io/fusion) is enabled here. Consequently, an API token from [Seqera Platform](https://cloud.seqera.io) is required to run with this profile. 

Once generated, this token can be defined in your `.bashrc` file for convenience:

```bash
export TOWER_ACCESS_TOKEN="your_token_here"
```

For reference, here is a minimal example command for execution with awsbatch:

```bash
nextflow run uwvirology-ngs/bk-polyomavirus-wgs -r main -latest \
    --input your_samplesheet.csv \
    --output your_output_directory \
    -profile awsbatch \
    -c your_nextflow_aws.config
```

## CLI Options

### Required Parameters
| Parameter | About | Example |
|---------|---------|---------|
| `--input` | samplesheet | /assets/samplesheet.csv |
| `--output` | output directory | results |

### Optional Parameters - Twist UMI Workflow
| Parameter | About | Default |
|---------|---------|---------|
| `--ref` | reference genome | null |

### Optional Parameters - Reference Selection
| Parameter | About | Default |
|---------|---------|---------|
| `--db` | multifasta database of reference genomes | /assets/bkv_multi.fa |
| `--ref_min_depth` | minimum sequencing depth | 3 |
| `--ref_min_cov` | minimum coverage | 30 |

## Acknowledgements

Reference selection logic incorporated from [revica-strm](https://github.com/epiliper/revica-strm) by the [Greninger Lab](https://github.com/greninger-lab).

## Contact

For feature suggestions and bug reports, please file an issue on the project [GitHub](https://github.com/aidantshea/nf_refvar-strict/issues).

For other inquiries, please reach out to aidants@uw.edu. 