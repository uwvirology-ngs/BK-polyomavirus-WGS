# BK-polyomavirus-WGS

[![UW DLMP](https://img.shields.io/badge/dev-UW_Medicine-32006e)](https://dlmp.uw.edu/)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.6-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![GitHub Actions CI Status](https://github.com/uwvirology-ngs/BK-polyomavirus-WGS/actions/workflows/ci-tests.yml/badge.svg)](https://github.com/uwvirology-ngs/BK-polyomavirus-WGS/actions/workflows/ci-tests.yml)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Description

Analyze BK Polyomavirus sequencing libraries prepared with the Twist Unique Molecular Identifier (UMI) Adapter System.

Here we produce alignments, consensus sequences, tabulated variants, and UMI-aware duplex consensus sequences given FASTQ files created by Illumina short-read sequencing experiments on BK Polyomavirus samples. This implementation follows and builds upon the guideline provided by [Twist Bioscience](https://www.twistbioscience.com/) for processing sequencing data prepared with their [UMI adapter system](https://www.twistbioscience.com/products/ngs/library-preparation/twist-umi-adapter-system). 

## Limitations

This pipeline is actively being developed and is midway though validation for clinical research purposes. We currently only support analysis of BK Polyomavirus genomes with a particular database of reference sequences, however support for custom databases and generalization to any pathogen is intended for future release.

## Workflow (UMI-aware steps)

<img src="./assets/img/bkv_workflow.drawio.png">

## Requirements

This pipeline is tested with [nf-test](https://www.nf-test.com/) on Ubuntu 26.04 LTS with the latest release of Nextflow and version 26.04.6. 

Install [`Nextflow`](https://www.nextflow.io/docs/latest/install.html)

Docker v29.4.2+ is supported for containerization. 

Install [`Docker`](https://docs.docker.com/engine/install/)

## Usage

### Example Run

To run an example dataset, clone this repository and execute the following command in the project root directory: 

```bash
nextflow run main.nf \
    -params-file assets/example_run/params.json
```

### Integration with Amazon Web Services

A profile is included for execution with [awsbatch](https://docs.seqera.io/nextflow/executor#aws-batch). This profile is tailored to users within the UW's Divison of Infectious Disease Diagnostics, but can be used freely by tuning the resource limits to your HPC specs. 

For reference, here is a minimal example command for execution with awsbatch:

```bash
nextflow run uwvirology-ngs/bk-polyomavirus-wgs -r main -latest \
    --input your_samplesheet.csv \
    --output your_output_directory \
    --sampleTo 2000000 \
    -profile awsbatch \
    -c your_nextflow_aws.config
```

## CLI Options

### Required Parameters
| Parameter | About | Example |
|---------|---------|---------|
| `--input` | samplesheet | /assets/samplesheet.csv |
| `--output` | output directory | results |

### Optional Parameters
| Parameter | About | Default |
|---------|---------|---------|
| `--sampleTo` | maximal reads per FASTQ after downsampling | 1000000 |

### Optional Parameters - Reference Selection
| Parameter | About | Default |
|---------|---------|---------|
| `--db` | multifasta database of reference genomes | /assets/bkv_multi.fa |
| `--ref_min_depth` | minimum sequencing depth | 3 |
| `--ref_min_cov` | minimum coverage | 30 |

### Optional Parameters - Variant Calling
| Parameter | About | Default |
|---------|---------|---------|
| `--ivar_variants_t` | minimum frequency threshold to call variants | 0.01 |
| `--ivar_variants_q` | minimum quality score threshold to count base | 20 |
| `--ivar_variants_m` | minimum read depth to call variants | 10 |

## Acknowledgements

The reference selection and consensus sequence assembly workflows were adapted from [revica-strm](https://github.com/epiliper/revica-strm) by the [Greninger Lab](https://github.com/greninger-lab). Variant calling was adapted from [nf_mpxv_f13l](https://github.com/greninger-lab/nf_mpxv_f13l) by the same.

## Contact

For feature suggestions and bug reports, please file an issue on the project [GitHub](https://github.com/aidantshea/nf_refvar-strict/issues).

For other inquiries, please reach out to aidants@uw.edu. 