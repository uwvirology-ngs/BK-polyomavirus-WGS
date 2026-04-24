# BK-polyomavirus-WGS

## Description

Whole genome sequencing for BK polyomavirus

## Requirements

This pipeline is developed and tested on Ubuntu 24.04 LTS with Nextflow version 25.10.4+ and Docker v29.3.1+. Compatibility is expected with releases greater than or equal in version number.

Install [`Nextflow`](https://www.nextflow.io/docs/latest/install.html)

Install [`Docker`](https://docs.docker.com/engine/install/)

## Usage

### Run Test:

```bash
nextflow run main.nf \
    -params-file test.json
```

## CLI Options

### Required Parameters
| Parameter | About | Example |
|---------|---------|---------|
| `--samplesheet` | samplesheet | /assets/samplesheet.csv |

### Optional Parameters
| Parameter | About | Default |
|---------|---------|---------|
| `--db` | multifasta with varied BK Polyomavirus genotypes | /assets/ref.fa |
| `--ref_min_depth` | minimum sequencing depth | 3 |
| `--ref_min_cov` | minimum coverage | 30 |

## Acknowledgements

Reference selection logic incorporated from [revica-strm](https://github.com/epiliper/revica-strm) by the [Greninger Lab](https://github.com/greninger-lab).

## Contact

For feature suggestions and bug reports, please file an issue on the project [GitHub](https://github.com/aidantshea/nf_refvar-strict/issues).

For other inquiries, please reach out to aidants@uw.edu. 