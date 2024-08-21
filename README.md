# Preprocess scRNAseq data

This workflow is specific to the CRUK CI computational infrastructure.

The workflow downloads the single cell data as fastq from the Genomics server using
Clarity tools, it then renames the fastq files to conform with CellRangers
expected input file name format. The workflow then runs CellRanger count on the
fastq files to generate the gene expression matrix.

Currently it is only configured to run CellRanger count; the intention is to
add functionality for other assays/chemistries as needed.

## Provide the SLX id

The SLX id should be specified in a parameters yaml file or on the command line:

* `slx_id`: The SLX ID of the sequencing run

All of the workflow outputs will be published into a directory named after the
SLX id.

## CellRanger reference

To specify the reference you can just provide the species in the parameters:

* `species` - "mus_musculus" or "homo_sapiens"

This will cause the workflow to download the reference data for the
specified species from 10X:

*homo_sapiens*: https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz  
*mus_musculus*: https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz  

The reference data will be published into the directory:

* **_${launchDir}_/references**

Alternatively, if you wish to use an existing reference specify:

* `reference_dir` - The path to the CellRanger reference_dir

## CellRanger software

The workflow has a singularity container that contains the CellRanger software
version 8.0.1 and an R installation. This is currenly located at:

* /home/bioinformatics/software/containers/cruk_ci_preprocess_scSeq-8.0.1.sif

The intention is to keep the container up to date as new versions of CellRanger
are released.

If you wish to use a different version of CellRanger specify:

* `cellranger_dir` - The path to the CellRanger software directory

## Outputs

The workflow will generate the following outputs in the SLX directory:

* **fastq** - the raw fastq files and accompanying files as downloaded
              from the genomics server using Clarity tools
* **reports/**:
    * ./*\<Barcode\>*.web_summary.html - The CellRanger web summary reports for
                                     each sample, with the barcode added to file name
    * ./collated_metrics_file.csv - The summary metrics for all samples collated
                                    into a single file
    * ./summary_metrics.pdf - Bar plots of read depth and cell count per sample
* **_\<Barcode\>_** - One directory of CellRanger count output for each
                  sample, named according to the sample barcode

## Running the pipeline

The pipeline can be run using the following command:

```bash
nextflow run crukci-bioinformatics/nf_preprocess_scSeq_CRUKCI \
    --slxid {slxID} \
    --species {species} \
    -profile epyc
```
