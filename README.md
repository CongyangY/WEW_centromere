# Wild emmer wheat centromere analyses

Custom R scripts used for the wild emmer wheat (WEW) centromere, CRW, CentT566, structural-variation and comparative analyses.

An English summary of the BinBash platform WGS workflow is provided in `BinBash_WGS_workflow.md`.

The sanitized genome assembly and annotation workflow is provided in `Genome_assembly_annotation_workflow.md`.

## Scope

This repository contains analysis scripts only. Large genomic inputs, source tables and generated figures are intentionally excluded. The scripts were preserved from the working analysis directory and grouped by analysis topic; they may require project-specific input files and software packages.

## Directory layout

```text
R/
  CENH3.R                         CENH3 signal and centromere analyses
  Cereba_new.R                    CRW insertion-time and phylogenetic analyses
  Response.R                      revised-response figure analyses
  WEW_durum_centromere.R          SV and centromere-position analyses
  cent566.R                       CentT566 abundance and similarity analyses
  cent_TE.R                       centromeric TE analyses
  centsSimpleRepeat.R             centromeric simple-repeat summaries
  flLTR_burst.R                   full-length LTR insertion-time analyses
  polyploidy.R                    polyploid comparison analyses
  run_fisher_for_clusters.R       local-cluster Fisher tests
  solo_LTR.R                      solo-LTR analyses
  supplementary.R                 supplementary-figure analyses
```

## Reproducibility

The original scripts use relative paths such as `../data/` and `../figure/`. To reproduce a result, place the required project data and figure directories at the expected relative locations, or adapt the input/output paths in a local copy. The repository does not include raw assemblies, FASTA/BAM/VCF files, large alignment files, or generated PDFs.

All statistical analyses were conducted in R v4.3.2, as specified in the manuscript Methods. The Methods do not specify exact versions for ggplot2 or ape; package versions required by individual scripts should therefore be checked against the execution environment before reproduction.

## Data and code availability

Sequencing data generated in the study, including PacBio HiFi, ultra-long ONT, Illumina, RNA-seq and histone ChIP-seq data, are deposited in the National Genomics Data Center under BioProject accession `PRJCA051525`. The genome assembly is available under accession `GWHHHUI00000000.1`.

The genome assembly, annotation files, gene-by-sample TPM matrix and SNP VCF files are archived at Zenodo: [10.5281/zenodo.19882165](https://doi.org/10.5281/zenodo.19882165). This repository contains the custom analysis scripts and is archived separately at Zenodo: [10.5281/zenodo.22177298](https://doi.org/10.5281/zenodo.22177298).

## License

No open-source license is specified in the manuscript Methods for this code archive. Reuse, redistribution or modification should follow the authors' terms until a license is added by the authors.
