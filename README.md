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

The main observed R environment used during preparation was R 4.3.3 with ggplot2 3.5.0 and ape 5.8.1. Package versions required by individual scripts should be checked against their `library()` calls before execution.

## Data and code availability

The complete raw datasets and intermediate results are available from the authors on reasonable request or through the project data archive. This code repository is intended to document the custom analysis code used for the study.

## License

Add the license selected by the authors before public release.
