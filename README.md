# README for Cancer Research journal: CFL1 and HIF1A DepMap Analysis

This R script generates the figures included in the supplementary data of the manuscript 
"IL-9 induces metastatic migration of malignant T-cells by activating a downstream HIF-1-Cofilin-1 axis in cutaneous T-cell lymphoma."

## Dependencies:
- R (version 3.4.1)
- Required packages: dplyr, tidyr, ggplot2, ggrepel

## Instructions:
1. Set the working directory to the folder containing the data files.
2. Install the required R packages (if not already installed).
3. Run the script in R or RStudio.

## Generated Figures:
- **Figure S3A**: HIF1A CRISPR KO Gene Effect vs Expression.
- **Figure S3B**: CFL1 CRISPR KO Gene Effect vs Expression.
- **Figure S3C**: Density plot for CFL1 Gene Effect.

## Data Files:
- "CFL1 CRISPR.csv": DepMap gene effect data for CFL1: 
Link: https://depmap.org/portal/partials/entity_summary/download?entity_id=4505&dep_enum_name=Chronos_Combined&size_biom_enum_name=expression&color=mutations_prioritized
- "HIF1A CRISPR.csv": DepMap gene effect data for HIF1A.
Link:
https://depmap.org/portal/partials/entity_summary/download?entity_id=11600&dep_enum_name=Chronos_Combined&size_biom_enum_name=expression&color=mutations_prioritized
- "OmicsExpression.csv": Omics expression data. 
Link:
https://depmap.org/portal/partials/entity_summary/download?entity_id=11600&dep_enum_name=Chronos_Combined&size_biom_enum_name=expression&color=mutations_prioritized
- "T cell cancers.csv": Metadata for T cell cancer lines.
See file: T cell cancers.csv

