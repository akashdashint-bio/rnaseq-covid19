Differential gene expression analysis of COVID-19 vs healthy whole blood using DESeq2 and R
RNA-seq Analysis: COVID-19 vs Healthy Whole Blood
Differential gene expression analysis of whole blood transcriptomics from COVID-19 patients and healthy controls, using publicly available data from GSE152641.

**Background**
COVID-19 (caused by SARS-CoV-2) triggers a broad and complex transcriptional response in blood cells, involving immune activation, antiviral defence, and neurological disruption. This project performs an end-to-end RNA-seq analysis to identify differentially expressed genes (DEGs) and enriched biological pathways in COVID-19 patients compared to healthy controls.
The dataset comprises whole blood samples from 62 COVID-19 patients and 24 healthy controls, profiling 19,936 genes after quality filtering.

**Key Findings**
9,168 significant DEGs identified (padj < 0.05, |log2FC| > 1)
6,294 upregulated in COVID-19
2,900 downregulated in COVID-19


Top upregulated gene: IFI27 — a classical interferon-stimulated gene, consistent with the known COVID-19 interferon response
APOBEC3A/APOBEC3A_B strongly upregulated — antiviral RNA-editing enzymes directly targeting the virus
Cell cycle disruption: RRM2, CDC6, CLSPN, MCM10, POLQ — COVID-19 hijacks DNA replication machinery, with implications for cancer biology (several are oncogenes)
Pathway enrichment reveals anosmia signature: the most significantly enriched GO terms were sensory perception of smell pathways — directly reflecting COVID-19's well-known neurological symptom of smell loss at a transcriptomic level

**Methods**
StepToolPurposeData retrievalGEOqueryDownload count matrix and metadata from NCBI GEODE analysisDESeq2Normalisation, dispersion estimation, negative binomial testingVisualisationggplot2, ggrepel, pheatmapVolcano plot, PCA, heatmapGene annotationorg.Hs.eg.db, AnnotationDbiEntrez ID to gene symbol conversionPathway enrichmentclusterProfilerGO biological process over-representation analysis

**Relevance to Research Fields**
Infectious Disease: Characterises the whole-blood transcriptional response to SARS-CoV-2, identifying key antiviral (APOBEC3A) and interferon (IFI27) mediators.
Neuroscience: Enrichment of olfactory sensory perception pathways in peripheral blood provides a systemic transcriptomic signature of COVID-19-associated anosmia — a clinically important neurological symptom.
Cancer Biology: Multiple top DEGs (RRM2, CDC6, POLQ, CLSPN, MCM10) are established oncogenes involved in DNA replication and repair. COVID-19 appears to dysregulate the same pathways exploited by cancer cells.

_> Reproducing This Analysis_

**Requirements**
R ≥ 4.3.0
RStudio

Install dependencies
rif (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "GEOquery"))
install.packages(c("ggplot2", "ggrepel", "pheatmap", "R.utils"))
Download data
rlibrary(GEOquery)
gse <- getGEO("GSE152641", GSEMatrix = TRUE)
getGEOSuppFiles("GSE152641")
Run analysis
Open scripts/analysis.R in RStudio and run section by section, or source the entire script:
rsource("scripts/analysis.R")

**Project Structure**
rnaseq-covid19/
├── data/
│   ├── GSE152641_counts.csv       # raw count matrix
│   └── metadata.csv               # sample condition labels
├── scripts/
│   └── analysis.R                 # full analysis pipeline
├── results/
│   ├── figures/
│   │   ├── pca_plot.png
│   │   ├── volcano_plot.png
│   │   └── pathway_enrichment.png
│   └── tables/
│       ├── DESeq2_results.csv
│       └── GO_enrichment.csv
├── .gitignore
├── LICENSE
└── README.md

**Dataset**
GEO Accession: GSE152641
Publication: Thair et al. (2020) Transcriptomic similarities and differences in host response between SARS-CoV-2 and other viral infections. iScience.

**Author**
Akash Dash | @AkashDash-bio
Independent bioinformatics project — RNA-seq differential expression analysis
