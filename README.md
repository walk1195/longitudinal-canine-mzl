## Longitudinal single-cell and spatial transcriptomics reveals intratumor heterogeneity, therapeutic response, and comparative value of canine marginal zone lymphoma.

Walker GE, Macchietto M, Reid K, Burt LE, Winter A, Pracht S, Buettner M, Penza V, Yung C, Dicovitsky R, Kuzmik A, Vallera DA, Demos-Davies K, Seelig DM, Borgatti A, Henson M, Feiock C, Modiano JF, Naik S, Sarver AL, Treeful AE.

**Abstract:** Diffuse B-cell lymphomas are the most prevalent canine hematologic malignancies, and their clinical presentation resembles that of human B-cell non-Hodgkin lymphomas (NHLs). However, a lack of a clear subtyping framework perpetuates imprecise treatment approaches that fail to address mechanisms of therapy resistance. Here, we employed deep phenotyping via serial sampling and longitudinal transcriptomics to evaluate intratumoral composition and therapy response in a 6-year-old neutered male goldendoodle with stage 5A marginal zone lymphoma (MZL). Following sequential treatment consisting of a single IV dose of oncolytic vesicular stomatitis virotherapy (VSV) and standard CHOP chemotherapy 30 days later, we performed scRNAseq on seven serial lymph node biopsies and spatial sequencing on a resection taken 10 days post-VSV treatment. B cells (>90%) comprised three transcriptionally distinct and temporally stable subpopulations, with specific copy number profiles and robust spatial organization despite disrupted lymph node architecture. Analysis across subpopulations revealed downregulation of the quiescence program including KLF2, suggesting similarity to KLF2-deficient human MZL. While VSV successfully reached target B cells, downstream transcriptional responses were predominantly localized to the T-cell compartment. T cells (≤6%) showed increased proliferation and upregulation of cytotoxicity markers post-VSV administration, enriching for leukocyte-mediated cytotoxicity pathways consistent with an anti-viral immune response. Ultimately, these results provide an in vivo proof-of-concept for the treatment paradigm in canine lymphoma, while emphasizing the need for improved subtyping frameworks and multidimensional treatment strategies targeting intratumoral heterogeneity.

---

This repository contains all code written for the analysis of ORv-01 longitudinal scRNAseq and spatial transcriptomics data presented in the manuscript "_Longitudinal single-cell and spatial transcriptomics reveals intratumor heterogeneity, therapeutic response, and comparative value of canine marginal zone lymphoma_".

### Environment
Analysis was performed using `R/4.4.0-openblas-rocky8` and `python/3.12.12`. 

Renv.lock and .yml files to configure single-cell and spatial analyses environments will be provided soon. 

### Usage
The scripts available in this repository are organized into the following subfolders:
- **single-cell/**
- **spatial/**

for respective analyses. The majority are written to be ran interactively, unless otherwise indicated (see script header). Note that all direct file paths have been removed.

### Data Availability
Raw data generated and analyzed in this study will be made publicly available in NCBI GEO upon manuscript acceptance.

---
### Contacts
Any other relevant data and code relating to this study is available from the authors upon request. For any inquiries, please raise an issue in GitHub or email Grace Walker (walk1195@umn.edu, First Author) or Amy Treeful (tree0002@umn.edu, Last Author).
