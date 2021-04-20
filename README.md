# scTPA local application
scTPA is the local tool for scTPA (http://sctpa.bio-data.cn/sctpa) single-cell transcriptome analysis and annotation based on biological pathway activation in human and mice. We collected a large number of biological pathways with different functional and taxonomic classifications, which facilitates the identification of key pathway signatures for cell type annotation and interpretation.

### What can scTPA do
* Calculating pathway activity score of single cell
* Dimension reduction
* Clustering of cell population by different methods
* Identifying significantly activated pathways of cell clusterings
* Comparison analysis of the associated gene expression profiles of pathways

# Installation
### **Step 1. Download scTPA**
scTPA local application can be download directly by
```
wget http://sctpa.bio-data.cn:8888/sctpa/resources/scTPA_local-v7.zip
unzip scTPA_local-v7.zip
cd scTPA_local-v7
```
### **Step 2. Install dependent R packages**
To install this packages, start "R" and enter:
```
install.packages('p2data', repos='https://kharchenkolab.github.io/drat/', type='source')
install.packages("devtools")
install.packages("rJava")
library(rJava)
library(devtools)
library(usethis)
devtools::install_local("dependencies/pagoda2-v0.1.1-master.zip")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("parallel")
BiocManager::install("optparse")
BiocManager::install("scImpute")
BiocManager::install("scran")
BiocManager::install("scater")
BiocManager::install("dplyr")
BiocManager::install("Seurat")
BiocManager::install("cluster")
BiocManager::install("fpc")
BiocManager::install("SIMLR")
BiocManager::install("this.path")
BiocManager::install("scales")
BiocManager::install("ggplot2")
BiocManager::install("cowplot")
BiocManager::install("pheatmap")
BiocManager::install("AUCell")
BiocManager::install("Cairo")
BiocManager::install("scde")
```
### **Step 3. Install dependent Python packages**
```
pip install pandas==0.25.3
pip install numpy
pip install json
pip install clustergrammer
pip install seaborn
pip install multiprocessing
```
# Usage
### Example

```
Rscript src/scTPA.R -f example/e1_scRNA_UMI_count.csv --cellType example/e1_cell_type.csv --data_type count --species homo -o test/test_output_e1
```
```
Rscript src/scTPA.R -f example/e2_melanoma.csv --cellType NULL --species homo --data_type TPM -o test/test_output_e2

```
Once the program has run successfully, a series of results files and folders will appear in the results folder.
The results can be found at the directory **test/test_output_e1** or **test/test_output_e2**,  The file organization is as follows: 
```
+--app.exe          //double click to visualize the results within Windows system.
+--data
   +--app.js
   +--bin
   +--node_modules
   +--content       //Result files and pictures where they are stored.
```

### Help Information
```
Rscript src/scTPA.R -h
```

##### **Required parameters:**
```
    -f FILE, --file=FILE
       Gene expression profile, genes X cells. The processed gene expression profile can be generated using different platforms, such as 10X genomics and Smart-seq. The values in this profile should be non-negative, and this file can be uploaded depending on data types of UMI count, read count, RPKM, FPKM, CPM or TPM. [default= NULL]
    --data_type=FILE
        Data type of gene expression profile，Available options are 'TPM' or 'count'. 'count' indicate that the expression profile is non-negative UMI or read count. 'TPM' indicate that the expression profile is normalized FPKM, RPKM, CPM or TPM. [default= TPM]
    --species=SPECIES
        "Species. Available options are 'homo' or 'mus'. [default= homo]
    --cellType=CELLTYPE
        Optional. Cell type file. First column is cell name (same as the colnames of gene expression profile), second column is cell type. No header names. [default= NULL]
```

##### **Optional parameters:**
```

    --normalize=NORMALIZE_METHOD
        Methods used for normalization. Available options are 'none', 'log', 'CLR', 'RC' 'sctrancform' or 'scran'. 'log', 'CLR' 'RC' and 'sctrancform': The normalization methods from Seurat R package. 'scran': The normalization strategy for scRNA-seq is implemented based on the deconvolutional size factor using the scran R package. "log", "CLR", "RC" or "scran"[default= none]
    --min_cells=MIN_CELLS
        Genes must be detected within a minimum number of cells. Used for filtering genes. [default= 3]
    --min_features=MIN_FEATURES
        Cells must have at least the minimum number of genes. Used for filtering cells. [default= 200]
    --imputation=IMPUTATION
        Imputation method. Available options are 'scImpute' or 'none'. 'scImpute': impute scRNA-seq profile using scImpute R package. [default= none]
    --pathway_database=PATHWAY_DATABASE
        Pathway database. Avalible database are avalible on https://github.com/sulab-wmu/scTPA [default= kegg]
    --user_pathway=USER_PATHWAY
        Optional. User defined pathway file in gmt format. [default = NULL]
    --pas_method=PAS_METHOD
        PAS (pathway activation signatures) transformation method. Available options are 'pagoda2', 'Vision', 'AUCell', 'gsva', 'ssgsea', 'zscore' or 'plage'. [default= ssgsea]
    --para_size=PARA_SIZE
        Number of kernels used for parallel computation. [default= 4]
    --cluster_method=CLUSTER_METHOD
        Clustering method. Available options are 'seurat', 'hclust', 'simlr', 'kmedoids', 'kmeans' or 'dbscan'. [default= seurat]
    --seurat_dims=SEURAT_DIMS
        Dimensions of PCA used in Seurat FindNeighbors method. [default= 8]
    --seurat_resolution=SEURAT_RESOLUTION
        Resolution used in Seurat FindClusters method. [default= 0.5]
    --k_cluster=K_CLUSTER
        Number of clusters. Used for clustering methods except Seurat and dbscan. [default= 5]
    --min_pts=MIN_PTS
        Number of nearest neighbors used in dbscan clustering. [default= 3]
    --dims=DIMS
        Number of PCA dimensions used for TSNE or UMAP. [default= 20]
    --marker_method=FIND_MAKER_METHOD
        Method for finding siginificant markers. [default= wilcox]
    --logFC_thre=THRESHOLD_LOGFC
        logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. This parameter was same as the 'logfc.threshold' of FindAllMarkers in Seurat R package. [default= 0.25]
    --min_pct=MIN_PCT
        Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations will be retain. This parameter was same as the 'min.pct' of FindAllMarkers in Seurat R package. [default= 0.1]
    --shown_markers=SHOWN_MAEKERS
        The number of markers for each cell type used for heatmap visualization. [default= 3]
    -o OUT_DIR, --out_dir=OUT_DIR
        Output folder. [default= NULL]
    -h, --help
        Show this help message and exit
```

#### Details for Specific Parameters
**`--normalize`:**
***log:*** Log transform. Feature counts output for each cell is divided by the total counts for that cell and multiplied by 1e4. This is then natural-log transformed.
***CLR:*** Centered log ratio. A commonly used Compositional Data Analysis (CoDA) transformation method.
***RC:*** Relative counts. Feature counts output for each cell are is divided by the total counts for that cell and multiplied by 1e4 (for TPM/CPM/FPKM/RPKM this value is 1e6).
***scran:*** The normalization strategy for scRNA-seq is implemented based on the deconvolutional size factor using the scran R package. Detials see [scran](https://github.com/MarioniLab/scran)
***none***: Do not implement normalization

**`--imputation`:**
***scImpute***: Imputing missing value of data matrix following filtering and normalization steps and this function is performed using scImpute R package
***none***: Do not implement imputation.

**`--data_type`:**
***count***: Discrete Data.
***TPM***: Continuous data.

**`--pathway_database`:**
when "--species" is "homo", "--pathway_database" can be select as follow:
***kegg***: An encyclopaedia for genes reaction and regulation. [KEGG](https://www.genome.jp/kegg/). 
***reactome***: A curated database for biomolecular pathways. [Reactome](https://reactome.org/). 
***biocarta***: A pathway database for gene regulation. [BioCarta](https://www.liebertpub.com/doi/pdf/10.1089/152791601750294344). 
***smpdb***: A small molecules pathway database. [SMPDB](https://smpdb.ca/). 
***humancyc***: A curated pathway database of human metabonomics. [HumanCyc](https://humancyc.org/). 
***panther***: A curated pathway database for protein annotation through evolutionary relationship. [PANTHER](http://www.pantherdb.org/). 
***pharmgkb***: A curated pathway database for pharmacogenomics. [pharmGKB](https://www.pharmgkb.org/). 
***acsn2***: A web-based resource depicting signalling and regulatory molecular processes in cancer cell and tumor microenvironment. [ACSN v2.0](https://acsn.curie.fr/ACSN2/ACSN2.html). 
***rb***: A curated map of molecular interactions about retinoblastoma protein (RB/RB1). [RB-Pathways](http://bioinfo-out.curie.fr/projects/rbpathway/). 
***h.all***: Hallmark gene sets. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c2.cgp***: Chemical and genetic perturbations. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c2.cp***: 
***c4.cgn***: Cancer gene neighborhoods. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c4.cm***: Cancer modules. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c5.bp***: GO biological process. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c5.mf***: GO cellular component. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c5.cc***: GO molecular function. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c6.all***: Oncogenic signatures. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 
***c7.all***: Immunologic signatures. [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp). 

when "--species" is mus, "--pathway_database" can be select as follow:
***kegg***: An encyclopaedia for genes reaction and regulation. [KEGG](https://www.genome.jp/kegg/). 
***reactome***: A curated database for biomolecular pathways. [Reactome](https://reactome.org/). 
***smpdb***: A small molecules pathway database. [SMPDB](https://smpdb.ca/). 
***c5.bp***: GO biological process. [GSKB](http://ge-lab.org/gskb/). 
***c5.mf***: GO cellular component. [GSKB](http://ge-lab.org/gskb/). 
***c5.cc***: GO molecular function. [GSKB](http://ge-lab.org/gskb/). 
***other***: Including "Location", "HPO", "STITCH", "MPO", "T3DB", "PID", "MethyCancer" and "MethCancerDB*, details see [table](http://ge-lab.org/gskb/Table%201-sources.pdf). [GSKB](http://ge-lab.org/gskb/). 
