# Open4Gene [![DOI](https://zenodo.org/badge/748825666.svg)](https://zenodo.org/doi/10.5281/zenodo.12768471)
Hurdle Model-based Method for Peak-to-Gene Linkage Analysis.
<img src="https://hbliulab.org/image/Resources/Open4Gene.Detail.png" alt="isolated" width="600"/>

## Main features of Open4Gene
- Accounting for excess zeros in single-nucleus RNA data based on a two-component mixture Hurdle model
- Modeling linkages between peak open chromatin (ATAC) and gene expression (RNA) using regression model with covariates
- Flexible specification of analysis using cells from a given cell type of interest, each cell type or all cells

## How to build & install in R (>= 4.1.0)
```r
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("GenomicRanges", quietly = TRUE)) BiocManager::install("GenomicRanges")
devtools::install_github("hbliu/Open4Gene")
library(Open4Gene)
```

## Run Open4Gene with example data
1. Load example data (./Open4Gene/inst/extdata/Open4Gene.Test.Data)
```r
load(system.file("extdata", "Open4Gene.Test.Data", package = "Open4Gene"))

ls()
# "ATAC.Counts", "RNA.Counts", "Gene.Annotation", "Meta.Data", "Peak.Gene"

head(Peak.Gene)
#Peak                   Gene
#chr5-39400433-39402082 DAB2
#chr5-39400433-39402082 DAB2
```


2. Preparing the object for Open4Gene analysis
```r
Open4Gene.obj <- CreateOpen4GeneObj(RNA = RNA.Counts,
                                    ATAC = ATAC.Counts,
                                    Meta.data = Meta.Data,
                                    Peak2Gene.Pairs = Peak.Gene,
                                    Covariates = c("lognCount_RNA","percent.mt"),
                                    Celltypes = "Cell_Type")
```

3. Run Open4Gene analysis
```r
Open4Gene.obj <- Open4Gene(object = Open4Gene.obj,
                          Celltype = "All", # Other options: Cell type name, e.g., "PT"; or "Each" to analyze each cell type
                          Binary = FALSE,   # Binarize ATAC data if binary = TRUE
                          MinNum.Cells = 5)   # Minimal number of cells with both RNA > 0 and ATAC > 0 for association test
```

4. Output Open4Gene result
```r
write.table(Open4Gene.obj@Res, file = "Open4Gene.obj.All.res.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
```

## Output
On output, Open4Gene.obj@Res provides the following values for each gene~peak pair:
1. Peak
2. Gene
3. Celltype (Cell type name, "All" means all cell in the input data)
4. TotalCellNum (Total number of cells used for the analysis)
5. ExpressCellNum (Number of cells expressing given gene, RNA read count > 1)
6. OpenCellNum (Number of cells with open chromatin in given peak, ATAC read count > 1)
7. hurdle.res.zero.beta (beta value of zero hurdle model, binomial with logit link) *
8. hurdle.res.zero.se (standard error value of zero hurdle model)
9. hurdle.res.zero.z (z value of zero hurdle model)
10. hurdle.res.zero.p (p value of zero hurdle model) *
11. hurdle.res.count.beta (beta value of count model, truncated negbin with log link) *
12. hurdle.res.count.se (standard error value of count model)
13. hurdle.res.count.z (z value of count model)
14. hurdle.res.count.p (p value of count model) *
15. hurdle.AIC (Akaike information criterion of hurdle model)
16. hurdle.BIC (Bayesian information criterion of hurdle model)
17. spearman.rho (Spearman's rank correlation coefficient between RNA and ATAC)
18. spearman.p (Spearman's p value between RNA and ATAC)

You may need columns with (*) for the downstream analysis.


## Input
This section describes how to prepare the input files for Open4Gene.
Open4Gene needs following input data and parameters:
- **RNA** [dgCMatrix] Sparse matrix of scRNAseq read count, gene (like DAB2) in row and cell in column
- **ATAC** [dgCMatrix] Sparse matrix of scATACseq read count, peak (like chr5-39400433-39402082) in row and cell in column
- **Meta.data** [data.frame] Meta data table with covariates; with cell ID in the rownames
- **Covariates** [character] Assign covariates that are needed for the analysis
- **Celltypes** [character] Assign celltype column from Meta.data
- **Peak2Gene.Pairs** [data.frame] Dataframe that contains Peak~Gene pairs to analyze using Open4Gene
- **Gene.Annotation** [GRanges] Gene annotation, e.g. EnsDb.Hsapiens.v75
- **Peak2Gene.Dis** [integer] Distance (peak to gene body), default is 100000 (bp)


## Data preparation

**1. RNA and ATAC read count matrix**

Code for extracting count matrix and cell information from Seurat object with both RNA and ATAC assays:
```r
RNA.Counts <- Seurat.object@assays$RNA@counts
ATAC.Counts <- Seurat.object@assays$ATAC@counts
Meta.Data <- Seurat.object@meta.data
```
Note here that, the cell IDs from different matrix should match with each other, e.g. columns of RNA matrix, columns of ATAC matrix, and rows of Meta.Data.


**2. Meta.Data**

This is a dataframe that contains cell information from the single cell RNA and ATAC, as following. 
|                           | orig.ident | lognCount_RNA | percent.mt | Cell_Type  |
| ------------------------- | ---------- | ------------- | ---------- | ---------- |
| HK2888_GTTTAACCAGCTCAAC-1 | HK2888     | 6.61          | 5.70       | PT         |
| HK2888_GGCTAGTGTCATGCCC-1 | HK2888     | 7.66          | 2.36       | Injured_PT |
| HK2888_GGCTAGTGTAAGCTCA-1 | HK2888     | 7.49          | 0.04       | LOH        |

Note here that, the cell IDs should be provided in the rownames. Make and check it using
```r
Meta.Data <- Seurat.object@meta.data
head(rownames(Meta.Data))
```

**3. Covariates**

The covariates are used for the regression analysis.
Features of cells in metadata can be used as covariates, e.g. lognCount_RNA, percent.mt.


**4. Peak2Gene.Pairs**

This is a dataframe that contains Peak~Gene pairs for Open4Gene, Peak (first column) and Gene (second column), as following.
| Peak | Gene                   |
| ---------------------- | ---- |
| chr5-39400433-39402082 | DAB2 |
| chr5-39369336-39370159 | DAB2 |


**5. Preparing the object for Open4Gene analysis using Gene.Annotation**

Open4Gene can pick up the Peak2Gene.Pairs based on peak and gene distance (Peak2Gene.Dis) based on the input data.
Here, the Gene.Annotation is a gene annotation in GRanges object, e.g. EnsDb.Hsapiens.v75.
Code for preparing an object for Open4Gene analysis using Gene.Annotation of EnsDb.Hsapiens.v75.

Note: bedtools should be installed and in your PATH for this function.

```r
if (!require("bedtoolsr", quietly = TRUE)) devtools::install_github("PhanstielLab/bedtoolsr")
if (!require("EnsDb.Hsapiens.v75", quietly = TRUE)) BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
Gene.Annotation <- genes(EnsDb.Hsapiens.v75)
Open4Gene.obj <- CreateOpen4GeneObj(RNA = RNA.Counts, ATAC = ATAC.Counts, Meta.data = Meta.Data,
                            Gene.Annotation = Gene.Annotation,
                            Peak2Gene.Dis = 100000,
                            Covariates = c("lognCount_RNA","percent.mt"), 
                            Celltypes = "Cell_Type")
```

Then, Open4Gene.obj will be input to run the Open4Gene analysis for all potential Peak2Gene.Pairs based on peak and gene distance.
```r
Open4Gene.obj <- Open4Gene(object = Open4Gene.obj,
                          Celltype = "All",
                          Binary = FALSE,
                          MinNum.Cells = 5)
```

The Peak2Gene.Pairs can be extracted before running Open4Gene as follows, and used to control the number of pairs analyzed in each chunk after splitting:
```r
Open4Gene.obj <- Extract.Peak2Gene.Pairs(Open4Gene.obj)
Peak2Gene.Pairs <- Open4Gene.obj@Peak2Gene.Pairs
colnames(Peak2Gene.Pairs) <- c("Peak","Gene")
head(Peak2Gene.Pairs)
#Peak                   Gene
#chr5-39185890-39186570 C9
#chr5-39187099-39187586 C9
```

## Run Open4Gene

**1. Run Open4Gene for a given cell type, e.g. PT**

```r
Open4Gene.obj <- Open4Gene(object = Open4Gene.obj,
                          Celltype = "PT",  
                          Binary = FALSE,
                          MinNum.Cells = 5)
```

**2. Run Open4Gene for each cell type**

```r
Open4Gene.obj <- Open4Gene(object = Open4Gene.obj,
                          Celltype = "Each",  
                          Binary = FALSE,
                          MinNum.Cells = 5)
```

**3. Run Open4Gene using all cells**

```r
Open4Gene.obj <- Open4Gene(object = Open4Gene.obj,
                          Celltype = "All",  
                          Binary = FALSE,
                          MinNum.Cells = 5)
```


## Warning
It is time-consuming to run genome-wide peak-to-gene linkage analysis using Open4Gene.
Open4Gene analysis on 3000 pairs takes about 5 hours.
To perform genome-wide analysis, we recommend using Peak2Gene.Pairs to control the number of pairs analyzed in each chunk.

## Other useful links
:Kidney Disease Genetic Scorecard: https://susztaklab.com/GWAS2M/index.php

## Citation
Hongbo Liu et al., Kidney multiome-based genetic scorecard reveals convergent coding and regulatory variants. Science (2025). PMID: 39913582 (https://www.science.org/doi/full/10.1126/science.adp4753)

## Contact
For any question, you are welcome to report your issue in Github or contact us hongbo919@gmail.com and ksusztak@pennmedicine.upenn.edu




