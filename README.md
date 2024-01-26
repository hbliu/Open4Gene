# Open4Gene
Hurdle Model-based Method for Peak-to-Gene Linkage Analysis.

## Main features of Open4Gene
- Accounting for excess zeros in single-nucleus RNA data based on a two-component mixture Hurdle model
- Modeling linkages between peak open chromatin (ATAC) and gene expression (RNA) using regression model with covariates
- Flexible specification of analysis using cells from a given cell type of interest, each cell type or all cells

## How to build & install
```r
devtools::install_github("hbliu/Open4Gene")
library(Open4Gene)
```

## Run Open4Gene with example data
1. Load example data (./Open4Gene/inst/extdata/Open4Gene.Test.Data)
```r
load(system.file("extdata", "Open4Gene.Test.Data", package = "Open4Gene"))
ls()
# [1] "ATAC.counts"   "RNA.counts"   "gene.annotation"   "gene_peak"   "meta"
head(gene_peak)
#Gene                   Peak
#1 DAB2 chr5-39400433-39402082
#2 DAB2 chr5-39369336-39370159
```


2. Preparing the object for Open4Gene analysis
```r
Open4Gene.obj <- CreateOpen4GeneObj(RNA = RNA.counts,
                                    ATAC = ATAC.counts,
                                    meta.data = meta,
                                    gene.peak.pair = gene_peak,
                                    covariates = c("lognCount_RNA","percent.mt"),
                                    celltypes = "Celltype")
```


3. Run Open4Gene analysis
```r
Open4Gene.obj <- Open4Gene(object = Open4Gene.obj,
                          celltype = "All", # Other options: Cell type name, e.g., "PT"; or "Each" to analyze each cell type
                          binary = FALSE,   # Binarize ATAC data if binary = TRUE
                          MinCellNum = 5)   # Minimal number of cells with both RNA > 0 and ATAC > 0 for association test
```

4. Output Open4Gene result
```r
write.table(Open4Gene.obj@res, file = "Open4Gene.obj.All.res.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
```

## Input
This section describes how to prepare the input files for Open4Gene.
Open4Gene needs following input data and parameters:
- RNA dgCMatrix. scRNAseq matrix read as a sparse matrix
- ATAC dgCMatrix. scATACseq matrix read as a sparse matrix
- meta.data data.frame. Metadata table with covariates and a cell ID column ("cell")
- gene.peak.pair data.frame. Dataframe that contains gene-peak pairs for Open4Gene to search through
- covariates character. Assign covariates that are needed for the analysis. Must be names that are in the columns of meta.data
- celltypes character. Assign celltype column from meta.data


## Output
On output, Open4Gene.obj@res provides the following values for each gene~peak pair:
1. gene
2. peak
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
