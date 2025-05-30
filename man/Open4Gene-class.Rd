% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Open4Gene.R
\docType{class}
\name{Open4Gene-class}
\alias{Open4Gene-class}
\alias{CreateOpen4GeneObj}
\title{Open4Gene Class}
\value{
Open4Gene object to use for further analysis
}
\description{
Open4Gene Class
}
\section{Slots}{

\describe{
\item{\code{RNA}}{dgCMatrix. A sparse matrix for RNA read count}

\item{\code{ATAC}}{dgCMatrix. A sparse matrix for ATAC read count}

\item{\code{Meta.data}}{data.frame. A meta data table with Covariates; and Cell IDs is in the rownames}

\item{\code{Meta.data}}{data.frame. A meta data table with Covariates; and Cell IDs is in the rownames}

\item{\code{Covariates}}{character. Assign Covariates that are needed for the analysis. Must be names that are in the columns of Meta.data}

\item{\code{Celltypes}}{character. Assign Celltype column from Meta.data. Must be a name that is in the columns of Meta.data}

\item{\code{Peak2Gene.Pairs}}{data.frame. A table including Peak~Gene Pairs for analysis, with Peak (1st column) and Gene (2nd column)}

\item{\code{Peak2Gene.Dis}}{integer. Maximal distance (Peak to Gene body) used to extract Peak2Gene.Pairs, only used if no Peak2Gene.Pairs is given}

\item{\code{Gene.Annotation}}{GRanges. Gene annotation (E.x. EnsDb.Hsapiens.v75) used to extract Peak2Gene.Pairs, only used if no Peak2Gene.Pairs is given}

\item{\code{Res}}{data.frame. Table for result of association test, which is initialized as empty.}
}}

