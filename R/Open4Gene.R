#library(pscl)  ## hurdle
#library(progress)
#library(GenomicRanges)
#library(Matrix)
#library(dplyr)
#library(bedtoolsr)

#' DataCheck and Type Checking for CreateOpen4GeneObject Constructor
DataCheck_function <- function(object){
	errors <- character()

	#Check number of cells
	RNA_Cell_Num <- ncol(object@RNA)
	ATAC_Cell_Num <- ncol(object@ATAC)

	Gene_Num <- lengths(object@RNA@Dimnames)[1]
	num_peaks <- lengths(object@ATAC@Dimnames)[1]

	#Check if the number of cells match between RNA and ATAC matrix.
	if(ncol(object@RNA) != ncol(object@ATAC)){
		msg <- paste("Error: The num of cells in RNA matrix is different with the number of cells in ATAC matrix.", "Please check the data.")
		errors <- c(errors, msg)
	}

	#Check gene.peak.pair
	if(!(length(object@gene.peak.pair) == 0)){
		if(!all(object@gene.peak.pair[[1]] %in% rownames(object@RNA))){
			msg <- paste("Some genes in the gene.peak.pair is NOT available in RNA matrix.", "Please check the data.")
			errors <- c(errors, msg)
		}

		if(!all(object@gene.peak.pair[[2]] %in% rownames(object@ATAC))){
			msg <- paste("Some peaks in the gene.peak.pair is NOT available in ATAC matrix.", "Please check the data.")
			errors <- c(errors, msg)
		}
	}

	#If gene.peak.pair is not present, use the gene.annotation and peak position to extract the gene.peak.pair:
	if(length(object@gene.peak.pair) == 0){
		if(length(object@gene.annotation) == 0){
			msg <- paste("The gene annotation is needed to extract the gene~peak pairs.")
			errors <- c(errors, msg)
		}
	}
	
	#### Summary of checking
	if (length(errors) == 0) TRUE else errors
}


#' Open4Gene Class Constructor
#'
#' @slot RNA dgCMatrix. scRNAseq matrix read as a sparse matrix
#' @slot ATAC dgCMatrix. scATACseq matrix read as a sparse matrix
#' @slot meta.data data.frame. Metadata table with covariates and a cell ID column ("cell")
#' @slot gene.peak.pair data.frame. Dataframe that contains gene-peak pairs for Open4Gene to search through
#' @slot gene.peak.dis integer. Distance (peak to gene body) used to extract gene.peak.pair
#' @slot gene.annotation GRanges. Gene annotation used to extract gene.peak.pair
#' @slot covariates character. Assign covariates that are needed for the analysis. Must be names that are in the columns of meta.data
#' @slot celltypes character. Assign celltype column from meta.data
#' @slot res data.frame. Table for result of association test, which is initialized as empty.
#'
#' @return Open4Gene object to use for further analysis
#' @export
CreateOpen4GeneObj <- setClass(
	Class = "Open4Gene",
	slots = c(
		RNA = 'dgCMatrix',
		ATAC = 'dgCMatrix',
		meta.data = 'data.frame',
		gene.peak.pair = 'data.frame',	### A table including gene~peak pairs for analysis, with gene (1st column) and peak (2nd column).
		gene.peak.dis = 'numeric',		### Distance (peak to gene body) used to extract gene.peak.pair. This is only used if no gene.peak.pair is detected. Please provide this and gene.annotation if you don't provide gene.peak.pair for test.
		gene.annotation = 'GRanges',	### Gene annotation used to extract gene.peak.pair
		covariates = 'character',
		celltypes = 'character',
		res = 'data.frame'
	),
	prototype = list(
		gene.peak.pair = data.frame(),
		gene.peak.dis = 100000,
		celltypes = "All"
	),
	validity = DataCheck_function
)


#' Open4Gene: main function
#'
#' @param object Open4Gene object
#' @param celltype User specified cell type defined in celltypes column of meta.data. "All" to test using all cells. "Each" to run test for each cell type. Others, the test will be perform for a given cell type.
#' @param binary If TRUE, the ATAC > 1 will be converted as 1.
#' @param method Statistic method used to calculate the correlation between ATAC and RNA. "hurdle" for Zero-inflated Negative Binomial Regression based on hurdle model.
#' @param MinCellNum Minimal number of cells with expression (RNA > 0) and cells with open chromatin (ATAC > 0) for association test.
#'
#' @return Open4Gene object with updated field ress
#' @export
Open4Gene <- function(object, celltype = "All", binary = FALSE, method = "hurdle", MinCellNum = 5){
	# Extract gene~peak pairs used for regression
	print('Prepare gene~peak pairs for regression analysis...', quote = FALSE)
	object <- Extract.gene.peak.pair(object)
	print(c(paste("Start regression analysis for",nrow(object@gene.peak.pair), "gene~peak pairs...")), quote = FALSE)
	res <- data.frame()
	
	# Extract RNA and ATAC used for regression
	object@RNA <- object@RNA[unique(object@gene.peak.pair[,1]),]
	object@ATAC <- object@ATAC[unique(object@gene.peak.pair[,2]),]

	#### Define progress bar
	pb <- progress_bar$new(total=nrow(object@gene.peak.pair), format = 'Processing [:bar] :current/:total (:percent) eta: :eta', clear = FALSE, width = 80)
	
	#### Define formula for regression
	res_var <- "RNA"
	pred_var <- c("ATAC", object@covariates) 
	formula <- as.formula(paste(res_var, paste(pred_var, collapse = "+"), sep = "~"))
	for (i in 1:length(object@covariates)){
		object@meta.data[,object@covariates[i]] <- as.integer(object@meta.data[,object@covariates[i]])
	}
	
	for (n in 1:nrow(object@gene.peak.pair)){
		pb$tick()
		gene <- object@gene.peak.pair[n,1] #Gene is in the first column in gene.peak.pair
		peak <- object@gene.peak.pair[n,2] #Peak is in the second column in gene.peak.pair
		ATAC_data <- data.frame(cell = colnames(object@ATAC), ATAC = object@ATAC[peak,])	### This step takes time if extracting peak from full dataset
		meta.data <- object@meta.data
		if(!("cell" %in% colnames(meta.data))){       ## use rownames as cell id if "cell" column is not available.
			meta.data$cell <- rownames(meta.data)
		}
		
		#Binarize ATAC data if binary = TRUE
		if(binary){
			ATAC_data[ATAC_data$ATAC > 0,]$ATAC <- 1
		}

		Gene_exp <- object@RNA[gene,]
		DataM <- data.frame(cell=names(Gene_exp),RNA=as.integer(Gene_exp))
		DataM <- merge(DataM, ATAC_data, by="cell")
		DataM <- merge(DataM, meta.data, by="cell")
		Celltypes <- as.character(unique(DataM[[object@celltypes]]))
	
		############# Run test according to the selection of cell types
		if(celltype == "All"){	# Run test using all cells
			if(length( DataM$RNA[ DataM$RNA == 0] ) > 0 & length( DataM$RNA[ DataM$RNA > 0] ) >= MinCellNum & length( DataM$ATAC[ DataM$ATAC > 0] ) >= MinCellNum){
				res <- AssociationTest(DataM, gene, peak, method, formula, celltype, res)
			}
		}else if(celltype == "Each"){# Run test for each cell type
			for (c in 1:length(Celltypes)){
				SubDataM <- DataM[DataM[[object@celltypes]] == Celltypes[c],]
				if(length( SubDataM$RNA[ SubDataM$RNA == 0] ) > 0 & length( SubDataM$RNA[ SubDataM$RNA > 0] ) >= MinCellNum & length( SubDataM$ATAC[ SubDataM$ATAC > 0] ) >= MinCellNum){
					res <- AssociationTest(SubDataM, gene, peak, method, formula, Celltypes[c], res)
				}
			}
		}else{ # Run test only for a given cell type
			SubDataM <- DataM[DataM[[object@celltypes]] == celltype,]
			if(length( SubDataM$RNA[ SubDataM$RNA == 0] ) > 0 & length( SubDataM$RNA[ SubDataM$RNA > 0] ) >= MinCellNum & length( SubDataM$ATAC[ SubDataM$ATAC > 0] ) >= MinCellNum){
				res <- AssociationTest(SubDataM, gene, peak, method, formula, celltype, res)
			}
		}
	}## for

	#Add new result into object
	object@res <- res
	return(object)
}



#' AssociationTest:
#' @param DataM Data Matrix including RNA count, ATAC count and meta data
#' @param method Statistic method used for association Test
#' @return res Statistic result
AssociationTest <- function(DataM, gene, peak, method, formula, celltype, res){
	Spearman <- cor.test(DataM$ATAC, DataM$RNA, method = "spearman",exact=FALSE)

	############# Zero-inflated Negative Binomial Regression based on hurdle model ##############
	if(method == "hurdle"){
		### hurdle is faster than zeroinfl, with similar performance
		hurdle.test <- hurdle(formula, data = DataM, link = "logit", dist = "negbin")
		hurdle.res <- summary(hurdle.test)
		hurdle.res.count <- hurdle.res$coefficients$count["ATAC",]
		hurdle.res.zero <- hurdle.res$coefficients$zero["ATAC",]
		# Count model coefficients (truncated negbin with log link):
		# Zero hurdle model coefficients (binomial with logit link):

		out <- data.frame(gene=gene,peak=peak, Celltype = celltype, TotalCellNum=nrow(DataM), ExpressCellNum=length( DataM$RNA[ DataM$RNA > 0] ),OpenCellNum = length( DataM$ATAC[ DataM$ATAC > 0] ),
							hurdle.res.zero.beta=hurdle.res.zero[1],hurdle.res.zero.se=hurdle.res.zero[2],hurdle.res.zero.z=hurdle.res.zero[3],hurdle.res.zero.p=hurdle.res.zero[4],
							hurdle.res.count.beta=hurdle.res.count[1],hurdle.res.count.se=hurdle.res.count[2],hurdle.res.count.z=hurdle.res.count[3],hurdle.res.count.p=hurdle.res.count[4],
							hurdle.AIC=AIC(hurdle.test),hurdle.BIC=BIC(hurdle.test),
							spearman.rho = Spearman$estimate, spearman.p = Spearman$p.value
							)
		out[,c(7:9,11:13,17)] <- round(out[,c(7:9,11:13,17)],digits = 6)
		out[,c(10,14,18)] <- signif(out[,c(10,14,18)],digits = 6)
	}
	
	#### Combine result
	res<-rbind(res,out)
	return(res)
}




#' Extract.gene.peak.pair: code for extracting gene~peak pairs based on gene annotation and peak locations
#'
#' @param object Open4Gene object
#' @param object gene.annotation GRanges
#' @return Open4Gene object with updated gene.peak.pair
#' @export
Extract.gene.peak.pair <- function(object){
	if(length(object@gene.peak.pair) > 0){
			colnames(object@gene.peak.pair) <- c("gene","peak")
			object@gene.peak.pair <- subset(object@gene.peak.pair, gene %in% rownames(object@RNA) & peak %in% rownames(object@ATAC))
			if(length(object@gene.peak.pair) == 0){
				msg <- paste("The gene~peak pairs provided are not available in RNA or ATAC data")
					errors <- c(errors, msg)
			}
	}else{
		gene.annotation <- object@gene.annotation
		gene.peak.dis <- object@gene.peak.dis
		ATAC <- object@ATAC
		gene.peak.pair <- object@gene.peak.pair
		if(length(gene.annotation) == 0){
					msg <- paste("The gene annotation is needed to extract the gene~peak pairs")
					errors <- c(errors, msg)
		}else if(length(gene.annotation) > 0 & gene.peak.dis >= 0){	### 
					### Gene to bed
					options(dplyr.summarise.inform = FALSE)
					seqlevelsStyle(gene.annotation) <- "UCSC"
					gene <- data.frame(gene.annotation)
					gene <- subset(gene, gene_name %in% rownames(object@RNA))
					gene <- gene %>% group_by(gene_name) %>% summarise(seqnames = unique(seqnames), start = min(start, na.rm=TRUE), end = max(end, na.rm=TRUE)) %>% as.data.frame
					if(length(gene) == 0){
							msg <- paste("No genes available in RNA dataset")
							errors <- c(errors, msg)
					}else{
						gene.ext <- gene[,c(2:4,1)]
						gene.ext$start <- gene.ext$start - gene.peak.dis
						gene.ext$start[gene.ext$start < 0] <- 0
						gene.ext$end <- gene.ext$end + gene.peak.dis
			
						### Peak to bed
						peak <- rownames(ATAC)
						peak <- data.frame(do.call(rbind, strsplit(peak, "-", fixed=TRUE)))
						peak$Peak <- paste(peak$X1, peak$X2, peak$X3, sep = '-')
			
						### bedtools
						gene.peak <- bedtoolsr::bt.intersect(gene.ext, peak, wa = TRUE, wb = TRUE)
						object@gene.peak.pair <- gene.peak[,c(4,8)]
						if(length(gene.peak) == 0){
								msg <- paste("No gene~peak pairs identified using given distance.")
								errors <- c(errors, msg)
						}
					}
			}
		}
	return(object)
}


