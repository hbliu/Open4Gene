#library(pscl)  ## hurdle
#library(progress)
#library(GenomicRanges)
#library(Matrix)
#library(dplyr)
#library(bedtoolsr)

#' Data Checking for CreateOpen4GeneObject
DataCheck_function <- function(object){
	Error <- character()

	#Check number of cells
	RNA_Cell_Num <- ncol(object@RNA)
	ATAC_Cell_Num <- ncol(object@ATAC)

	Gene_Num <- lengths(object@RNA@Dimnames)[1]
	num_Peaks <- lengths(object@ATAC@Dimnames)[1]

	#Check the number of cells in RNA and ATAC matrix.
	if(ncol(object@RNA) != ncol(object@ATAC)){
		Error.term <- paste("Error: The numbers of cells in ATAC matrix and RNA matrix are different.", "Please check and input the Paired data")
		Error <- c(Error, Error.term)
	}

	#Check Peak2Gene.Pairs
	if(!(length(object@Peak2Gene.Pairs) == 0)){
		if(!all(object@Peak2Gene.Pairs[[1]] %in% rownames(object@ATAC))){
			Error.term <- paste("Warning: The Peak2Gene.Pairs include Peaks which are NOT available in ATAC matrix.", "Please take care of these Peaks.")
			Error <- c(Error, Error.term)
		}
		if(!all(object@Peak2Gene.Pairs[[2]] %in% rownames(object@RNA))){
			Error.term <- paste("Warning: The Peak2Gene.Pairs include Genes which are NOT available in RNA matrix.", "Please take care of these Genes.")
			Error <- c(Error, Error.term)
		}

	}

	#Use the Gene.Annotation and Peak position to extract the Peak2Gene.Pairs if Peak2Gene.Pairs is not given:
	if(length(object@Peak2Gene.Pairs) == 0){
		if(length(object@Gene.Annotation) == 0){
			Error.term <- paste("Please provide matched Gene annotation (E.x. EnsDb.Hsapiens.v75) which is needed to extract the Peak~Gene Pairs.")
			Error <- c(Error, Error.term)
		}
	}
	
	#### Data Checking Report
	if (length(Error) == 0) TRUE else Error
}


#' Open4Gene Class
#'
#' @slot RNA dgCMatrix. A sparse matrix for RNA read count
#' @slot ATAC dgCMatrix. A sparse matrix for ATAC read count
#' @slot Meta.data data.frame. A meta data table with Covariates; and Cell IDs is in the rownames
#' @slot Meta.data data.frame. A meta data table with Covariates; and Cell IDs is in the rownames
#' @slot Covariates character. Assign Covariates that are needed for the analysis. Must be names that are in the columns of Meta.data
#' @slot Celltypes character. Assign Celltype column from Meta.data. Must be a name that is in the columns of Meta.data
#' @slot Peak2Gene.Pairs data.frame. A table including Peak~Gene Pairs for analysis, with Peak (1st column) and Gene (2nd column)
#' @slot Peak2Gene.Dis integer. Maximal distance (Peak to Gene body) used to extract Peak2Gene.Pairs, only used if no Peak2Gene.Pairs is given
#' @slot Gene.Annotation GRanges. Gene annotation (E.x. EnsDb.Hsapiens.v75) used to extract Peak2Gene.Pairs, only used if no Peak2Gene.Pairs is given
#' @slot Res data.frame. Table for result of association test, which is initialized as empty.
#'
#' @return Open4Gene object to use for further analysis
#' @export
CreateOpen4GeneObj <- setClass(
	Class = "Open4Gene",
	slots = c(
		RNA = 'dgCMatrix',
		ATAC = 'dgCMatrix',
		Meta.data = 'data.frame',
		Covariates = 'character',
		Celltypes = 'character',
		Peak2Gene.Pairs = 'data.frame',
		Peak2Gene.Dis = 'numeric',
		Gene.Annotation = 'GRanges',
		Cell.IDs = 'character',
		Res = 'data.frame'
	),
	prototype = list(
		Peak2Gene.Pairs = data.frame(),
		Peak2Gene.Dis = 100000,
		Celltypes = "Cell_Type"
	),
	validity = DataCheck_function
)


#' Open4Gene: Main function
#'
#' @param object Open4Gene object
#' @param Celltype User specified cell type defined in Celltypes column of Meta.data. "All" to test using all cells. "Each" to run test for each cell type. Others, the test will be perform for a given cell type.
#' @param Binary If TRUE, the ATAC > 1 will be converted as 1.
#' @param Method "hurdle" (default) or "fasthurdle" for Zero-inflated Negative Binomial Regression based on hurdle model. Statistical method used to calculate the correlation between ATAC and RNA. 
#' @param MinNum.Cells Minimal number of cells with expression (RNA > 0) and open chromatin (ATAC > 0) for association test.
#' @return Open4Gene object with Results from hurdle model
#' @export
Open4Gene <- function(object, Celltype = "All", Binary = FALSE, Method = "hurdle", MinNum.Cells = 5){
	#### Extract Peak~Gene Pairs used for regression
	print('Prepare Peak~Gene Pairs for regression analysis...', quote = FALSE)
	object <- Extract.Peak2Gene.Pairs(object)
	print(c(paste("Start regression analysis for", nrow(object@Peak2Gene.Pairs), "Peak~Gene Pairs using", Method)), quote = FALSE)
	Res <- data.frame()
	
	#### Extract RNA and ATAC used for regression
	object@RNA <- object@RNA[unique(object@Peak2Gene.Pairs[,2]), , drop = FALSE]
	object@ATAC <- object@ATAC[unique(object@Peak2Gene.Pairs[,1]), , drop = FALSE]

	#### Define progress bar
	Probar <- progress_bar$new(total=nrow(object@Peak2Gene.Pairs), format = 'Processing [:bar] :current/:total (:percent) eta: :eta', clear = FALSE, width = 80)
	
	#### Define Formula for regression
	Var.Res <- "RNA"
	Var.Factor <- c("ATAC", object@Covariates) 
	Formula <- as.formula(paste(Var.Res, paste(Var.Factor, collapse = "+"), sep = "~"))
	for (i in 1:length(object@Covariates)){
		object@Meta.data[,object@Covariates[i]] <- as.integer(object@Meta.data[,object@Covariates[i]])
	}
	
	for (n in 1:nrow(object@Peak2Gene.Pairs)){
		Probar$tick()
		Peak <- object@Peak2Gene.Pairs[n,1] ### Peak is in the first column in Peak2Gene.Pairs
		Gene <- object@Peak2Gene.Pairs[n,2] ### Gene is in the second column in Peak2Gene.Pairs
		ATAC_data <- data.frame(Cell.ID = colnames(object@ATAC), ATAC = object@ATAC[Peak,])
		Meta.data <- object@Meta.data
		Meta.data$Cell.ID <- rownames(Meta.data)
		
		#ATAC data will be binarized if Binary = TRUE
		if(Binary){
			ATAC_data[ATAC_data$ATAC > 0,]$ATAC <- 1
		}

		Gene_exp <- object@RNA[Gene,]
		DatMat <- data.frame(Cell.ID=names(Gene_exp),RNA=as.integer(Gene_exp))
		DatMat <- merge(DatMat, ATAC_data, by="Cell.ID")
		DatMat <- merge(DatMat, Meta.data, by="Cell.ID")
		Celltypes <- as.character(unique(DatMat[[object@Celltypes]]))
	
		############# Run test according to the selection of cell types
		if(Celltype == "All"){	# Run test using all cells
			if(length(DatMat$RNA[DatMat$RNA == 0]) > 0 & length(DatMat$RNA[DatMat$RNA > 0]) >= MinNum.Cells & length(DatMat$ATAC[ DatMat$ATAC > 0]) >= MinNum.Cells){
				Res <- AssociationTest(DatMat, Gene, Peak, Method, Formula, Celltype, Res)
			}
		}else if(Celltype == "Each"){# Run test for each cell type
			for (c in 1:length(Celltypes)){
				SubDatMat <- DatMat[DatMat[[object@Celltypes]] == Celltypes[c],]
				if(length(SubDatMat$RNA[SubDatMat$RNA == 0] ) > 0 & length(SubDatMat$RNA[ SubDatMat$RNA > 0]) >= MinNum.Cells & length(SubDatMat$ATAC[ SubDatMat$ATAC > 0]) >= MinNum.Cells){
					Res <- AssociationTest(SubDatMat, Gene, Peak, Method, Formula, Celltypes[c], Res)
				}
			}
		}else{ # Run test only for a given cell type
			SubDatMat <- DatMat[DatMat[[object@Celltypes]] == Celltype,]
			if(length(SubDatMat$RNA[ SubDatMat$RNA == 0]) > 0 & length(SubDatMat$RNA[ SubDatMat$RNA > 0]) >= MinNum.Cells & length(SubDatMat$ATAC[ SubDatMat$ATAC > 0]) >= MinNum.Cells){
				Res <- AssociationTest(SubDatMat, Gene, Peak, Method, Formula, Celltype, Res)
			}
		}
	}## for

	#Add new Result into object
	object@Res <- Res
	return(object)

}



#' AssociationTest:
#' @param DatMat Data Matrix including RNA count, ATAC count and meta data
#' @param Method Statistic Method used for the association Test
#' @return Res Statistic Result
AssociationTest <- function(DatMat, Gene, Peak, Method, Formula, Celltype, Res){

	### Association analysis using Spearman correlation
	Spearman <- cor.test(DatMat$ATAC, DatMat$RNA, method = "spearman",exact=FALSE)

	### Association analysis using Zero-inflated Negative Binomial Regression based on hurdle model
	if(Method == "hurdle"){
		# print(Method)
		### hurdle is faster than zeroinfl, with similar performance
		hurdle.test <- hurdle(Formula, data = DatMat, link = "logit", dist = "negbin")
		hurdle.Res <- summary(hurdle.test)
		hurdle.Res.count <- hurdle.Res$coefficients$count["ATAC",]
		hurdle.Res.zero <- hurdle.Res$coefficients$zero["ATAC",]
		out <- data.frame(Peak=Peak, Gene=Gene, Celltype = Celltype, TotalCellNum=nrow(DatMat), 
							ExpRessCellNum=length(DatMat$RNA[DatMat$RNA > 0]),OpenCellNum = length(DatMat$ATAC[ DatMat$ATAC > 0]),
							hurdle.Res.zero.beta=hurdle.Res.zero[1],hurdle.Res.zero.se=hurdle.Res.zero[2],hurdle.Res.zero.z=hurdle.Res.zero[3],hurdle.Res.zero.p=hurdle.Res.zero[4],
							hurdle.Res.count.beta=hurdle.Res.count[1],hurdle.Res.count.se=hurdle.Res.count[2],hurdle.Res.count.z=hurdle.Res.count[3],hurdle.Res.count.p=hurdle.Res.count[4],
							hurdle.AIC=AIC(hurdle.test),hurdle.BIC=BIC(hurdle.test),
							spearman.rho = Spearman$estimate, spearman.p = Spearman$p.value
							)
		out[,c(7:9,11:13,17)] <- round(out[,c(7:9,11:13,17)],digits = 6)
		out[,c(10,14,18)] <- signif(out[,c(10,14,18)],digits = 6)
	}

	### Association analysis using Zero-inflated Negative Binomial Regression based on hurdle model, implemented by fasthurdle (https://github.com/mkanai/fasthurdle)
	if(Method == "fasthurdle"){
		#print(Method)
		if (!require("devtools", quietly = TRUE)) install.packages("devtools")
		if (!require("fastglm", quietly = TRUE)) install.packages("fastglm")
		if (!require("fasthurdle", quietly = TRUE)) devtools::install_github("mkanai/fasthurdle")
		library(fasthurdle)
		### fasthurdle is faster than hurdle, with similar performance
		hurdle.test <- fasthurdle(Formula, data = DatMat, link = "logit", dist = "negbin")
		hurdle.Res <- summary(hurdle.test)
		hurdle.Res.count <- hurdle.Res$coefficients$count["ATAC",]
		hurdle.Res.zero <- hurdle.Res$coefficients$zero["ATAC",]
		out <- data.frame(Peak=Peak, Gene=Gene, Celltype = Celltype, TotalCellNum=nrow(DatMat), 
							ExpRessCellNum=length(DatMat$RNA[DatMat$RNA > 0]),OpenCellNum = length(DatMat$ATAC[ DatMat$ATAC > 0]),
							hurdle.Res.zero.beta=hurdle.Res.zero[1],hurdle.Res.zero.se=hurdle.Res.zero[2],hurdle.Res.zero.z=hurdle.Res.zero[3],hurdle.Res.zero.p=hurdle.Res.zero[4],
							hurdle.Res.count.beta=hurdle.Res.count[1],hurdle.Res.count.se=hurdle.Res.count[2],hurdle.Res.count.z=hurdle.Res.count[3],hurdle.Res.count.p=hurdle.Res.count[4],
							hurdle.AIC=AIC(hurdle.test),hurdle.BIC=BIC(hurdle.test),
							spearman.rho = Spearman$estimate, spearman.p = Spearman$p.value
							)
		out[,c(7:9,11:13,17)] <- round(out[,c(7:9,11:13,17)],digits = 6)
		out[,c(10,14,18)] <- signif(out[,c(10,14,18)],digits = 6)
	}
	
	### Combine Result
	Res<-rbind(Res,out)
	return(Res)
}



#' Extract.Peak2Gene.Pairs: Extracting Peak~Gene Pairs using bedtoolsr
#'
#' @param object Open4Gene object
#' @param object Gene.Annotation GRanges
#' @return Open4Gene object with Peak2Gene.Pairs extracted
#' @export
Extract.Peak2Gene.Pairs <- function(object){
	if(length(object@Peak2Gene.Pairs) > 0){
			colnames(object@Peak2Gene.Pairs) <- c("Peak", "Gene")
			object@Peak2Gene.Pairs <- subset(object@Peak2Gene.Pairs, Peak %in% rownames(object@ATAC) & Gene %in% rownames(object@RNA))
			if(length(object@Peak2Gene.Pairs) == 0){
				Error.term <- paste("Please provide the Peak~Gene pairs available in RNA and ATAC data.")
					Error <- c(Error, Error.term)
			}
	}else{
		Gene.Annotation <- object@Gene.Annotation
		Peak2Gene.Dis <- object@Peak2Gene.Dis
		ATAC <- object@ATAC
		Peak2Gene.Pairs <- object@Peak2Gene.Pairs
		if(length(Gene.Annotation) == 0){
			Error.term <- paste("Please provide a valid gene annotation to extract Peak~Gene pairs.")
			Error <- c(Error, Error.term)
		}else if(length(Gene.Annotation) > 0 & Peak2Gene.Dis >= 0){
			options(dplyr.summarise.inform = FALSE)
			seqlevelsStyle(Gene.Annotation) <- "UCSC"
			Gene <- data.frame(Gene.Annotation)
			Gene <- subset(Gene, gene_name %in% rownames(object@RNA))
			Gene <- Gene %>% group_by(gene_name) %>% summarise(seqnames = unique(seqnames), start = min(start, na.rm=TRUE), end = max(end, na.rm=TRUE)) %>% as.data.frame
			if(length(Gene) == 0){
					Error.term <- paste("Please provide genes available in RNA data")
					Error <- c(Error, Error.term)
			}else{
				Gene.Pos <- Gene[,c(2:4,1)]
				Gene.Pos$start <- Gene.Pos$start - Peak2Gene.Dis
				Gene.Pos$start[Gene.Pos$start < 0] <- 0
				Gene.Pos$end <- Gene.Pos$end + Peak2Gene.Dis
			
				Peak <- rownames(ATAC)
				Peak <- data.frame(do.call(rbind, strsplit(Peak, "-", fixed=TRUE)))
				Peak$Peak <- paste(Peak$X1, Peak$X2, Peak$X3, sep = '-')
						
				Gene.Peak <- bedtoolsr::bt.intersect(Peak, Gene.Pos, wa = TRUE, wb = TRUE)
				object@Peak2Gene.Pairs <- Gene.Peak[,c(4,8)]
				if(length(Gene.Peak) == 0){
						Error.term <- paste("There is no Peak~Gene pairs within given distance. Please check data and enlarger Peak2Gene.Dis.")
						Error <- c(Error, Error.term)
				}
			}
		}
	}
	return(object)
}


