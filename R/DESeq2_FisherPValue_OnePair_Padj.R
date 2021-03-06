#' Combine Adjusted P values for one pair of input and output
#' 
#' This function is a variation of function \link[mmcc]{DESeq2_FisherPvalue_Onepair}
#' 
#' This function should be used on a raw data with the following layout:
#' 1. The first three columns must be "Feature", "Locus_Tag", "Strand"
#'    Here "Feture" can be barcode (preferred) or anything that unique. Be sure that the colnames are right!
#' 2. The next columns are for input and output successively. 
#'    Note that there is no constraint for the names of these input and
#'    output columns, but there can only be 1 kind of output. The number
#'    of columns in each output can be any number.
#' 
#' @param filepath  The path of the raw data (excluding the data file itself). Remember to double quote the filepath.
#' @param rawdatafile The name of the raw data file. Remember to double quote the rawdatafile, and add the ".xlsx" after the name of rawdatafile (It must be excel file)!
#' @param num.in The number of input columns in the raw data file.
#' @param num.out The number of output columns in the raw data file.
#' @param name.out The prefix of columns in the output. For example, "A vs B" or "B vs C", in which "A". "B" and "C" represent the names of input and output in the raw data file. The default is "Input vs Output".
#' @keywords mmcc
#' @export
#' @examples DESeq2_FisherPvalue_Onepair_Padj(filepath = "C:/Users/", rawdatafile = "rawdata.xlsx", num.in = 3, num.out = 3, name.out = "Input vs Output")

DESeq2_FisherPvalue_Onepair_Padj <- function(filepath, rawdatafile,
                                              num.in, num.out, name.out = "Input vs Output"){
  options(scipen = 999)
  setwd(filepath)
  outputname.feature <- substr(rawdatafile, 1, gregexpr("\\.xlsx", rawdatafile)[[1]][1] - 1)
  outputname.oligo <- substr(rawdatafile, 1, gregexpr("\\.xlsx", rawdatafile)[[1]][1] - 1)
  if ("DESeq2" %in% rownames(installed.packages()) == FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite("DESeq2")
	library(DESeq2)
  }else{
  	library(DESeq2)
  }
  if ("readxl" %in% rownames(installed.packages()) == FALSE){
	install.packages("readxl")
	library(readxl)
  }else{
  	library(readxl)
  }
  if ("survcomp" %in% rownames(installed.packages()) == FALSE){
	source("https://bioconductor.org/biocLite.R")
	biocLite("survcomp")
	library(survcomp)
  }else{
	library(survcomp)
  }

  data <- read_excel(rawdatafile, col_names = TRUE)
  data <- data[colSums(!is.na(data)) > 0]
  data <- as.data.frame(data)
  data$FeatureID <- paste0(as.character(data[,1]),"_", as.character(data$Locus_Tag))
  data <- data[,c(1,2,ncol(data),4:ncol(data)-1)]


  countData <- data[,5:ncol(data)]
  countData <- sapply(countData, as.numeric)
  rownames(countData) <- data[,3]
  countData <- as.data.frame(countData)
  colData <- data.frame(c(rep("A", num.in), rep("B", num.out)))
  rownames(colData) <- colnames(data)[5:length(data)]
  colnames(colData) <- c("V1")
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~V1)
  
  dds <- DESeq(dds)
  ### for output
  res <- results(dds)
  output <- cbind(data[,2], res$log2FoldChange, res$lfcSE, 
                          res$pvalue, res$padj)
  colnames(output) <- c("Locus_Tag", "log2FoldChange", "lfcSE", "pvalue", "padj")
  output <- data.frame(output)
  output[,2:5] <- sapply(sapply(output[,2:5], as.character), as.numeric)
  output$weight <- 1/(output$lfcSE^2)
  
  Oligo <- data[,1]
  output.oligo <- cbind(Oligo,output)
  
  write.csv(output.oligo, paste0(outputname.oligo, "_oligo", "_in", as.character(num.in), "_out", as.character(num.out), ".csv"))
  
  ######
  ## 1/p
  weightFC1 <- data.frame(output[,c(1,2,4,5)])
  weightFC1 <- na.omit(weightFC1)
  weightFC1 <- split(weightFC1, weightFC1[,1])
  
  weight.fc <- function (x){
    p.sum <- sum(1/x[,3])
    p.adj.sum <- sum(1/x[,4])
    x$PvalueFC <- ((1/x[,3]) / p.sum) * x[,2]
    x$PvalueAFC <- ((1/x[,4]) / p.adj.sum) * x[,2]
    return (x)
  }
  
  weightFC1 <- lapply(weightFC1, weight.fc)
  
  helper <- function(x) {
    return (x[,5:6])
  }
  
  weightFC1 <- lapply(weightFC1, helper)
  
  PFC.1 <- unlist(lapply(weightFC1, function(x) sum(x[,1])))
  APFC.1 <- unlist(lapply(weightFC1, function(x) sum(x[,2])))
  
  FC.1 <- cbind(PFC.1, APFC.1)
  rownames(FC.1) <- c(1:nrow(FC.1))
  
  ### 1 - p
  
  weightFC1_1minp <- data.frame(output[,c(1,2,4,5)])
  weightFC1_1minp <- na.omit(weightFC1_1minp)
  weightFC1_1minp <- split(weightFC1_1minp, weightFC1_1minp[,1])
  
  weight.fc_1minp <- function (x){
    p.sum <- sum(1 - x[,3])
    p.adj.sum <- sum(1 - x[,4])
    x$PvalueFC_1minp <- ((1 - x[,3]) / p.sum) * x[,2]
    x$PvalueAFC_1minp <- ((1 - x[,4]) / p.adj.sum) * x[,2]
    return (x)
  }
  
  weightFC1_1minp <- lapply(weightFC1_1minp, weight.fc_1minp)
  
  weightFC1_1minp <- lapply(weightFC1_1minp, helper)
  
  PFC.1_1minp <- unlist(lapply(weightFC1_1minp, function(x) sum(x[,1])))
  APFC.1_1minp <- unlist(lapply(weightFC1_1minp, function(x) sum(x[,2])))
  
  FC.1_1minp <- cbind(PFC.1_1minp, APFC.1_1minp)
  rownames(FC.1_1minp) <- c(1:nrow(FC.1_1minp))
  
  ######
  
  df <- output

  df$padj <- df$padj / 2
  df.trans_plus <- df
  df.trans_minus <- df
  for (j in 1:nrow(df)){
    if (!is.na(df[j,2])){
      if (df[j,2] > 0 & !is.na(df[j,5])){
        df.trans_plus[j, 5] <- 1 - df.trans_plus[j, 5]
      }else if (df[j,2] < 0 & !is.na(df[j,5])){
        df.trans_minus[j, 5] <- 1 - df.trans_minus[j, 5]
      }
    }
  }
    df2.trans_plus <- split(df.trans_plus, f = df.trans_plus[, 1])
    df2.trans_minus <- split(df.trans_minus, f = df.trans_minus[, 1])
    
    #####
    
    p.df.trans_plus <- data.frame()
    
    for (k in 1:length(df2.trans_plus)){
      set <- df2.trans_plus[[k]]
	p.df.trans_plus[k,1] <- set[1,1]
	set <- na.omit(set)
      p.df.trans_plus[k,2] <- combine.test(set[,5], method = "fisher")
      p.df.trans_plus[k,3] <- combine.test(set[,5], method = "z.transform")
      p.df.trans_plus[k,4] <- combine.test(set[,5], weight = set[,6], method = "z.transform")
    }
    
    #########

    p.df.trans_minus <- data.frame()
    
    for (k in 1:length(df2.trans_minus)){
      set <- df2.trans_minus[[k]]
      p.df.trans_minus[k,1] <- set[1,1]
	set <- na.omit(set)
      p.df.trans_minus[k,2] <- combine.test(set[,5], method = "fisher")
      p.df.trans_minus[k,3] <- combine.test(set[,5], method = "z.transform")
      p.df.trans_minus[k,4] <- combine.test(set[,5], weight = set[,6], method = "z.transform")
    }
    
    #########
  p.df <- cbind(p.df.trans_plus, p.df.trans_minus)
  p.df <- p.df[,-5]
  colnames(p.df) <- c("Locus_Tag", paste(name.out,"Fisher (trans p with fold change >0)"), paste(name.out,"Z (trans p with fold change >0)"),
				paste(name.out,"Z with weight (trans p with fold change >0)"), paste(name.out,"Fisher (trans p with fold change <0)"), paste(name.out,"Z (trans p with fold change <0)"),
				paste(name.out,"Z with weight (trans p with fold change <0)"))

  anno <- as.data.frame(unique(data[,2]))
  colnames(anno) <- "Locus_Tag"
  anno <- merge(anno, p.df, all.x = TRUE)

  anno <- anno[,c(1,2,5,3,6,4,7)]
  anno1 <- anno
  for (i in 1:nrow(anno)){
	if(!is.na(anno[i,2])){
		minFisher <- min(anno[i,2],anno[i,3])
		minZ <- min(anno[i,4],anno[i,5])
		minZweight <- min(anno[i,6],anno[i,7])
		if(minFisher == anno[i,2]){
			anno[i,3] <- "R"
		}else{
			anno[i,2] <- minFisher
			anno[i,3] <- "L"
		}
		if(minZ == anno[i,4]){
			anno[i,5] <- "R"
		}else{
			anno[i,4] <- minZ
			anno[i,5] <- "L"
		}
		if(minZweight == anno[i,6]){
			anno[i,7] <- "R"
		}else{
			anno[i,6] <- minZweight
			anno[i,7] <- "L"
		}
	}
  }
  
  colnames(anno) <- c("Locus_Tag", paste(name.out,"Fisher"), "Tail", paste(name.out,"Z test"), "Tail", paste(name.out,"Z with weight"), "Tail")

  #### Now we have get the result of input vs output1, input vs output2, output1
  #### vs output2 for every gene

  #### I will sum the raw data and process it to DESeq2.
  data2 <- aggregate(data[,5:ncol(data)], by = list(Locus_Tag = data$Locus_Tag), FUN = sum)
  
  countData <- data2[,2:ncol(data2)]
  countData <- sapply(countData, as.numeric)
  rownames(countData) <- data2[,1]
  countData <- as.data.frame(countData)
  colData <- data.frame(c(rep("A", num.in), rep("B", num.out)))
  rownames(colData) <- colnames(data2)[2:ncol(data2)]
  colnames(colData) <- c("V1")
  
  dds2 <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~V1)
  
  dds2 <- DESeq(dds2)
  ######################
  res2 <- results(dds2)

  output2 <- cbind(data2[,1], res2$log2FoldChange, res2$lfcSE, 
                          res2$pvalue, res2$padj)
  colnames(output2) <- c("Locus_Tag", paste(name.out,"log2FoldChange"), paste(name.out,"lfcSE"), paste(name.out,"pvalue"), paste(name.out,"padj"))
  output2 <- data.frame(output2)
  output2[,2:5] <- sapply(sapply(output2[,2:5], as.character), as.numeric)

  Locus_Tag <- anno$Locus_Tag
  anno <- cbind(Locus_Tag, FC.1, FC.1_1minp, anno[,2:length(anno)])
  names(anno)[2:3] <- c("Input vs Output1 weighted (1/p) based log2FoldChange", "Input vs Output1 weighted adjusted (1/p) based log2FoldChange")
  names(anno)[4:5] <- c("Input vs Output1 weighted (1-p) based log2FoldChange", "Input vs Output1 weighted adjusted (1-p) based log2FoldChange")
  result <- merge(anno, output2, all.x = TRUE)
  

  write.csv(result,paste0(outputname.feature,"_feature", "_in", as.character(num.in), "_out", as.character(num.out), "_P_adj", ".csv")) ## You can replace the "output" with another name!
}








