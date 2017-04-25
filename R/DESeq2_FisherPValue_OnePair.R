#' Combine P values for one pair of input and output
#' 
#' This function is used to: 
#' 1. To the raw data, using DESeq2 to give the p values of input vesus output
#' 2. Combine the p values in the output of DESeq2
#' 3. To the raw data, pre-sum the count data group by locus_tag, and using DESeq2 to give the p values of input vesus output
#' 
#' This function should be used on a raw data with the following layout:
#' 1. The first three columns must be "Feature", "Locus_Tag", "Strand"
#'    Here "Feture" is like "A000XXX". Be sure that the colnames are right!
#' 2. The next columns are for input and output successively. 
#'    Note that there is no constraint for the names of these input and
#'    output columns, but there can only be 1 kind of output. The number
#'    of columns in each output can be any number.
#' 
#' @param filepath  The path of the raw data (excluding the data file itself). Remember to double quote the filepath.
#' @param rawdatafile The name of the raw data file. Remember to double quote the rawdatafile, and add the ".xlsx" after the name of rawdatafile (It must be excel file)!
#' @param num.in The number of input columns in the raw data file.
#' @param num.out The number of output columns in the raw data file.
#' @param name.out The prefix of columns in the output. For example, "A vs B" or "B vs C", in which "A". "B" and "C" represent the names of input and output in the raw data file.
#' @keywords mmcc
#' @export
#' @examples DESeq2_FisherPvalue_Onepair(filepath = "C:/Users/", rawdatafile = "rawdata.xlsx", outputname = "output", num.in = 3, num.out = 3, name.out = "A vs B")

DESeq2_FisherPvalue_Onepair <- function(filepath, rawdatafile, 
                                        num.in, num.out, name.out = "Input vs Output"){
  options(scipen = 999)
  outputname.feature <- substr(rawdatafile, 1, gregexpr("\\.xlsx", rawdatafile)[[1]][1] - 1)
  outputname.oligo <- substr(rawdatafile, 1, gregexpr("\\.xlsx", rawdatafile)[[1]][1] - 1)
  setwd(filepath)
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
  
  write.csv(output, paste0(outputname.oligo, "_in", as.character(num.in), "_out", as.character(num.out), ".csv"))
  
  df <- output

  df$pvalue <- df$pvalue / 2
  df.trans_plus <- df
  df.trans_minus <- df
  for (j in 1:nrow(df)){
    if (!is.na(df[j,2])){
      if (df[j,2] > 0 & !is.na(df[j,4])){
        df.trans_plus[j, 4] <- 1 - df.trans_plus[j, 4]
      }else if (df[j,2] < 0 & !is.na(df[j,4])){
        df.trans_minus[j, 4] <- 1 - df.trans_minus[j, 4]
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
      p.df.trans_plus[k,2] <- combine.test(set[,4], method = "fisher")
      p.df.trans_plus[k,3] <- combine.test(set[,4], method = "z.transform")
      p.df.trans_plus[k,4] <- combine.test(set[,4], weight = set[,6], method = "z.transform")
    }
    
    #########

    p.df.trans_minus <- data.frame()
    
    for (k in 1:length(df2.trans_minus)){
      set <- df2.trans_minus[[k]]
      p.df.trans_minus[k,1] <- set[1,1]
	set <- na.omit(set)
      p.df.trans_minus[k,2] <- combine.test(set[,4], method = "fisher")
      p.df.trans_minus[k,3] <- combine.test(set[,4], method = "z.transform")
      p.df.trans_minus[k,4] <- combine.test(set[,4], weight = set[,6], method = "z.transform")
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

  result <- merge(anno, output2, all.x = TRUE)

  write.csv(result,paste0(outputname.feature,"_in", as.character(num.in), "_out", as.character(num.out), ".csv")) ## You can replace the "output" with another name!
  }

