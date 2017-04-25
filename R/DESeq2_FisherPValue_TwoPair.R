#' Combine P values for raw data with one input and two outputs
#' 
#' This function is a variation of function \link[mmcc]{DESeq2_FisherPvalue_Onepair}. Now two different types of output are allowed.
#' 
#' This function should be used on a raw data with the following layout:
#' 1. The first three columns must be "Feature", "Locus_Tag", "Strand"
#'    Here "Feture" is like "A000XXX". Be sure that the colnames are right!
#' 2. The next columns are for input, output1, output2 successively. 
#'    Note that there is no constraint for the names of these input and
#'    output columns, but there can only be 2 kinds of outputs. The number
#'    of columns in each output can be any number.
#'    
#' @param filepath  The path of the raw data (excluding the data file itself). Remember to double quote the filepath.
#' @param rawdatafile The name of the raw data file. Remember to double quote the rawdatafile, and add the ".xlsx" after the name of rawdatafile (It must be excel file)!
#' @param outputname.feature The name of output file for each gene you want.
#' @param outputname.oligo The name of output file for each oligo you want.
#' @param num.in The number of input columns in the raw data file.
#' @param num.out1 The number of output1 columns in the raw data file.
#' @param num.out2 The number of output2 columns in the raw data file.
#' @param name.out1 The prefix of output1 columns in the output file. For example, if "A" is input, "B1" is output1, "B2" is output2, then here "B1" can be used.
#' @param name.out2 The prefix of output2 columns in the output file. For example, if "A" is input, "B1" is output1, "B2" is output2, then here "B2" can be used.
#' @keywords mmcc
#' @export
#' @examples DESeq2_FisherPvalue_Twopair(filepath = "C:/Users/", rawdatafile = "rawdata.xlsx", outputname = "output", num.in = 3, num.out1 = 3, num.out2 = 3, name.out1 = "B1", name.out2 = "B2")

DESeq2_FisherPvalue_Twopair <- function(filepath, rawdatafile, outputname.feature, outputname.oligo, num.in, num.out1, num.out2, name.out1, name.out2){
  options(scipen = 999)
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
  colData <- data.frame(c(rep("A", num.in), rep("B", num.out1), rep("C", num.out2)))
  rownames(colData) <- colnames(data)[5:length(data)]
  colnames(colData) <- c("V1")
  
  dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~V1)
  
  dds <- DESeq(dds)
  ### for output1
  res.out1 <- results(dds, contrast = c("V1", "B", "A"))
  output.out1 <- cbind(data[,2], res.out1$log2FoldChange, res.out1$lfcSE, 
                          res.out1$pvalue, res.out1$padj)
  colnames(output.out1) <- c("Locus_Tag", "log2FoldChange", "lfcSE", "pvalue", "padj")
  output.out1 <- data.frame(output.out1)
  output.out1[,2:5] <- sapply(sapply(output.out1[,2:5], as.character), as.numeric)
  output.out1$weight <- 1/(output.out1$lfcSE^2)
  
  ### for output2
  res.out2 <- results(dds, contrast = c("V1", "C", "A"))
  output.out2 <- cbind(data[,2], res.out2$log2FoldChange, res.out2$lfcSE, 
                            res.out2$pvalue, res.out2$padj)
  colnames(output.out2) <- c("Locus_Tag", "log2FoldChange", "lfcSE", "pvalue", "padj")
  output.out2 <- data.frame(output.out2)
  output.out2[,2:5] <- sapply(sapply(output.out2[,2:5], as.character), as.numeric)
  output.out2$weight <- 1/(output.out2$lfcSE^2)
  
  ### for output1 vs output2
  res.out3 <- results(dds, contrast = c("V1", "B", "C"))
  output.out3 <- cbind(data[,2], res.out3$log2FoldChange, res.out3$lfcSE, 
                            res.out3$pvalue, res.out3$padj)
  colnames(output.out3) <- c("Locus_Tag", "log2FoldChange", "lfcSE", "pvalue", "padj")
  output.out3 <- data.frame(output.out3)
  output.out3[,2:5] <- sapply(sapply(output.out3[,2:5], as.character), as.numeric)
  output.out3$weight <- 1/(output.out3$lfcSE^2)
  
  oligo <- cbind(output.out1[])
  alldata <- list(data1 = output.out1, data2 = output.out2,
                  data3 = output.out3)
  
  p.df <- list()

  for (i in 1:3){
    df <- alldata[[i]]
    df[,4] <- df[,4] / 2
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
    p.df[[i]] <- cbind(p.df.trans_plus, p.df.trans_minus)
    p.df[[i]] <- p.df[[i]][,-5]
    if (i == 1){
    	colnames(p.df[[i]]) <- c("Locus_Tag", paste(name.out1,"Fisher (trans p with fold change >0)"), paste(name.out1,"Z (trans p with fold change >0)"),
				paste(name.out1,"Z with weight (trans p with fold change >0)"), paste(name.out1,"Fisher (trans p with fold change <0)"), paste(name.out1,"Z (trans p with fold change <0)"),
				paste(name.out1,"Z with weight (trans p with fold change <0)"))
    }else if (i == 2){
	colnames(p.df[[i]]) <- c("Locus_Tag", paste(name.out2,"Fisher (trans p with fold change >0)"), paste(name.out2,"Z (trans p with fold change >0)"),
				paste(name.out2,"Z with weight (trans p with fold change >0)"), paste(name.out2,"Fisher (trans p with fold change <0)"), paste(name.out2,"Z (trans p with fold change <0)"),
				paste(name.out2,"Z with weight (trans p with fold change <0)"))
    }else if (i == 3){
	colnames(p.df[[i]]) <- c("Locus_Tag", paste(name.out1,"vs",name.out2,"Fisher (trans p with fold change >0)"), paste(name.out1,"vs",name.out2,"Z (trans p with fold change >0)"),
				paste(name.out1,"vs",name.out2,"Z with weight (trans p with fold change >0)"), paste(name.out1,"vs",name.out2,"Fisher (trans p with fold change <0)"), paste(name.out1,"vs",name.out2,"Z (trans p with fold change <0)"),
				paste(name.out1,"vs",name.out2,"Z with weight (trans p with fold change <0)"))
    }
  }
  anno <- as.data.frame(unique(data[,2]))
  colnames(anno) <- "Locus_Tag"
  for (i in 1:3){
    anno <- merge(anno, p.df[[i]], all.x = TRUE)
  }
  
  anno <- anno[,c(1,2,5,3,6,4,7,8,11,9,12,10,13,14,17,15,18,16,19)]
  for (i in 1:nrow(anno)){
    if(!is.na(anno[i,2])){
      minFisherOut1 <- min(anno[i,2], anno[i,3])
      minZOut1 <- min(anno[i,4], anno[i,5])
      minZweightOut1 <- min(anno[i,6], anno[i,7])
      
      if(minFisherOut1 == anno[i,2]){
        anno[i,3] <- "R"
      }else{
        anno[i,2] <- minFisherOut1
        anno[i,3] <- "L"
      }
      if(minZOut1 == anno[i,4]){
        anno[i,5] <- "R"
      }else{
        anno[i,4] <- minZOut1
        anno[i,5] <- "L"
      }
      if(minZweightOut1 == anno[i,6]){
        anno[i,7] <- "R"
      }else{
        anno[i,6] <- minZweightOut1
        anno[i,7] <- "L"
      }
    }
    if (!is.na(anno[i,8])){
      minFisherOut2 <- min(anno[i,8], anno[i,9])
      minZOut2 <- min(anno[i,10], anno[i,11])
      minZweightOut2 <- min(anno[i,12], anno[i,13])

      if(minFisherOut2 == anno[i,8]){
        anno[i,9] <- "R"
      }else{
        anno[i,8] <- minFisherOut2
        anno[i,9] <- "L"
      }
      if(minZOut2 == anno[i,10]){
        anno[i,11] <- "R"
      }else{
        anno[i,10] <- minZOut2
        anno[i,11] <- "L"
      }
      if(minZweightOut2 == anno[i,12]){
        anno[i,13] <- "R"
      }else{
        anno[i,12] <- minZweightOut2
        anno[i,13] <- "L"
      }
    }
    if(!is.na(anno[i,14])){
      minFisherOut1Out2 <- min(anno[i,14], anno[i,15])
      minZOut1Out2 <- min(anno[i,16], anno[i,17])
      minZweightOut1Out2 <- min(anno[i,18], anno[i,19])
      
      if(minFisherOut1Out2 == anno[i,14]){
        anno[i,15] <- "R"
      }else{
        anno[i,14] <- minFisherOut1Out2
        anno[i,15] <- "L"
      }
      if(minZOut1Out2 == anno[i,16]){
        anno[i,17] <- "R"
      }else{
        anno[i,16] <- minZOut1Out2
        anno[i,17] <- "L"
      }
      if(minZweightOut1Out2 == anno[i,18]){
        anno[i,19] <- "R"
      }else{
        anno[i,18] <- minZweightOut1Out2
        anno[i,19] <- "L"
      }
    }
  }
  
  colnames(anno) <- c("Locus_Tag", paste(name.out1,"Fisher"), "Tail", paste(name.out1,"Z test"), "Tail", paste(name.out1,"Z with weight"), "Tail"
                      , paste(name.out2,"Fisher"), "Tail", paste(name.out2,"Z test"), "Tail", paste(name.out2,"Z with weight"), "Tail"
                      , paste(name.out1,"vs",name.out2,"Fisher"), "Tail", paste(name.out1,"vs",name.out2,"Z test"), "Tail", paste(name.out1,"vs",name.out2,"Z with weight"), "Tail")
  
  #### Now we have get the result of input vs output1, input vs output2, output1
  #### vs output2 for every gene

  #### I will sum the raw data and process it to DESeq2.
  data2 <- aggregate(data[,5:ncol(data)], by = list(Locus_Tag = data$Locus_Tag), FUN = sum)
  
  countData <- data2[,2:ncol(data2)]
  countData <- sapply(countData, as.numeric)
  rownames(countData) <- data2[,1]
  countData <- as.data.frame(countData)
  colData <- data.frame(c(rep("A", num.in), rep("B", num.out1), rep("C", num.out2)))
  rownames(colData) <- colnames(data2)[2:ncol(data2)]
  colnames(colData) <- c("V1")
  
  dds2 <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData,
                                design = ~V1)
  
  dds2 <- DESeq(dds2)
  ######################
  res2.ab <- results(dds2, contrast = c("V1", "B", "A"))
  res2.ac <- results(dds2, contrast = c("V1", "C", "A"))
  res2.bc <- results(dds2, contrast = c("V1", "B", "C"))

  output2.out1 <- cbind(data2[,1], res2.ab$log2FoldChange, res2.ab$lfcSE, 
                          res2.ab$pvalue, res2.ab$padj)
  colnames(output2.out1) <- c("Locus_Tag", paste(name.out1,"log2FoldChange"), paste(name.out1,"lfcSE"), paste(name.out1,"pvalue"), paste(name.out1,"padj"))
  output2.out1 <- data.frame(output2.out1)
  output2.out1[,2:5] <- sapply(sapply(output2.out1[,2:5], as.character), as.numeric)

  output2.out2 <- cbind(data2[,1], res2.ac$log2FoldChange, res2.ac$lfcSE, 
                          res2.ac$pvalue, res2.ac$padj)
  colnames(output2.out2) <- c("Locus_Tag", paste(name.out2,"log2FoldChange"), paste(name.out2,"lfcSE"), paste(name.out2,"pvalue"), paste(name.out2,"padj"))
  output2.out2 <- data.frame(output2.out2)
  output2.out2[,2:5] <- sapply(sapply(output2.out2[,2:5], as.character), as.numeric)

  output2.out3 <- cbind(data2[,1], res2.bc$log2FoldChange, res2.bc$lfcSE, 
                          res2.bc$pvalue, res2.bc$padj)
  colnames(output2.out3) <- c("Locus_Tag", paste(name.out1, "vs", name.out2,"log2FoldChange"),  paste(name.out1, "vs", name.out2,"lfcSE"),  paste(name.out1, "vs", name.out2,"pvalue"),  paste(name.out1, "vs", name.out2,"padj"))
  output2.out3 <- data.frame(output2.out3)
  output2.out3[,2:5] <- sapply(sapply(output2.out3[,2:5], as.character), as.numeric)
  
  result.presum <- cbind(output2.out1, output2.out2, output2.out3)
  result.presum <- result.presum[,c(1:5,7:10,12:15)]

  result <- merge(anno, result.presum, all.x = TRUE)

  write.csv(result,paste0(outputname.feature,".csv"))
}

