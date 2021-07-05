suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sys"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("minfi")) 
parser <- ArgumentParser(description='Limma two-group (Unpaired)')

parser$add_argument("-i", "--input", 
                    help="Input  file")
parser$add_argument("-d", "--design", 
                    help="Input design matrix. NOTE: The levels for contrast matrix must be in the second column and be Tumor, Normal.")
parser$add_argument("-o", "--output", 
                    help="output prefix")

args= parser$parse_args()

input<- file_path_as_absolute(args$input)
design<- file_path_as_absolute(args$design)
output<- args$output
setwd(basename(dirname(output)))

#quality Control
clinical <- read.delim(design ,header = TRUE, sep = "\t", row.names = 1)
if (length(unique(clinical[,1])) != 2) {
  stop("There are more than 2 or less than 2 groups in second column of design matrix for comparison.")
}

methyl <- read.delim(input,header = TRUE, sep = "\t",check.names=FALSE)
rowname<- colnames(methyl)[1]
methyl <- methyl%>%column_to_rownames(rowname)
#DE analysis
CheckROws2<- all(row.names(clinical) == colnames(methyl))

if (isTRUE(CheckROws2)){
  ToSubstituate<- factor(clinical[,1], levels = unique(clinical[,1]))
  # Creating design matrix
  design<- model.matrix(~0+ToSubstituate)
  colnames(design)<- c(gsub("ToSubstituate","",colnames(design)[1]),gsub("ToSubstituate","",colnames(design)[2]))
  
  # Fitting linear model and extracting DEGs with FDR<0.05 by topTable  
  fit <- lmFit(methyl, design)
  contrastmatrix<- makeContrasts(TvsN = Tumor-Normal, levels=design)
  contrastfit <- contrasts.fit(fit, contrastmatrix)
  DM<- eBayes(contrastfit)
  Limmaresults<- topTable(DM, number = Inf)
  Limmaresults<- Limmaresults%>%rownames_to_column("Name")
  write.table(Limmaresults, file = paste(output,"_LimmaResults.txt",sep=""),row.names = F, quote = F,sep = "\t")
  print("Job Finished Successfully")  
} else {
  stop("The rows of Design Matrix are not Equal to the columns of data matrix")
}


