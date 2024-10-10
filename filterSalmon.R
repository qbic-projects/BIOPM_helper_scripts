rm(list = ls(all = TRUE)); graphics.off()
path <- getwd(); setwd(path)
path

args = commandArgs(trailingOnly=TRUE)

#args[1] <- "star_salmon"
#args[2] <- "final_DE_gene_list.tsv"
#args[3] <- 5
#args[4] <- "metadata.tsv"

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} else if (length(args)==1) {
}

system("rm -rf results")
dir.create("results")
dir.create("results/input")
dir.create("results/output")
dir.create(paste("results/output",args[1],sep = "/"))

# prepare and copy to input
system(paste("cp -r",args[1], "results/input/",sep = " "))
sf <- list.files("results/input",recursive = T,full.names = T)
sf <- sf[grepl("quant",sf)]
sf <- sf[!grepl("logs",sf)]
file.copy(args[2], "results/input/")
fn <- "results/input/baseMeanValue.txt"
sink(fn)
print(paste("baseMean value used for filtering:",as.numeric(args[3]),sep = " "))
sink()
file.copy(args[4], "results/input/")


#prepare DE list to filter
de <- read.table(args[2],
                 header = T,sep = "\t",na.strings =c(""," ","NaN"),
                 quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
de$baseMean <- as.numeric(de$baseMean)


de <- subset(de, outcome == "DE")
de <- subset(de, baseMean < as.numeric(args[3]))

if (nrow(de) == 0) {
  print("No DEGs with baseMean < 5, stop")
}else{
  print(paste(nrow(de),"DEGs found with baseMean < 5, continue...",sep = " "))
  Sys.sleep(3)
}
#as salmon folders has Ensembl_ID, work with those...
de <- unique(de$Ensembl_ID)
#write to file
write.table(de, file = "results/output/filtered_Ensembl_IDs.txt" , append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = T,  
            col.names = F, qmethod = c("escape", "double"))

#metadata
m <- read.table(args[4],
                 header = T,sep = "\t",na.strings =c(""," ","NaN"),
                 quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
codes <- m$QBiC.Code
if (length(codes) == 0) {
  print("No QBiC Ids, stop")
}else{
  print(paste(length(codes),  " QBiC Ids found, continue...",sep =""))
  Sys.sleep(3)    
}


for (i in 1:length(codes)) {
  #i=1
  fn <- codes[i]
  print(paste("QBiC Barcode:",fn,sep = " "))
  dir.create(paste("results/output",args[1],fn,sep = "/"))
  wd <- paste("results/output",args[1],fn,sep = "/")
  fn <- sf[grepl(fn,sf)]
  if (length(fn) == 0) {
    print(paste("No salmon input files with barcode ", codes[i],sep=""))
  }else{
    print(paste("Process now files ",fn,sep = ""))
    print("####")
    fn <- unlist(strsplit(fn," "))
    for (f in 1:length(fn)) {
      #f=1
      tmp <- read.table(fn[f],
                        header = T,sep = "\t",na.strings =c(""," ","NaN"),
                        quote=NULL,stringsAsFactors=F,dec=".",fill=TRUE)
      tmp <- subset(tmp,!Name %in% de)
      nn <- basename(fn[f])
      write.table(tmp, file = paste(wd,nn,sep="/"), append = FALSE, quote = FALSE, sep = "\t",
                  eol = "\n", na = "NA", dec = ".", row.names = T,  
                  col.names = NA, qmethod = c("escape", "double"))
      rm(tmp,nn)
  }
  }
  rm(fn)
}

print(paste("Salmon Folder with filtered files to be used for resubmitting with qbic-pipelines/rnadeseq can be found here: ",
            dirname(wd),sep=""))
#clean-up










