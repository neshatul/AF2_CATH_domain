


sedb_file <- file.path(CATH_DB_SEQ_PATH, "v4_3_0_seq.fasta")


fl <- list.files(CATH_DB_PATH, pattern = "\\.pdb")

if(!file.exists(sedb_file)){  
  cat("", file=sedb_file, append = F)
  domSeq <- function(i){
    dom <- strsplit(fl[i], split="\\.")[[1]][1]
    po <- bio3d::read.pdb(file.path(CATH_DB_PATH, fl[i])) #pdb object
    seq <- po$atom[po$calpha,]$resid %>% bio3d::aa321() %>% paste0(., collapse = "")
    cat(sprintf(">%s
%s\n", dom, seq), file=sedb_file, append=T)
  }
  
 tmp <- pbmcapply::pbmclapply(1: length(fl), domSeq, mc.cores = NUM_CORES)
}


system(str_interp("${CDHIT_PATH}/cd-hit -i ${sedb_file} -o ${CATH_DB_SEQ_PATH}/v4_3_0_seq_100hit.fasta -c 1.00 -n 5 -M 2000"))
system(str_interp("${CDHIT_PATH}/cd-hit -i ${sedb_file} -o ${CATH_DB_SEQ_PATH}/v4_3_0_seq_099hit.fasta -c 1.00 -n 5 -M 2000"))
system(str_interp("${CDHIT_PATH}/cd-hit -i ${sedb_file} -o ${CATH_DB_SEQ_PATH}/v4_3_0_seq_098hit.fasta -c 1.00 -n 5 -M 2000"))
system(str_interp("${CDHIT_PATH}/cd-hit -i ${sedb_file} -o ${CATH_DB_SEQ_PATH}/v4_3_0_seq_095hit.fasta -c 1.00 -n 5 -M 2000"))
#/hpc/apps/cd-hit/4.8.1
# Run cd-hit on v4_3_0_seq.fasta using following command
# cd-hit -i v4_3_0_seq.fasta -o v4_3_0_seq_100hit.fasta -c 1.00 -n 5 -M 2000
# cd-hit -i v4_3_0_seq.fasta -o v4_3_0_seq_098hit.fasta -c 0.98 -n 5 -M 2000
# cd-hit -i v4_3_0_seq.fasta -o v4_3_0_seq_095hit.fasta -c 0.95 -n 5 -M 2000

#list.files(CATH_DB_SEQ_PATH, pattern = "clstr")
#sequence cluster file
# scf <- c("v4_3_0_seq_095hit.fasta.clstr", 
#          "v4_3_0_seq_098hit.fasta.clstr",
#          "v4_3_0_seq_099hit.fasta.clstr", 
#          "v4_3_0_seq_100hit.fasta.clstr")

scf <- list.files(CATH_DB_SEQ_PATH, pattern = "clstr")

# j=1
# seqClusterFile <- scf[j]                     # cd-hit generated seq cluster file
# CXX <- c("C95", "C98", "C100")[j]
# nCXX <- c("nC95", "nC98", "nC100")[j]
# 
# sc <- readLines(file.path(CATH_DB_SEQ_PATH, seqClusterFile))
# ci <- grep("^>", sc)                                 # cluster index in the file

parseHitClstr <- function(i){
  cb <- ci[i]+1                          #cluster begin line in the cluster file
                                           #cluster end line in the cluster file
  if(i == length(ci)) {
    ce <- length(sc)
  } else{
    ce <- ci[i+1] - 1 
  }
  
  clstr <- sc[cb:ce]
  
  #cluster member
  cm <- lapply(clstr, function(x) strsplit(x, "\\.\\.\\.")[[1]] %>% 
                 .[1] %>% strsplit(. ,">") %>% .[[1]] %>% .[2]) %>% 
        unlist
  
  
  return(cm)
}

#prepare cluster member df
prepClstrMemDf <-  function(i) {
  cmdf <- list()
  cm <- parseHitClstr(i)
  cmdf[["Dom"]] <- cm
  cmdf[[CXX]] <- rep(i, length(cm))
  cmdf[[nCXX]] <- rep(length(cm), length(cm))
  
  return(as.data.frame(cmdf))
  
}

seqClusterFile <- scf[1]                     # cd-hit generated seq cluster file
CXX <- c("C95", "C98", "C99", "C100")[1]
nCXX <- c("nC95", "nC98", "nC99", "nC100")[1]

sc <- readLines(file.path(CATH_DB_SEQ_PATH, seqClusterFile))
ci <- grep("^>", sc) 
Clstr95 <- do.call(rbind, pbmcapply::pbmclapply(1:length(ci), prepClstrMemDf, mc.cores = NUM_CORES))



seqClusterFile <- scf[2]                     # cd-hit generated seq cluster file
CXX <- c("C95", "C98", "C99", "C100")[2]
nCXX <- c("nC95", "nC98", "nC99", "nC100")[2]

sc <- readLines(file.path(CATH_DB_SEQ_PATH, seqClusterFile))
ci <- grep("^>", sc) 
Clstr98 <- do.call(rbind, pbmcapply::pbmclapply(1:length(ci), prepClstrMemDf, mc.cores = NUM_CORES))



seqClusterFile <- scf[3]                     # cd-hit generated seq cluster file
CXX <- c("C95", "C98", "C99", "C100")[3]
nCXX <- c("nC95", "nC98", "nC99", "nC100")[3]

sc <- readLines(file.path(CATH_DB_SEQ_PATH, seqClusterFile))
ci <- grep("^>", sc) 
Clstr99 <- do.call(rbind, pbmcapply::pbmclapply(1:length(ci), prepClstrMemDf, mc.cores = NUM_CORES))



seqClusterFile <- scf[4]                     # cd-hit generated seq cluster file
CXX <- c("C95", "C98", "C99", "C100")[4]
nCXX <- c("nC95", "nC98", "nC99", "nC100")[4]

sc <- readLines(file.path(CATH_DB_SEQ_PATH, seqClusterFile))
ci <- grep("^>", sc) 
Clstr100 <- do.call(rbind, pbmcapply::pbmclapply(1:length(ci), prepClstrMemDf, mc.cores = NUM_CORES))

df1 <- merge(Clstr95, Clstr98, by.y = "Dom")
df2 <- merge(df1, Clstr99, by.y="Dom")
df3 <- merge(df2, Clstr100, by.y="Dom")


clstrFile <- file.path(data, "data05_SeqClust.tab")
write.table(df3, file = clstrFile, row.names = F, col.names = T, quote = F, sep = "\t")








