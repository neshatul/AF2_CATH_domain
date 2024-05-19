cath_str_db_path <- "/scratch/g/mtzimmermann/data/Cath_Local_db/v4_3_0"
cth_af2_str_db_path <- "/scratch/g/mtzimmermann/data/Cath_Local_db/v4_3_0_af2_dom_aliOk"
af2_str_db_path <- "/scratch/g/mtzimmermann/data/AF2_local_db/af2_general"


kras_unip <- 'P01116'
# cath enties for K-Ras
cth_kras <- df1$Dom[df1$SP_PRIMARY == kras_unip]





library(Biostrings)
library(ShortRead)
seq_list <- c()
for(j in seq_along(cth_kras)){
  #j <- which(row.names(KRAS_RMSD1) == '3gftC00')
  
  cthdom <- cth_kras[j]
  cth_p <- read.pdb(file.path(cath_str_db_path, paste0(cthdom, '.pdb')))
  cth_seq <- aa31( cth_p$atom[cth_p$calpha,]$resid) %>% paste0(., collapse = '')
  
  af2_p <- read.pdb(file.path(cth_af2_str_db_path, paste0(cthdom, '_af2.pdb')))
  af2_seq <- aa31( af2_p$atom[af2_p$calpha,]$resid) %>% paste0(., collapse = '')
  
  names(cth_seq) <- cthdom
  seq_list <- c(seq_list, cth_seq)
  
}


# KRAS full sequence
kras_af2_FN <- paste0(kras_unip, "-F1-model_v3.pdb")
kras_p <- read.pdb(file.path(af2_str_db_path, kras_af2_FN))
Kras_full_seq <- aa31( kras_p$atom[kras_p$calpha,]$resid) %>% paste0(., collapse = '')
names(Kras_full_seq) <- 'KRAS'

kras_allseq <- AAStringSet(c(Kras_full_seq, seq_list))
if(0){
  writeXStringSet(kras_allseq, 'kras.fasta')
}

require(msa)
# kras_msa <- msa(kras_allseq, method="Muscle", gapOpening=12, gapExtension=3, maxiters=16,
#     cluster="upgmamax", SUEFF=0.4, brenner=FALSE,
#     order="input", verbose=FALSE)
kras_msa <- msa(kras_allseq, method="ClustalOmega",
                order="input", verbose=FALSE)

mymsa_file <- 'kras.alan.fasta'
kras_msa1 <- readAAStringSet(mymsa_file)



# The ref requence analysis
# here ref Seqis K-Ras
# If the msa is among the sequence of same protein then gaps in the middle of the seq is not accepted

# CHECK for gaps in the middle
# following code is incomplete
if(0){
  require(dplyr)
  RefSeq_chr <- as.character(kras_msa@unmasked[1]) %>% strsplit(.,'') %>% .[[1]]
  # RefSeq_chr <- c("-", "-", "-", "-", "M", "T", "E", "Y", "K", "L", "V", "V", 
  #                 "V", "G", "A", "-", "-", "V", "G", "K", "S", "A", "L", "T", "I", 
  #                 "Q", "L", "I", "Q", "N", "H", "F", "-", "D", "E", "Y", "-", "-")
  Chechk4Gaps <- function(RefSeq_chr=NULL){
    gaps_inds <- which(RefSeq_chr=='-')
    
    # Check if gaps are in begining or at end
    gapsRinBegEnd <- c()
    for(i in gaps_inds){
      if(i==1){
        gapsRinBegEnd <- c(gapsRinBegEnd, 'TRUE')
      }else if(i == length(RefSeq_chr)){
        gapsRinBegEnd <- c(gapsRinBegEnd, 'TRUE')
      } else if(i != 1){
        backward_check <- RefSeq_chr[i-1] == '-'
        forward_check <- RefSeq_chr[i+1] == '-'
        gapsRinBegEnd <- c(gapsRinBegEnd, all(backward_check, forward_check))
      }
    }
    return(gapsRinBegEnd)
  }
  
  if(all(Chechk4Gaps(RefSeq_chr))){
    stop('The RefSequence has gap in the middle')
  }
}
################################################################################

RefSeq_chr <- as.character(kras_msa@unmasked[1]) %>% strsplit(.,'') %>% .[[1]]
gaps_inds <- which(RefSeq_chr=='-')
ref_resid_offset <- max(gaps_inds)

dfw_kras1 <- list()
fl <- 1:164 # KRAS 'f'ull 'l'ength  residue list
s1 <- 26:46 # KRAS switch I residue list
s2 <- 51:78 # KRAS switch II residue list


as.character(kras_msa@unmasked[2]) %>% strsplit(.,'') %>% .[[1]]


xr <- as.character(kras_msa@unmasked[2]) %>% strsplit(.,'') %>% .[[1]]
xr2 <- which(xr!='-')

x1 <- as.character(kras_msa@unmasked[4]) %>% strsplit(.,'') %>% .[[1]]
x2 <- which(x1!='-')

Correct_order_of_res <- list()

for(i in 2:length(kras_msa@unmasked)){
  #i <- which(names(kras_msa@unmasked) == '6bofA00')
  ref_resid_offset <- 4
  Seq_ref <- as.character(kras_msa@unmasked[1]) %>% strsplit(.,'') %>% .[[1]]
  Seq_ref_ind <- which(Seq_ref!='-') - ref_resid_offset
  
  Seq_X <- as.character(kras_msa@unmasked[i]) %>% strsplit(.,'') %>% .[[1]]
  Seq_X_ind <- which(Seq_X!='-') - ref_resid_offset
  
  cthdom <- names(kras_msa@unmasked)[i]
  cthdom_file <- file.path(cath_str_db_path, paste0(cthdom, ".pdb"))
  cthdom_pdb <- read.pdb(cthdom_file, verbose = F)
  
  Corr_order <- all(cthdom_pdb$atom[cthdom_pdb$calpha,]$resno == Seq_X_ind)
  inCorrResOrder_count <- sum(cthdom_pdb$atom[cthdom_pdb$calpha,]$resno != Seq_X_ind)
  Correct_order_of_res$Dom[i] <- cthdom
  Correct_order_of_res$CorrResOrder[i] <- Corr_order
  Correct_order_of_res$inCorrResOrder_count[i] <- inCorrResOrder_count
  
}
Correct_order_of_res <- data.frame(Correct_order_of_res)
Correct_order_of_res <- Correct_order_of_res[-1,] # Remove position of REF from '1'
message('K-Ras structures from CATH DB not having correct order of residues')
Correct_order_of_res[Correct_order_of_res$CorrResOrder == FALSE,]
dim(Correct_order_of_res[Correct_order_of_res$CorrResOrder == FALSE,])

################################## MSA from server #############################
Correct_order_of_res <- list()
for(i in 2:length(kras_msa1)){
  #i <- which(names(kras_msa@unmasked) == '5tb5A00')
  ref_resid_offset <- 4
  Seq_ref <- as.character(kras_msa1[1]) %>% strsplit(.,'') %>% .[[1]]
  Seq_ref_ind <- which(Seq_ref!='-') - ref_resid_offset
  
  Seq_X <- as.character(kras_msa1[i]) %>% strsplit(.,'') %>% .[[1]]
  Seq_X_ind <- which(Seq_X!='-') - ref_resid_offset
  
  cthdom <- names(kras_msa1)[i]
  cthdom_file <- file.path(cath_str_db_path, paste0(cthdom, ".pdb"))
  cthdom_pdb <- read.pdb(cthdom_file, verbose = F)
  
  corr_order <- all(cthdom_pdb$atom[cthdom_pdb$calpha,]$resno == Seq_X_ind)
  inCorrResOrder_count <- sum(cthdom_pdb$atom[cthdom_pdb$calpha,]$resno != Seq_X_ind)
  Correct_order_of_res$Dom[i] <- cthdom
  Correct_order_of_res$CorrResOrder[i] <- corr_order
  Correct_order_of_res$inCorrResOrder_count[i] <- inCorrResOrder_count
  
}
Correct_order_of_res <- data.frame(Correct_order_of_res)
Correct_order_of_res <- Correct_order_of_res[-1,] # Remove position of REF from '1'
message('K-Ras structures from CATH DB not having correct order of residues')
Correct_order_of_res[Correct_order_of_res$CorrResOrder == FALSE,]
dim(Correct_order_of_res[Correct_order_of_res$CorrResOrder == FALSE,])
################################################################################

################################################################################

for(i in 2:length(kras_msa@unmasked)){
  #i <- which(names(kras_msa@unmasked) == '6gqtD00')
  ref_resid_offset <- 4
  Seq_ref <- as.character(kras_msa@unmasked[1]) %>% strsplit(.,'') %>% .[[1]]
  Seq_ref_ind <- which(Seq_ref!='-') - ref_resid_offset
  
  Seq_X <- as.character(kras_msa@unmasked[i]) %>% strsplit(.,'') %>% .[[1]]
  Seq_X_ind <- which(Seq_X!='-') - ref_resid_offset
  
  cthdom <- names(kras_msa@unmasked)[i]
  cthdom_file <- file.path(cath_str_db_path, paste0(cthdom, ".pdb"))
  cthdom_pdb <- read.pdb(cthdom_file, verbose = F)
  
  cthdom_pdb$atom[cthdom_pdb$calpha,]$resno == Seq_X_ind
}

