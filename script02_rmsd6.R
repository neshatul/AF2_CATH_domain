#'@description The script extracts domains from AF2 structures. 
#'RMSD6 
#'pairwisealignment with local-global produced satisfactory result
#'terminal truncation done using isCore


isCore <- function(pdbobj=NULL, nResCheck=10, distCutoff=4.5){
  
  pdbobj$atom$eleno <- 1:dim(pdbobj$atom)[1]
  rl <- pdbobj$atom[pdbobj$calpha,]$resno # residue list
  
  # N-ter residue list
  nt_resL <- rl[1:nResCheck]
  nt_atmL_cols <- pdbobj$atom[pdbobj$atom$resno %in% nt_resL,]$eleno
  nt_atmL_rows <- pdbobj$atom[!pdbobj$atom$resno %in% nt_resL,]$eleno
  ct_resL <- rl[(length(rl) - nResCheck + 1): length(rl)]
  ct_atmL_cols <- pdbobj$atom[pdbobj$atom$resno %in% ct_resL,]$eleno
  ct_atmL_rows <- pdbobj$atom[!pdbobj$atom$resno %in% ct_resL,]$eleno
  
  
  cth_dist <- pdbobj$atom[c('x', 'y', 'z')] %>% dist(., method="euclidean", upper=F, diag=F)
  cth_dist <- as.matrix(cth_dist)
  
  nt_dist_logical <- cth_dist[nt_atmL_rows, nt_atmL_cols] <=distCutoff
  nt_atm_int_list <- apply(nt_dist_logical, 2, function(x){any(x==TRUE)}) 
  
  nt_atm_TRUE <- nt_atm_int_list[nt_atm_int_list==TRUE]
  ntTRUE <- pdbobj$atom[pdbobj$atom$eleno %in% as.numeric(names(nt_atm_TRUE)),]$resno %>% unique()
  if(length(ntTRUE)>4){isNTerCore <- TRUE} else {isNTerCore <- FALSE}
  
  ct_dist_logical <- cth_dist[ct_atmL_rows, ct_atmL_cols] <=distCutoff
  ct_atm_int_list <- apply(ct_dist_logical, 2, function(x){any(x==TRUE)}) 
  
  ct_atm_TRUE <- ct_atm_int_list[ct_atm_int_list==TRUE]
  ctTRUE <- pdbobj$atom[pdbobj$atom$eleno %in% as.numeric(names(ct_atm_TRUE)),]$resno %>% unique()
  if(length(ctTRUE)>4){isCTerCore <- TRUE} else {isCTerCore <- FALSE}
  ter_assess <- list()
  ter_assess$isNTerCore <- isNTerCore
  ter_assess$isCTerCore <- isCTerCore
  return(ter_assess)
}
#isCore(bio3d::read.pdb("/scratch/g/mtzimmermann/data/Cath_Local_db/v4_3_0/1e18A01.pdb"))
#d=0

GetAF2Cols <- function(d, outputdir, outputfile, 
                       SAVE_TRUNCATED=FALSE, 
                       EXTRACT_AF2_DOM=FALSE){
  
  cthdom <- CATHUNP2$Dom[d]
  cthdom_file <- file.path(CATH_DB_PATH, paste0(cthdom, ".pdb"))
  
  error01 <- FALSE
  tryCatch({
    cthdom_pdb <- read.pdb(cthdom_file)
    if(dim(cthdom_pdb$atom[cthdom_pdb$calpha,])[1] > 20){
      # Number of residues to check to be logical part of the domain
      nResCheck <- 10
      rl <- cthdom_pdb$atom[cthdom_pdb$calpha,]$resno # residue list
      
      # N-ter residue list
      nt_resL <- rl[1:nResCheck]
      ct_resL <- rl[(length(rl) - nResCheck + 1): length(rl)]
      
      
      #No of N-ter & C-ter residue to check
      tercheck_NC <- isCore(cthdom_pdb)
      tercheck_N <- tercheck_NC$isNTerCore
      tercheck_C <- tercheck_NC$isCTerCore
      
      if(!tercheck_N){
        atmvector <- cthdom_pdb$atom[!cthdom_pdb$atom$resno %in% nt_resL,]$eleno
        cthres_inds <- atom.select(pdb=cthdom_pdb, eleno = atmvector)
        NterChecked_pdb <- trim.pdb(cthdom_pdb, cthres_inds)
        Nter_trunc <- TRUE
      } else{
        NterChecked_pdb <- cthdom_pdb
        Nter_trunc <- FALSE
      }
      
      if(!tercheck_C){
        atmvector <- NterChecked_pdb$atom[!NterChecked_pdb$atom$resno %in% ct_resL,]$eleno
        cthres_inds <- atom.select(pdb=NterChecked_pdb, eleno = atmvector)
        CterChecked_pdb <- trim.pdb(NterChecked_pdb, cthres_inds)
        Cter_trunc <- TRUE
      } else{
        CterChecked_pdb <- NterChecked_pdb
        Cter_trunc <- FALSE
      }
      
      cthdom_pdb <- CterChecked_pdb
      
      
      AF2df=list()
      cthatom_vector <- cthdom_pdb$atom[cthdom_pdb$calpha,]$eleno
      cthdom_seq <- cthdom_pdb$atom[cthdom_pdb$calpha,]$resid %>% aa321() %>% paste0(., collapse = "")
      cthdom_seq_length <- cthdom_pdb$atom[cthdom_pdb$calpha,]$resid %>% aa321() %>% length()
      
      cthres_inds <- atom.select(pdb=cthdom_pdb, eleno = cthatom_vector)
      cthCaPDB <- trim.pdb(cthdom_pdb, cthres_inds)
      
      uid <- CATHUNP2$SP_PRIMARY[d]
      af2_fn <- paste0(uid, "-F1-model_v3.pdb")
      af2_file <- file.path(AF2_DB_PATH, af2_fn)
      af2_pdb <- read.pdb(af2_file)
      af2_seq <- af2_pdb$atom[af2_pdb$calpha,]$resid %>% aa321() %>% paste0(., collapse = "")
      
      seqlist <- c(af2_seq, cthdom_seq)
      names(seqlist) <- c("af2", "cth")
      
      #seq AAstring set
      seq_ss <- AAStringSet(seqlist)
      
      #pairwise alignment
      pali <- pairwiseAlignment(subject = cthdom_seq, pattern = af2_seq, substitutionMatrix = "BLOSUM62", 
                                gapOpening = 10, gapExtension = 0.2, 
                                type = "local-global",
                                #type = "global",
                                scoreOnly = FALSE)
      #writePairwiseAlignments(pali, block.width = 90)
      
      #z1 <- msa(seq_ss, gapOpening = 10, gapExtension = 0.2, method = "ClustalW", type = "protein")
      
      
      seq1aln <- pattern(pali) # Get the alignment for the first sequence
      seq2aln <- subject(pali) # Get the alignment for the second sequence
      alnlen  <- nchar(seq1aln)     # Find the number of columns in the alignment
      seq1aln_st <- seq1aln@range@start
      
      af2_seq_vector_incorrect_resno <- seq2aln %>% as.character(.) %>%  strsplit(., split='') %>% .[[1]] %>% grep(pattern = '[A-Z]', x = .) # Initial gaps not considered
      af2_seq_vector <- af2_seq_vector_incorrect_resno  + seq1aln_st-1
      
      
      # percent of residues exactly matched in alignment out of total cath dom res no
      # In general PID , percent identity, means  percentage of identity out of identity + matched + gapped
      PID1 <- (nmatch(pali) * (100/cthdom_seq_length))
      PID <- pid(pali)
      
      
      af2atom_vector <- af2_pdb$atom[af2_pdb$calpha,][af2_seq_vector,]$eleno
      af2_ca_inds <- atom.select(pdb=af2_pdb, eleno =af2atom_vector)
      af2CaPDB <- trim.pdb(af2_pdb, af2_ca_inds)
      
      
      cth_af2_rmsd <- rmsd(cthCaPDB$xyz, af2CaPDB$xyz, fit = T)
      if(EXTRACT_AF2_DOM){
        cth_af2_str_fn <- paste0(cthdom, "_af2.pdb")
        if(!file.exists(file.path(AF2_DOM_DB_PATH, cth_af2_str_fn))){
          af2_resno_inds <- atom.select(pdb=af2_pdb, resno =af2_seq_vector)
          af2resPDB <- trim.pdb(af2_pdb, af2_resno_inds)
          
          write.pdb(af2resPDB, file.path(AF2_DOM_DB_PATH, cth_af2_str_fn))
          
        }
      }
      
      if(SAVE_TRUNCATED){
        # writting the cathdom pdb after terminal truncation check
        cthdom_truncFN <- file.path(CATH_DB_PATH_TRUNC, paste0(cthdom, ".pdb"))
        if(!file.exists(cthdom_truncFN)){
          write.pdb(cthdom_pdb, file.path(cthdom_truncFN))
        }
      }
      
      #sse <- stride.WF(cthdom_file)
      
      AF2df$Dom[1] <- CATHUNP2$Dom[d]
      AF2df$SP_PRIMARY[1] <- CATHUNP2$SP_PRIMARY[d]
      
      plddt_vector <- af2CaPDB$atom[af2CaPDB$calpha, ]$b
      plddt <- quantile(plddt_vector, c(.1, .2, .3))
      
      HSeg <- Vector2HyphenSeg(af2_seq_vector)
      # if(length(strsplit(HSeg, split=",")[[1]]) < 5){
      #   AF2df$af2ResBoundary[d] <- Vector2HyphenSeg(af2_seq_vector)
      # } else  AF2df$af2ResBoundary[d] <- "TooManySegs"
      
      AF2df$af2ResBoundary[1] <- HSeg
      AF2df$af2DomLen[1] <- length(af2_seq_vector)
      AF2df$PID[1] <- round(PID, 3)
      AF2df$PID1[1] <- round(PID1, 3)
      AF2df$pLDDTq10[1] <- plddt[[1]]
      AF2df$pLDDTq20[1] <- plddt[[2]]
      AF2df$pLDDTq30[1] <- plddt[[3]]
      AF2df$RMSD6[1] <- cth_af2_rmsd
      AF2df$Nter_trunc[1] <- Nter_trunc
      AF2df$Cter_trunc[1] <- Cter_trunc
      
    }else {
      AF2df$Dom[1] <-  CATHUNP2$Dom[d]
      AF2df$SP_PRIMARY[1] <- CATHUNP2$SP_PRIMARY[d]
      AF2df$af2ResBoundary[1] <- NA
      AF2df$af2DomLen[1] <- NA
      AF2df$PID[1] <- NA
      AF2df$PID1[1] <- NA
      AF2df$pLDDTq10[1] <- NA
      AF2df$pLDDTq20[1] <- NA
      AF2df$pLDDTq30[1] <- NA
      AF2df$RMSD6[1] <- NA
      AF2df$Nter_trunc[1] <- NA
      AF2df$Cter_trunc[1] <- NA
    }
  }, error = function(e) {error01 <<- TRUE})
  
  if (error01){
    AF2df$Dom[1] <-  CATHUNP2$Dom[d]
    AF2df$SP_PRIMARY[1] <- CATHUNP2$SP_PRIMARY[d]
    AF2df$af2ResBoundary[1] <- NA
    AF2df$af2DomLen[1] <- NA
    AF2df$PID[1] <- NA
    AF2df$PID1[1] <- NA
    AF2df$pLDDTq10[1] <- NA
    AF2df$pLDDTq20[1] <- NA
    AF2df$pLDDTq30[1] <- NA
    AF2df$RMSD6[1] <- NA
    AF2df$Nter_trunc[1] <- NA
    AF2df$Cter_trunc[1] <- NA
  }

  outdf <- as.data.frame(AF2df)
  if(d==1){
    write.table(outdf, file = file.path(outputdir, outputfile), row.names = F, col.names = T, quote = F, sep = "\t", append = F)
  }else{
    write.table(outdf, file = file.path(outputdir, outputfile), row.names = F, col.names = F, quote = F, sep = "\t", append = T)
  }
}



