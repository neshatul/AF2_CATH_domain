
# THE SCRIPT ANNOTATES THE CATH DOMAIN WITH RelASA scores. 
# It is run after all the RelASA file are generated.
# library(parallel)
# library(dplyr)
# source("~/repos/workflows/AF2_protein_dynamics/lib_RCSB.R")
#local databases


aa_list <- c("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", 
             "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", 
             "TYR", "VAL")


# Source: http://us.expasy.org/tools/pscale/Hphob.Eisenberg.html
# Amino acid scale: Normalized consensus hydrophobicity scale
# Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
# Reference: J. Mol. Biol. 179:125-142 (1984)
#
# Amino acid scale values:
#
# Ala:  0.620
# Arg: -2.530
# Asn: -0.780
# Asp: -0.900
# Cys:  0.290
# Gln: -0.850
# Glu: -0.740
# Gly:  0.480
# His: -0.400
# Ile:  1.380
# Leu:  1.060
# Lys: -1.500
# Met:  0.640
# Phe:  1.190
# Pro:  0.120
# Ser: -0.180
# Thr: -0.050
# Trp:  0.810
# Tyr:  0.260
# Val:  1.080
# definition of hydrophobic residues in classical terms
H1 <- c("ALA", "CYS", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL", "PRO", "GLY")
P1 <- c("SER",  "THR", "HIS", "LYS", "ARG", "ASP", "ASN", "GLU", "GLN")



# Definition of hydrophobic residues as defined in our lab finding of B75 residue distribution
# hydrophobic2 is the list hydrophobic residues
H2 <- c("ALA", "CYS", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "VAL")
# polar2 is a list of polar and charged residue
P2 <- c("PRO", "SER", "GLY", "HIS", "THR", "GLU", "LYS", "ASN", "ASP", "ARG",  "GLN")


GetRelASACols <- function(i){
  
  pdbid <- strsplit(CATHUNP2$Dom[i], split='')[[1]][1:4] %>% paste0(., collapse = '')
  chn <- strsplit(CATHUNP2$Dom[i], split='')[[1]][5]
  #resno_vect <- HyphenSeg2Vec(CATHUNP2$Boundary[i])
  
  #domain 3d
  d3d <-read.pdb(file.path(CATH_DB_PATH, paste0(CATHUNP2$Dom[i], ".pdb")))
  resno_vect <- d3d$atom[d3d$calpha,]$resno
  
  tdf <- data.frame(Dom=character(), SP_PRIMARY=character(), 
                    BE80_ratio=double(), 
                    
                    B80HP_ratio1=double(),
                    E80HP_ratio1=double(),
                    
                    B80HP_ratio2=double(),
                    E80HP_ratio2=double(), 
                    
                    nrHelix=double(),
                    l5LH=character(),
                    
                    nrStrand=double(),
                    l5LS=character(),
                    
                    nrCoil=double(),
                    l5LC=character()
                    )
  
  
  if(file.exists(file.path(RSA_PATH, paste0(pdbid, ".rsa")))) {
    cond1 <- length(read.table(file.path(RSA_PATH, paste0(pdbid, ".rsa")), header = T)[,1]) > 0
    if(cond1){
      rsa_df <- read.table(file.path(RSA_PATH, paste0(pdbid, ".rsa")), header = T) %>% 
        .[.[["Chain"]] == chn,] %>%
        .[.[["Resno"]] %in% resno_vect,]
      
      #sec str length cutoff
      sslc <- 10
      Helix <- rsa_df %>% filter(SSE_short=="H")
      if(dim(Helix)[1] != 0){
        nrHelix <- dim(Helix)[1]                  # no of residues in helix
        wrHelix <- Vector2HyphenSeg(as.numeric(Helix$Resno))  # which residues in helix
        
        ###list of 5 longest helix with no of residue
        l5LH <- lapply(strsplit(wrHelix,',')[[1]], function(X){
        ll <- strsplit(X, "-")[[1]][1] %>% as.numeric() #lower limit
        ul <- strsplit(X, "-")[[1]][2] %>% as.numeric() #upper limit
        return(ul-ll+1)
        }) %>% unlist() %>% .[order(., decreasing=T)] %>% .[1:sslc] %>% paste0(., collapse=",")

        
      } else {
        nrHelix <- 0
        l5LH <- 0
      }
      
      Strand <- rsa_df %>% filter(SSE_short=="E")
      if(dim(Strand)[1] != 0){
        nrStrand <- dim(Strand)[1]                # no of residues in helix
        wrStrand <- Vector2HyphenSeg(as.numeric(Strand$Resno)) # which residue in strand
        
        ###list of 5 longest strand with no of residue
        l5LS <- lapply(strsplit(wrStrand,',')[[1]], function(X){
          ll <- strsplit(X, "-")[[1]][1] %>% as.numeric() #lower limit
          ul <- strsplit(X, "-")[[1]][2] %>% as.numeric() #upper limit
          return(ul-ll+1)
        }) %>% unlist() %>% .[order(., decreasing=T)] %>% .[1:sslc] %>% paste0(., collapse=",")
        
      } else {
        nrStrand <- 0
        l5LS <- 0
      }
      
      Coil <- rsa_df %>% filter(SSE_short!="H")  %>% filter(SSE_short!="E")
      if(dim(Coil)[1] != 0){
        nrCoil <- dim(Coil)[1] # no of residues in helix
        wrCoil <- Vector2HyphenSeg(as.numeric(Coil$Resno))
        
        #List of 5 longest coil
        l5LC <- lapply(strsplit(wrCoil,',')[[1]], function(X){
          ll <- strsplit(X, "-")[[1]][1] %>% as.numeric() #lower limit
          ul <- strsplit(X, "-")[[1]][2] %>% as.numeric() #upper limit
          return(ul-ll+1)
        }) %>% unlist() %>% .[order(., decreasing=T)] %>% .[1:sslc] %>% paste0(., collapse=",")
        
      } else {
        nrCoil <- 0
        l5LC <- 0
      }
       
      
      
      if(length(rsa_df[,1]) >0){
        rsa_df$RelBSA <- 100 - rsa_df$RelASA
        
        tl <- length(rsa_df[,1]) # tot_length
        # buried vs exposed
        # BExx_ratio is the ratio of
        #BE10_ratio means above 10 till 100
        #BE20_ratio means above 20 till 100
        
        BE80_ratio <- sum(rsa_df$RelBSA >= 80 ) / (tl - sum(rsa_df$RelBSA >= 80 ))
        
        B80HP_ratio1 <- sum(rsa_df[rsa_df$RelBSA >= 80,][["Resid"]] %in% H1) /
          sum(rsa_df[rsa_df$RelBSA >= 80,][["Resid"]] %in% P1)
        E80HP_ratio1 <- sum(rsa_df[!rsa_df$RelBSA >= 80,][["Resid"]] %in% H1) /
          sum(rsa_df[!rsa_df$RelBSA >= 80,][["Resid"]] %in% P1)
        
        
        
        B80HP_ratio2 <- sum(rsa_df[rsa_df$RelBSA >= 80,][["Resid"]] %in% H2) /
          sum(rsa_df[!rsa_df$RelBSA >= 80,][["Resid"]] %in% P2)
        E80HP_ratio2 <- sum(!rsa_df[rsa_df$RelBSA >= 80,][["Resid"]] %in% H2) /
          sum(rsa_df[!rsa_df$RelBSA >= 80,][["Resid"]] %in% P2)
        
        
        tdf[1,1] <- CATHUNP2$Dom[i]
        tdf[1,2] <- CATHUNP2$SP_PRIMARY[i]
        
        tdf[1,3] <- BE80_ratio
        tdf[1,4] <- B80HP_ratio1
        tdf[1,5] <- E80HP_ratio1
        tdf[1,6] <- B80HP_ratio2
        tdf[1,7] <- E80HP_ratio2
        tdf[1,8] <- nrHelix
        tdf[1,9] <- l5LH
        tdf[1,10] <- nrStrand
        tdf[1,11] <- l5LS
        tdf[1,12] <- nrCoil
        tdf[1,13] <- l5LC
        
        # print(i)
        # cat(paste0('Running ', i, "\n"))
      } else {
        tdf[1,1] <- CATHUNP2$Dom[i]
        tdf[1,2] <- CATHUNP2$SP_PRIMARY[i]
        tdf[1,3:length(tdf)] <- rep(NA, length(tdf)-2)
        #print(i)
      }
      
    } else {
      tdf[1,1] <- CATHUNP2$Dom[i]
      tdf[1,2] <- CATHUNP2$SP_PRIMARY[i]
      tdf[1,3:length(tdf)] <- rep(NA, length(tdf)-2)
      #print(i)
    }
    
  } else {
    tdf[1,1] <- CATHUNP2$Dom[i]
    tdf[1,2] <- CATHUNP2$SP_PRIMARY[i]
    tdf[1,3:length(tdf)] <- rep(NA, length(tdf) -2 )
    #print(i)
  }
  
  return(tdf)
}





