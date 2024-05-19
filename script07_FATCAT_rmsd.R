

## data obtained from basic output of FATCAT
getFATCAT <- function(i){
  tryCatch(exp={
    FC_rmsd <- list()
    dom <- CATHUNP2$Dom[i]
    i1 <- CATH_DB_PATH
    p1 <- paste0(dom, ".pdb")
    
    
    i2 <- AF2_DOM_DB_PATH
    p2 <- paste0(dom, "_af2.pdb")
    
    fc_cmd <- stringr::str_interp('${FATCAT_EXE} -i1 ${i1} -i2 ${i2} -p1 ${p1} -p2 ${p2} -b')
    
    out <- system(fc_cmd, intern = T) %>% suppressWarnings()
    
    FC_rmsd$Dom[1] <- dom
    FC_rmsd$FCtwist[1] <- out[[1]] %>% strsplit(., "") %>% getElement(1) %>% .[54:55] %>% paste0(., collapse = "") %>% as.numeric()
    FC_rmsd$FCrmsd[1] <- out[[1]] %>% strsplit(., "") %>% getElement(1) %>% .[73:76] %>% paste0(., collapse = "") %>% as.numeric()
    FC_rmsd$FCrmsd_NoTrans[1] <- out[[1]] %>% strsplit(., "") %>% getElement(1) %>% .[81:85] %>% paste0(., collapse = "") %>% as.numeric()
    return(FC_rmsd)
  },
  error = function(e){
    FC_rmsd$Dom[1] <- dom
    FC_rmsd$FCtwist[1] <- NA
    FC_rmsd$FCrmsd[1] <- NA
    FC_rmsd$FCrmsd_NoTrans[1] <- NA
    return(FC_rmsd)
  }
  )
}

## data obtained from alignment output of FATCAT
getFATCAT2 <- function(i){
  tryCatch(exp={
    FC_rmsd <- list()
    dom <- CATHUNP2$Dom[i]
    i1 <- CATH_DB_PATH
    p1 <- paste0(dom, ".pdb")
    
    
    i2 <- AF2_DOM_DB_PATH
    p2 <- paste0(dom, "_af2.pdb")
    
    fc_cmd <- stringr::str_interp('${FATCAT_EXE} -i1 ${i1} -i2 ${i2} -p1 ${p1} -p2 ${p2} -q')
    
    out <- system(fc_cmd, intern = T) %>% suppressWarnings()
    
    FC_rmsd$Dom[1] <- dom
    FC_rmsd$FCtwist[1] <- out[2] %>% strsplit(., split=" ") %>% getElement(1) %>% getElement(2) %>% as.numeric()
    FC_rmsd$FCrmsd[1] <- out[2] %>% strsplit(., split=" ") %>% getElement(1) %>% getElement(10) %>% as.numeric()
    FC_rmsd$FCrmsd_NoTrans[1] <- out[2] %>% strsplit(., split=" ") %>% getElement(1) %>% getElement(12) %>% as.numeric()
    FC_rmsd$FCgap_percent[1] <- out[2] %>% strsplit(., split=" ") %>% getElement(1) %>% getElement(19) %>% gsub("\\(|\\)|\\%", "", .) %>% as.numeric()
    FC_rmsd$FCpvalue[1]  <- out[3] %>% strsplit(., split=" ") %>% getElement(1) %>% getElement(2) %>% as.numeric()
    FC_rmsd$FCidendity[1]  <- out[3] %>% strsplit(., split=" ") %>% getElement(1) %>% getElement(6) %>% gsub("%", "", .) %>% as.numeric()
    FC_rmsd$FCsimilarity[1]  <- out[3] %>% strsplit(., split=" ") %>% getElement(1) %>% getElement(8) %>% gsub("%", "", .) %>% as.numeric()
    return(FC_rmsd)
  },
  error = function(e){
    FC_rmsd$Dom[1] <- dom
    FC_rmsd$FCtwist[1] <- NA
    FC_rmsd$FCrmsd[1] <- NA
    FC_rmsd$FCrmsd_NoTrans[1] <- NA
    FC_rmsd$FCgap_percent[1] <- NA
    FC_rmsd$FCpvalue[1]  <- NA
    FC_rmsd$FCidendity[1]  <- NA
    FC_rmsd$FCsimilarity[1]  <- NA
    return(FC_rmsd)
  }
  )
}