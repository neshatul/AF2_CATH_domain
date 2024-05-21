
library(dplyr)


v3vol <- function(i){
  #i=85410
  pref <- CATHUNP2$Dom[i]
  volsurf=list()
  if(!dir.exists("temp")){dir.create("temp")}
  # temp file name random
  TEMPfn <- tempfile(pattern = "pdb2xyzr_", fileext=".pdb", tmpdir = "temp")
  #Remove Hetero Atom (Water, Salt, etc)
  #egrep "^ATOM  " 1A01.pdb > 1a01-noions.pdb
  system(paste0("egrep ", "'^", "ATOM  ",  "' ",file.path(CATH_DB_PATH, paste0(pref, ".pdb"))," > ", TEMPfn))
  
  TEMPfn2 <- tempfile(pattern = "pdb2xyzr_", fileext=".xyzr", tmpdir = "temp")
  #pdb_to_xyzr ../1a01-noions.pdb > ../1a01-noions.xyzr
  system(paste0("sh ", PDB_to_XYZR, " ", TEMPfn ," > ", TEMPfn2), ignore.stderr = T, intern = F)
  
  if(file.exists(TEMPfn)){
    invisible(file.remove(TEMPfn))
  }
  
  #Volume.exe -i 1a01-noions.xyzr -p 1.5 -g 0.5
  
  #vol and surf are output
  # PATH to Volume.exe is added to the $PATH so can be called directlys 
  vs <- system(paste0('Volume.exe -i ', TEMPfn2, " -p 1.5 -g 0.5"), ignore.stderr = T, intern = T)
  
  if(file.exists(TEMPfn2)){
    invisible(file.remove(TEMPfn2))
  }
  
  #try(log("not a number"), silent = TRUE)
  try(vs_split <- strsplit(vs, split = "\t")[[1]], silent = TRUE)
  
  if(length(vs_split) == 8){
    
    vol <- vs_split[3] %>% as.numeric()
    area <- vs_split[5] %>% as.numeric()
    volsurf$Dom[1] <- CATHUNP2$Dom[i]
    volsurf$SP_PRIMARY[1] <- CATHUNP2$SP_PRIMARY[i]
    volsurf$Vol[1] <- vol
    volsurf$Area[1] <- area
    volsurf <- as.data.frame(volsurf)
    return(volsurf)
  }else{

    volsurf$Dom[1] <- CATHUNP2$Dom[i]
    volsurf$SP_PRIMARY[1] <- CATHUNP2$SP_PRIMARY[i]
    volsurf$Vol[1] <- NA
    volsurf$Area[1] <- NA
    volsurf <- as.data.frame(volsurf)
    return(volsurf)
  }

}





