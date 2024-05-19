
readSTRIDE <- function(SF) {
  
  out <- tryCatch(
    {
      
      #message("This is the 'try' part")
      
      read.table(SF, header = TRUE) %>% is.data.frame()
      
    },
    error=function(cond) {
      
      
      return(NULL)
    },
    warning=function(cond) {
      
      return(NULL)
    }
  )
  return(out)
}



getSecStr <- function(i){
  dom <- CATHUNP2$Dom[i]
  
  # cath and af2 stride file check
  cthStrideFC <- readSTRIDE(file.path(CATH_DB_RSA_PATH, paste0(dom, ".rsa")))
  af2StrideFC <- readSTRIDE(file.path(AF2_DOM_DB_RSA_PATH, paste0(dom, "_af2.rsa")))
  
  sec_str <- list()
  sec_str$Dom[1] <- dom
  
  if((!is.null(cthStrideFC)) & (!is.null(af2StrideFC))){
    
    
    cdom_secstr <- read.table(file.path(CATH_DB_RSA_PATH, paste0(dom, ".rsa")), header = T)
    adom_secstr <- read.table(file.path(AF2_DOM_DB_RSA_PATH, paste0(dom, "_af2.rsa")), header = T)
    

    
    if(length(cdom_secstr[,1]) > 20){
      
      cth_Helix_count <- sum(cdom_secstr$SSE_short == "H")
      cth_Sheet_count <- sum(cdom_secstr$SSE_short == "E")
      cth_secstr_count <- cth_Helix_count + cth_Sheet_count
      cth_loop_count <- cdom_secstr %>% filter(SSE_short !="E") %>% filter(SSE_short !="H") %>% .$SSE_short %>% length()
      
      af2_Helix_count <- sum(adom_secstr$SSE_short == "H")
      af2_Sheet_count <- sum(adom_secstr$SSE_short == "E")
      af2_secstr_count <- af2_Helix_count + af2_Sheet_count
      af2_loop_count <- adom_secstr %>% filter(SSE_short !="E") %>% filter(SSE_short !="H") %>% .$SSE_short %>% length()
      

      sec_str$cth_Helix_count[1] <- cth_Helix_count
      sec_str$cth_Sheet_count[1] <- cth_Sheet_count
      sec_str$cth_secstr_count[1] <- cth_secstr_count
      sec_str$cth_loop_count[1] <- cth_loop_count
      
      sec_str$af2_Helix_count[1] <- af2_Helix_count
      sec_str$af2_Sheet_count[1] <- af2_Sheet_count
      sec_str$af2_secstr_count[1] <- af2_secstr_count
      sec_str$af2_loop_count[1] <- af2_loop_count
    } else{
      sec_str$cth_Helix_count[1] <- NA
      sec_str$cth_Sheet_count[1] <- NA
      sec_str$cth_secstr_count[1] <- NA
      sec_str$cth_loop_count[1] <- NA
      
      sec_str$af2_Helix_count[1] <- NA
      sec_str$af2_Sheet_count[1] <- NA
      sec_str$af2_secstr_count[1] <- NA
      sec_str$af2_loop_count[1] <- NA
    }
    
    
  } else{
    sec_str$cth_Helix_count[1] <- NA
    sec_str$cth_Sheet_count[1] <- NA
    sec_str$cth_secstr_count[1] <- NA
    sec_str$cth_loop_count[1] <- NA
    
    sec_str$af2_Helix_count[1] <- NA
    sec_str$af2_Sheet_count[1] <- NA
    sec_str$af2_secstr_count[1] <- NA
    sec_str$af2_loop_count[1] <- NA

  }
  sec_str <- as.data.frame(sec_str)
  return(sec_str)
}
