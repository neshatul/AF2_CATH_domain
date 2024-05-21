## Data used is our latest updated one developed using the filter like pLDDT>80, 
## discarding the Nter and Cter extension tagged and PID1>99



df02_110k <- TRUE
merge.keep.all <- TRUE

co_pos <- 10
co_neg <- -10  #25
CO <- 1000 # Limit the abs(diff) to this value


if(df02_110k) {
  dataframe <- df02
  GO_annotation.tsv <- "GO_df02_annotation_df02_110k.tsv"
  if(merge.keep.all){
    outfile <- "GO_df02_110k_MergeKeepAll02_long"
  } else {outfile <- "GO_df02_110k_NoMergeKeepAll02"}
}else{
  dataframe <- df01
  GO_annotation.tsv <- "GO_df01_annotation_130k.tsv"
  if(merge.keep.all){
    outfile <- "GO_df01_130k_MergeKeepAll02"
  } else {outfile <- "GO_df01_130k_NoMergeKeepAll02"}
}


#'@description
#'FATCAT median RMSD median for fold group/cluster
#'
FCrq50 <- dataframe %>% group_by(fold) %>% summarise(Rq50 = median(FCrmsd))

#'@description
#'Canonical median RMSD median for fold group/cluster
#'
Rq50 <- dataframe %>% group_by(fold) %>% summarise(Rq50 = median(RMSD))


#'@description Download GO annotations from Uniprot
#' Unique uniprot ID obtained from df01. One UniID can be mapped to multiple CATH 
#' domain obtained from one or multiple experimets.
#' Use Unique UniID,  to obtain GO term annotations.
if(!file.exists(file.path(data, GO_annotation.tsv))){
  UId_string <- df01$SP_PRIMARY %>% unique() %>% paste0(., collapse = ",") #string for ID mapping search
  
  files = list(
    ids = UId_string,
    from = "UniProtKB_AC-ID",
    to = "UniProtKB"
  )
  r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, encode = "multipart", accept_json())
  submission <- content(r, as = "parsed")
  
  if (isJobReady(submission[["jobId"]])) {
    url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
    r <- GET(url = url, accept_json())
    details <- content(r, as = "parsed")
    url <- getResultsURL(details[["redirectURL"]])
    # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
    url <- paste(url, "?fields=accession%2Creviewed%2Cprotein_name%2Corganism_name%2Clength%2Cgene_names%2Cgo%2Cgo_p%2Cgo_c%2Cgo_f%2Cgo_id&format=tsv", sep = "")
    r <- GET(url = url)
    bin <- content(r, "raw")
    writeBin(bin, file.path(data, GO_annotation.tsv))
  }
}



# unp_fun_anno_FN <- "GO_semantic_data_2024_04_26.tsv.gz"
# 
# GOdf <- read.csv(file.path(data, unp_fun_anno_FN), sep="\t")%>% # GO data
#   filter(grepl("Human", Organism)) # GO data for human

GOdf <- read.csv(file.path(data, GO_annotation.tsv), sep="\t") %>% # GO data
  filter(grepl("Human", Organism)) # GO data for human

head(GOdf); dim(GOdf)


#'@description
#' list of Unique GO terms and protein name
#' 
GO_bp_all_unq <- GOdf %>% select(Gene.Ontology..biological.process.) %>% getElement(1) %>% 
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>% unique() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>% unique()

GO_cc_all_unq <- GOdf %>% select(Gene.Ontology..cellular.component.) %>% getElement(1) %>% 
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>% unique() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>% unique()

GO_mf_all_unq <- GOdf %>% select(Gene.Ontology..molecular.function.) %>% getElement(1) %>% 
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>% unique() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>% unique()

#'@description
#'List of GO terms in Low 10% of domains/folds. Rq50 around 0.00
#'Proteins are not repeated.
#' Since we talking about Low 10% it remains same for median of RMSD and FCrmsd
#' Finally generate frequency table for different GO terms

Low_cutoff <- dim(FCrq50)[1] * 0.1
FCrq50_folds <- FCrq50 %>% arrange(., Rq50) %>% .[1:Low_cutoff,] %>% select(fold) %>% getElement(1)
RMSDq50_folds <- Rq50 %>% arrange(., Rq50) %>% .[1:Low_cutoff,] %>% select(fold) %>% getElement(1)


FCrq50_Uni <- dataframe %>% filter(fold %in% FCrq50_folds) %>% select(SP_PRIMARY) %>% unique() %>% getElement(1)
RMSDq50_Uni <- dataframe %>% filter(fold %in% RMSDq50_folds) %>% select(SP_PRIMARY) %>% unique() %>% getElement(1)


# GO terms frequency table for FCrmsd
FCLow10_GObp <- GOdf %>% filter(Entry %in% FCrq50_Uni) %>% 
  select(Gene.Ontology..biological.process.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

FCLow10_GOcc <- GOdf %>% filter(Entry %in% FCrq50_Uni) %>% 
  select(Gene.Ontology..cellular.component.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

FCLow10_GOmf <- GOdf %>% filter(Entry %in% FCrq50_Uni) %>% 
  select(Gene.Ontology..molecular.function.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()


# GO terms frequency table for RMSD
Low10_GObp <- GOdf %>% filter(Entry %in% RMSDq50_Uni) %>% 
  select(Gene.Ontology..biological.process.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

Low10_GOcc <- GOdf %>% filter(Entry %in% RMSDq50_Uni) %>% 
  select(Gene.Ontology..cellular.component.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

Low10_GOmf <- GOdf %>% filter(Entry %in% RMSDq50_Uni) %>% 
  select(Gene.Ontology..molecular.function.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()


#'@description
#'Frequency distribution table for High End RMSD and FCrmsd
#'


FCrq50_foldsHi <- FCrq50 %>% filter(Rq50 >= 2.5) %>% select(fold) %>% getElement(1)
RMSDq50_foldsHi <- Rq50 %>% filter(Rq50 >= 2.5) %>% select(fold) %>% getElement(1)

FCrq50_UniHi <- df01 %>% filter(fold %in% FCrq50_foldsHi) %>% select(SP_PRIMARY) %>% unique() %>% getElement(1)
RMSDq50_UniHi <- df01 %>% filter(fold %in% RMSDq50_foldsHi) %>% select(SP_PRIMARY) %>% unique() %>% getElement(1)

# Frequency table of GO terms for High end FCrmsd
FCHi_GObp <- GOdf %>% filter(Entry %in% FCrq50_UniHi) %>% 
  select(Gene.Ontology..biological.process.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

FCHi_GOcc <- GOdf %>% filter(Entry %in% FCrq50_UniHi) %>% 
  select(Gene.Ontology..cellular.component.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

FCHi_GOmf <- GOdf %>% filter(Entry %in% FCrq50_UniHi) %>% 
  select(Gene.Ontology..molecular.function.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

# Frequency table of GO terms for High end RMSD
RmsdHi_GObp <- GOdf %>% filter(Entry %in% RMSDq50_UniHi) %>% 
  select(Gene.Ontology..biological.process.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

RmsdHi_GOcc <- GOdf %>% filter(Entry %in% RMSDq50_UniHi) %>% 
  select(Gene.Ontology..cellular.component.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()

RmsdHi_GOmf <- GOdf %>% filter(Entry %in% RMSDq50_UniHi) %>% 
  select(Gene.Ontology..molecular.function.) %>% getElement(1) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){gsub("; ", ";", x)}, mc.cores = NUM_CORES ) %>% unlist() %>%
  pbmcapply::pbmclapply(., function(x){strsplit(., split=";")}, mc.cores = NUM_CORES) %>% unlist() %>%
  table() %>% as.data.frame()


#'@description
#'Differential distribution of GO terms 
#'Difference is calculated as difference between freq of GO terms in Low10 and Higher end of RMSD
#'
#'Diff  = FREQlow10   -  FREQhigh
#'


# merge.keep.all <- FALSE declared on the top
# if(merge.keep.all){
#   outfile <- "GO_MergeKeepAll00.pdf"
# } else {outfile <- "GO_MergeKeepAll00_NO.pdf"}
# co_pos <- 6
# co_neg <- -6


# FREQ difference in FCrmsd          BP
names(FCLow10_GObp) <- c("GO", "Freq")
names(FCHi_GObp) <- c("GO", "Freq")
DiffBP_FC <- merge(FCLow10_GObp, FCHi_GObp, by="GO", all=merge.keep.all) %>% 
  mutate(GO = gsub(" \\[GO:.*\\]", "", GO)) %>%
  mutate(Freq.x = if_else(is.na(Freq.x), 0, Freq.x)) %>%
  mutate(Freq.y = if_else(is.na(Freq.y), 0, Freq.y)) %>%
  mutate(Diff= (Freq.y - Freq.x)) %>% arrange(., -Diff) %>%  
  filter((Diff > co_pos | Diff < co_neg)) %>% # consider if the difference is above  and below 6
  mutate(Col=ifelse(Diff > co_pos, "High", ifelse(Diff < co_neg,  "Low", NA)))

# FREQ difference in FCrmsd          CC
names(FCLow10_GOcc) <- c("GO", "Freq")
names(FCHi_GOcc) <- c("GO", "Freq")
DiffCC_FC <- merge(FCLow10_GOcc, FCHi_GOcc, by="GO", all=merge.keep.all) %>% 
  mutate(GO = gsub(" \\[GO:.*\\]", "", GO)) %>%
  mutate(Freq.x = if_else(is.na(Freq.x), 0, Freq.x)) %>%
  mutate(Freq.y = if_else(is.na(Freq.y), 0, Freq.y)) %>%
  mutate(Diff= (Freq.y - Freq.x)) %>% arrange(., -Diff) %>%  
  filter((Diff > co_pos | Diff < co_neg)) %>% # consider if the difference is above  and below 6
  mutate(Col=ifelse(Diff > co_pos, "High", ifelse(Diff < co_neg,  "Low", NA)))

# FREQ difference in FCrmsd          MF
names(FCLow10_GOmf) <- c("GO", "Freq")
names(FCHi_GOmf) <- c("GO", "Freq")
DiffMF_FC <- merge(FCLow10_GOmf, FCHi_GOmf, by="GO", all=merge.keep.all) %>% 
  mutate(GO = gsub(" \\[GO:.*\\]", "", GO)) %>%
  mutate(Freq.x = if_else(is.na(Freq.x), 0, Freq.x)) %>%
  mutate(Freq.y = if_else(is.na(Freq.y), 0, Freq.y)) %>%
  mutate(Diff= (Freq.y - Freq.x)) %>% arrange(., -Diff) %>%  
  filter((Diff > co_pos | Diff < co_neg)) %>% # consider if the difference is above  and below 6
  mutate(Col=ifelse(Diff > co_pos, "High", ifelse(Diff < co_neg,  "Low", NA)))


# FREQ difference in RMSD          BP
names(Low10_GObp) <- c("GO", "Freq")
names(RmsdHi_GObp) <- c("GO", "Freq")
DiffBP <- merge(Low10_GObp, RmsdHi_GObp, by="GO", all=merge.keep.all) %>%
  mutate(GO = gsub(" \\[GO:.*\\]", "", GO)) %>%
  mutate(Freq.x = if_else(is.na(Freq.x), 0, Freq.x)) %>%
  mutate(Freq.y = if_else(is.na(Freq.y), 0, Freq.y)) %>%
  mutate(Diff= (Freq.y - Freq.x)) %>% arrange(., -Diff) %>%  
  filter((Diff > co_pos | Diff < co_neg)) %>% # consider if the difference is above  and below 6
  mutate(Col=ifelse(Diff > co_pos, "High", ifelse(Diff < co_neg,  "Low", NA)))

# FREQ difference in RMSD          CC
names(Low10_GOcc) <- c("GO", "Freq")
names(RmsdHi_GOcc) <- c("GO", "Freq")
DiffCC <- merge(Low10_GOcc, RmsdHi_GOcc, by="GO", all=merge.keep.all) %>%
  mutate(GO = gsub(" \\[GO:.*\\]", "", GO)) %>%
  mutate(Freq.x = if_else(is.na(Freq.x), 0, Freq.x)) %>%
  mutate(Freq.y = if_else(is.na(Freq.y), 0, Freq.y)) %>%
  mutate(Diff= (Freq.y - Freq.x)) %>% arrange(., -Diff) %>%  
  filter((Diff > co_pos | Diff < co_neg)) %>% # consider if the difference is above  and below 6
  mutate(Col=ifelse(Diff > co_pos, "High", ifelse(Diff < co_neg,  "Low", NA)))

# FREQ difference in RMSD          MF
names(Low10_GOmf) <- c("GO", "Freq")
names(RmsdHi_GOmf) <- c("GO", "Freq")
DiffMF <- merge(Low10_GOmf, RmsdHi_GOmf, by="GO", all=merge.keep.all) %>%
  mutate(GO = gsub(" \\[GO:.*\\]", "", GO)) %>%
  mutate(Freq.x = if_else(is.na(Freq.x), 0, Freq.x)) %>%
  mutate(Freq.y = if_else(is.na(Freq.y), 0, Freq.y)) %>%
  mutate(Diff= (Freq.y - Freq.x)) %>% arrange(., -Diff) %>%  
  filter((Diff > co_pos | Diff < co_neg)) %>% # consider if the difference is above  and below 6
  mutate(Col=ifelse(Diff > co_pos, "High", ifelse(Diff < co_neg,  "Low", NA)))






plot.go <- function(data=NULL, plot.title=NULL, plot.sub.title=NULL){
  theme_base_size <- 10
  ytextsize <- 6
  y_lab_wrap_length <- 30
  ggplot(data, aes(x=Diff, y=factor(GO, levels = GO))) +
    geom_point(alpha=1, size = 4, aes(col=Col)) +
    #scale_colour_gradientn(colours = redblue(20), values=seq(-10, 20, length.out=30)/30) +
    scale_color_manual(values=c("High"="red", "Low"="blue"))+#, limits=c(-5, 14))+
    theme_bw(base_size = theme_base_size)+
    labs(
      title = plot.title,
      subtitle = paste0("RMSD (High - Low) | Xlim ( -", plot.sub.title , " to ", plot.sub.title," )" ),
      caption = stringr::str_interp("Count diff |${co_pos} & ${co_neg}| are removed."),
      #tag = "Figure 1",
      x = "Difference of count",
      y = ""
    )+
    scale_y_discrete(labels = function(x) str_wrap(str_to_title(x), 
                                                   width = y_lab_wrap_length))
}


plot.go2 <- function(data=NULL, plot.title=NULL, plot.sub.title=NULL){
  theme_base_size <- 10
  ytextsize <- 6
  y_lab_wrap_length <- 30
  ggplot(data, aes(x=Diff, y=factor(GO, levels = GO))) +
    geom_point(alpha=1, size = 4, aes(col=Col)) +
    #scale_colour_gradientn(colours = redblue(20), values=seq(-10, 20, length.out=30)/30) +
    scale_color_manual(values=c("High"="red", "Low"="blue"))+#, limits=c(-5, 14))+
    theme_bw(base_size = theme_base_size)+
    labs(
      title = plot.title,
#      subtitle = paste0("RMSD (High - Low) | Xlim ( -", plot.sub.title , " to ", plot.sub.title," )" ),
#      caption = stringr::str_interp("Count diff |${co_pos} & ${co_neg}| are removed."),
      #tag = "Figure 1",
      x = "Difference of count",
      y = ""
    )
}

# when unwrap is not used longest name on y-axis defines the space occupied on y-acis name.
# longest word is "positive regulation of vascular associated smooth muscle cell differentiation involved in phenotypic switching"
# So we adjust one name in every data frame to have longest word
Add_space_to_first_last_GO_TERM <- function(dfr=NULL){
  row_last <- dim(dfr)[1]
  lw <- strsplit("positive regulation of vascular associated smooth
muscle cell differentiation involved in phenotypic switching", split="")[[1]] %>% 
    length()
  space_to_add1 <- lw - length(strsplit(dfr$GO[1], split="")[[1]])
  space_to_add2 <- lw - length(strsplit(dfr$GO[row_last], split="")[[1]])
  dfr$GO[1] <- paste0(paste0(rep(" ", space_to_add1+5), collapse = ""), dfr$GO[1], collapse = "")
  dfr$GO[row_last] <- paste0(paste0(rep(" ", space_to_add2+5), collapse = ""), dfr$GO[row_last], collapse = "")
  return(dfr)
}




DiffBP <- DiffBP %>% Add_space_to_first_last_GO_TERM() %>% filter((Diff < CO & Diff > -CO))
p.r.bp <- plot.go(DiffBP, "RMSD GO_BP", CO)

DiffCC <- DiffCC %>% Add_space_to_first_last_GO_TERM() %>% filter((Diff < CO & Diff > -CO))
p.r.cc <- plot.go(DiffCC, "RMSD GO_CC", CO)

DiffMF <- DiffMF %>% Add_space_to_first_last_GO_TERM() %>% filter((Diff < CO & Diff > -CO))
p.r.mf <- plot.go(DiffMF, "RMSD GO_MF", CO)



DiffBP_FC <- DiffBP_FC %>% Add_space_to_first_last_GO_TERM() %>% filter((Diff < CO & Diff > -CO))
p.f.bp <- plot.go(DiffBP_FC, "FCrmsd GO_BP", CO)

DiffCC_FC <- DiffCC_FC %>% Add_space_to_first_last_GO_TERM() %>% filter((Diff < CO & Diff > -CO))
p.f.cc <- plot.go(DiffCC_FC, "FCrmsd GO_CC", CO)

DiffMF_FC <- DiffMF_FC %>% Add_space_to_first_last_GO_TERM() %>% filter((Diff < CO & Diff > -CO))
p.f.mf <- plot.go(DiffMF_FC, "FCrmsd GO_MF", CO)



saveplot <- TRUE
if(saveplot){
  pdf(file.path(data, paste0(outfile, ".pdf")), height = 7.50, width = 14.25)
}

require(gridExtra)
grid.arrange(p.r.bp, p.r.cc, p.r.mf, p.f.bp, p.f.cc, p.f.mf, ncol=3)

if(saveplot){dev.off()}

saveplot <- TRUE
if(saveplot){
  pdf(file.path(data, paste0(outfile, "_02.pdf")), height = 12, width = 8.5)
}
DiffBP <- DiffBP %>% filter((Diff < CO & Diff > -CO))
p.r.bp <- DiffBP %>% plot.go2(., "RMSD GO_BP", CO)
p.r.bp_pos <- DiffBP %>% filter(Diff>0) %>% plot.go2(., "RMSD GO_BP positive", CO)
p.r.bp_neg <- DiffBP %>% filter(Diff<0) %>% plot.go2(., "RMSD GO_BP negative", CO)

DiffCC <- DiffCC %>% filter((Diff < CO & Diff > -CO))
p.r.cc <- DiffCC %>% plot.go2(., "RMSD GO_CC", CO)
p.r.cc_pos <- DiffCC %>% filter(Diff>0) %>% plot.go2(., "RMSD GO_CC positive", CO)
p.r.cc_neg <- DiffCC %>% filter(Diff<0) %>% plot.go2(., "RMSD GO_CC negative", CO)

DiffMF <- DiffMF %>% filter((Diff < CO & Diff > -CO))
p.r.mf <- DiffMF %>% plot.go2(., "RMSD GO_MF", CO)
p.r.mf_pos <- DiffMF %>% filter(Diff>0) %>% plot.go2(., "RMSD GO_MF", CO)
p.r.mf_neg <- DiffMF %>% filter(Diff<0) %>% plot.go2(., "RMSD GO_MF", CO)



DiffBP_FC <- DiffBP_FC %>% filter((Diff < CO & Diff > -CO))
p.f.bp <- plot.go2(DiffBP_FC, "FCrmsd GO_BP", CO)

DiffCC_FC <- DiffCC_FC %>% filter((Diff < CO & Diff > -CO))
p.f.cc <- plot.go2(DiffCC_FC, "FCrmsd GO_CC", CO)

DiffMF_FC <- DiffMF_FC %>% filter((Diff < CO & Diff > -CO))
p.f.mf <- plot.go2(DiffMF_FC, "FCrmsd GO_MF", CO)


print(p.r.bp)
print(p.r.cc)
print(p.r.mf)

# print(p.r.bp_pos)
# print(p.r.bp_neg)
# print(p.r.cc_pos)
# print(p.r.cc_neg)
# print(p.r.mf_pos)
# print(p.r.mf_neg)

print(p.f.bp)
print(p.f.cc)
print(p.f.mf)


if(saveplot){dev.off()}








