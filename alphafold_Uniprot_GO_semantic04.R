library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
# colorRampPalette "crp" bw is blue white and rw is red white
crp_bw <- colorRampPalette(c("blue","white"), space="rgb")
crp_rw <- colorRampPalette(c("red","white"), space="rgb")

# function section
ExpandGotermAsList <- function(unp_GoAnno = FALSE){
  #unp_GoAnno are usually like a list as follows
  #[1] "acetylcholine receptor activator activity [GO:0030549]; acetylcholine receptor binding [GO:0033130]"
  #[2] "amylin binding [GO:0097645]; apolipoprotein binding [GO:0034185]"
  #[3] "chemoattractant activity [GO:0042056]; chromatin binding [GO:0003682]; DNA binding [GO:0003677]"
  # The function clubs all these GO texts and GOterms as a dataframe.
  AnnoMerged <- c()
  for(i in unp_GoAnno){
    tmp <- strsplit(i, ';')[[1]]
    AnnoMerged <- c(AnnoMerged, stringr::str_trim(tmp))
  }
  
  AnnoList <- list()
  c=0
  for(i in AnnoMerged){
    c <- c+1
    tmp2_1 <- gsub('(.*)(\\[GO\\:[0-9].*\\])', '\\1\\',i)
    tmp2_2 <- stringr::str_trim(tmp2_1)
    AnnoList$Text[c] <- tmp2_2
    
    tmp3_1 <- gsub('(.*)(\\[GO\\:[0-9].*\\])', '\\2\\',i)
    tmp3_2 <- stringr::str_trim(tmp3_1)
    AnnoList$GOterm[c] <- tmp3_2
  }
  return(data.frame(AnnoList))
}

# Expand the go term in a dataframe
ExpandGotermDF <- function(unp_GoAnnoDF = FALSE, GO=FALSE){
  # GO are  Gene.Ontology..biological.process.  or  Gene.Ontology..cellular.component.
  # or Gene.Ontology..molecular.function.
  # unp_GoAnnoDF <- Last10_df
  # GO <- "Gene.Ontology..biological.process."
  # unp_GoAnnoDF is like a dataframe as follows
  # c("Entry", "Rq50", "Gene.Ontology..biological.process.", "Gene.Ontology..cellular.component.", 
  #   "Gene.Ontology..molecular.function.")
  
  temp_df <- as.data.frame(matrix(NA, ncol = 4, nrow = 0))
  
  for(i in seq_along(unp_GoAnnoDF[[GO]])){
    GOText <- c()
    gotext <- unp_GoAnnoDF[[GO]][i]
    tmp <- strsplit(gotext, ';')[[1]]
    GOText <- c(GOText, stringr::str_trim(tmp))
    n <- length(GOText)
    Entry <- rep(unp_GoAnnoDF$Entry[i], n)
    Rq50min <- rep(unp_GoAnnoDF$Rq50min[i], n)
    Rq50max <- rep(unp_GoAnnoDF$Rq50max[i], n)
    temp_df <- rbind(temp_df, data.frame(Entry, GOText, Rq50min, Rq50max))
    
  }
  temp_df <- temp_df %>% separate(GOText, c('GO_Text', 'GO_Term'), sep=" \\[")
  temp_df$GO_Term <- stringr::str_replace(temp_df$GO_Term, "]", "")
  return(temp_df)
}



#extract one half of a matrix object in R 
GetMatrixOneHalf <- function(mat){
  matrix_oneHalf <- c()
  for(i in 1:dim(mat)[1]){
    
    if(i == 1){
      matrix_oneHalf <- matrix_oneHalf
    }else{
      tmp <- unname(mat[i, 1:(i-1)])
      matrix_oneHalf <- c(matrix_oneHalf, tmp)
    }
    
  }
  return(matrix_oneHalf)
}
################################################################################
#proj dir
projdir <- "~/projects/af2pdb_comparision_05_analysis/"
#unp_fun_anno_FN <- "AF2_3540_proteins_UniProt_data.tsv"
unp_fun_anno_FN <- "idmapping_2023_08_11.tsv.gz"

unp_fun_anno <- read.csv(file.path(projdir, unp_fun_anno_FN), sep="\t")

dim(unp_fun_anno)
head(unp_fun_anno)

unp_fun_anno$Protein.families[1:10]


# only human protein
unp_fun_anno <- unp_fun_anno[grepl("Human", unp_fun_anno$Organism),]
dim(unp_fun_anno)

# list all the possible terms
GO_bp_all <- c()
GO_cc_all <- c()
GO_mf_all <- c()

for(i in seq_along(unp_fun_anno[,1])){
  if(unp_fun_anno$Gene.Ontology..biological.process.[i] != ""){
    tmpbp <- strsplit(unp_fun_anno$Gene.Ontology..biological.process.[i], split=";")[[1]]
    GO_bp_all <- c(GO_bp_all, stringi::stri_trim(tmpbp))
  }
  
  if(unp_fun_anno$Gene.Ontology..cellular.component.[i] != ""){
    tmpcc <- strsplit(unp_fun_anno$Gene.Ontology..cellular.component.[i], split=";")[[1]]
    GO_cc_all <- c(GO_cc_all, stringr::str_trim(tmpcc))
  }
  
  if(unp_fun_anno$Gene.Ontology..molecular.function.[i] != ""){
    tmpmf <- strsplit(unp_fun_anno$Gene.Ontology..molecular.function.[i], split=";")[[1]]
    GO_mf_all <- c(GO_mf_all, stringr::str_trim(tmpmf))
  }
  
}
GO_bp_all_unq <- unique(GO_bp_all); length(GO_bp_all); length(GO_bp_all_unq)
GO_cc_all_unq <- unique(GO_cc_all); length(GO_cc_all); length(GO_cc_all_unq)
GO_mf_all_unq <- unique(GO_mf_all); length(GO_mf_all); length(GO_mf_all_unq)

# associate biological function with 3540 proteins and ~170k domains in human
biologicaldf <- df[names(df) %in% c("Dom", "RMSD", "SP_PRIMARY", "fold")] ;dim(biologicaldf)

biologicaldf <- biologicaldf[biologicaldf$SP_PRIMARY %in%unp_fun_anno$Entry,] ;dim(biologicaldf)

for(i in seq_along(unp_fun_anno[,1])){
  biologicaldf$GObp[biologicaldf$SP_PRIMARY == unp_fun_anno$Entry[i]] <- unp_fun_anno$Gene.Ontology..biological.process.[i]
  biologicaldf$GOcc[biologicaldf$SP_PRIMARY == unp_fun_anno$Entry[i]] <- unp_fun_anno$Gene.Ontology..cellular.component.[i]
  biologicaldf$GOmf[biologicaldf$SP_PRIMARY == unp_fun_anno$Entry[i]] <- unp_fun_anno$Gene.Ontology..molecular.function.[i]
}

#dffold_top10pcent <- dffold[order(dffold$Rq50, decreasing = T),][1:75,] # ~10 % of 744 domain
dffold_top10pcent <- dffold[dffold$Rq50 >= 2.5,] # ~ 6.58% of 744
unp_index_hu <- lapply(unp_fun_anno$Entry, function(x){which(grepl(x,dffold_top10pcent$UnipId)) }) %>% unlist()
dffold_top10pcent_hu <- dffold_top10pcent[unique(unp_index_hu),] # at least one protein is human

#uni_top10pcent_hu <- lapply(unp_fun_anno$Entry, function(x){if(any(grepl(x,dffold_top10pcent$UnipId))){return(x)} }) %>% unlist()

uni_top10pcent_hu <- c()
uni_top10pcent_hu_RMSDq50min <- c()
uni_top10pcent_hu_RMSDq50max <- c()
for(i in unp_fun_anno$Entry){
  if(any(grepl(i,dffold_top10pcent$UnipId))){
    uni_top10pcent_hu <- c(uni_top10pcent_hu, i)
    tmp <- grepl(i,dffold_top10pcent$UnipId)
    #tmp <- grepl("P01008",dffold_top10pcent$UnipId)
    uni_top10pcent_hu_RMSDq50max <- c(uni_top10pcent_hu_RMSDq50max, max(dffold_top10pcent$Rq50[tmp]))
    uni_top10pcent_hu_RMSDq50min <- c(uni_top10pcent_hu_RMSDq50min, min(dffold_top10pcent$Rq50[tmp]))
  }
}
Top10UniRmsd <- data.frame(Entry=uni_top10pcent_hu, Rq50min=uni_top10pcent_hu_RMSDq50min, Rq50max=uni_top10pcent_hu_RMSDq50max)
dftmp <- unp_fun_anno[unp_fun_anno$Entry %in% Top10UniRmsd$Entry,][c(2,19, 20, 21)]

Top10_df <- merge(Top10UniRmsd, dftmp, by.x="Entry", sort = FALSE)



dffold_last10pcent <- dffold[order(dffold$Rq50, decreasing = F),][1:75,]
unp_index_hu <- lapply(unp_fun_anno$Entry, function(x){which(grepl(x,dffold_last10pcent$UnipId)) }) %>% unlist()
dffold_last10pcent_hu <- dffold_last10pcent[unique(unp_index_hu),] # at least one protein is human

uni_last10pcent_hu <- lapply(unp_fun_anno$Entry, function(x){if(any(grepl(x,dffold_last10pcent$UnipId))){return(x)} }) %>% unlist()

uni_last10pcent_hu <- c()
uni_last10pcent_hu_Rq50min <- c()
uni_last10pcent_hu_Rq50max <- c()
for(i in unp_fun_anno$Entry){
  if(any(grepl(i,dffold_last10pcent$UnipId))){
    uni_last10pcent_hu <- c(uni_last10pcent_hu, i)
    
    tmp <- grepl(i,dffold_last10pcent$UnipId)
    uni_last10pcent_hu_Rq50max <- c(uni_last10pcent_hu_Rq50max, max(dffold_last10pcent$Rq50[tmp]))
    uni_last10pcent_hu_Rq50min <- c(uni_last10pcent_hu_Rq50min, min(dffold_last10pcent$Rq50[tmp]))
  }
}
Last10UniRmsd <- data.frame(Entry=uni_last10pcent_hu, Rq50min=uni_last10pcent_hu_Rq50min, Rq50max=uni_last10pcent_hu_Rq50max)

dftmp <- unp_fun_anno[unp_fun_anno$Entry %in% Last10UniRmsd$Entry,][c(2,19, 20, 21)]

Last10_df <- merge(Last10UniRmsd, dftmp, by.x="Entry", sort = FALSE)
head(Last10_df)




SAVEIMAGE <- T
if(SAVEIMAGE) pdf(file.path(projdir, "OverRepresentation_05.pdf"), height = 9, width = 3.5)
if(SAVEIMAGE) {par(mar = c(5, 20, 4, 2) + 0.1)} else {par(mar = c(5, 17, 4, 2) + 0.1)}
theme_base_size <- 9
ytextsize <- 6
y_lab_wrap_length <- 30
################################################################################
#                                            BP
co_pos <- 4
co_neg <- -2
Top10_df_long <- ExpandGotermDF(Top10_df, "Gene.Ontology..biological.process.")
Last10_df_long <- ExpandGotermDF(Last10_df, "Gene.Ontology..biological.process.")
BP_toplast10 <- rbind(Top10_df_long, Last10_df_long)
head(BP_toplast10)


GOmerge <- c(BP_toplast10$GO_Text)
length(GOmerge)
GOmerge <- unique(GOmerge) 
length(GOmerge)

toplast10 <- list()
for(i in seq_along(GOmerge)){
  #i=120
  tmp_df <- BP_toplast10[BP_toplast10$GO_Text == GOmerge[i],]
  toplast10$GOterm[i] <- GOmerge[i]
  toplast10$Text[i] <- tmp_df$GO_Term[1]
  toplast10$Entry[i] <- tmp_df$Entry[1]
  topcount <- sum(Top10_df_long$GO_Text == GOmerge[i])
  toplast10$Top10_Count[i] <- topcount
  lastcount <- sum(Last10_df_long$GO_Text == GOmerge[i])
  toplast10$Last10_Count[i] <- lastcount
  
  if(topcount >lastcount){
    toplast10$Rq50_Av[i] <- median(tmp_df$Rq50max)
  }else toplast10$Rq50_Av[i] <- median(tmp_df$Rq50min)
  #toplast10$Rq50_Av[i] <- round(sum(tmp_df$Rq50)/(dim(tmp_df)[1]), 2)
}

toplast10 <- as.data.frame(toplast10)
toplast10$Diff_Count <- toplast10$Top10_Count - toplast10$Last10_Count

toplast10_BP <- toplast10[order(toplast10$Diff_Count, decreasing = T),]
toplast10_BP <- toplast10_BP[((toplast10_BP$Diff_Count > co_pos) | (toplast10_BP$Diff_Count < co_neg)),]

head(toplast10_BP)
tail(toplast10_BP)


toplast10_BP$Text <- factor(toplast10_BP$Text, levels = toplast10_BP$Text)
toplast10_BP$GOterm <- factor(toplast10_BP$GOterm, levels = toplast10_BP$GOterm)
#col_ryb <- colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(dim(toplast10_BP)[1])
col_list <- c(crp_rw(sum(toplast10_BP$Diff_Count>=0)), rev(crp_bw(sum(toplast10_BP$Diff_Count <0))))

size=toplast10_BP$Rq50_Av
ggplot(toplast10_BP, aes(x=Diff_Count, y=GOterm, colour=Diff_Count)) +
  geom_point(alpha=1, size = size) +
  scale_colour_gradientn(colours = redblue(20), values=seq(-10, 20, length.out=30)/30) +
  scale_colour_gradientn(colours = redblue(20), limits=c(-5, 14))+
  geom_point(aes(color =Diff_Count, size = size))+
  theme_bw(base_size = theme_base_size)+
  labs(
    title = "Biological Process",
    subtitle = "Rq50_Av is mean of Medians",
    caption = stringr::str_interp("Count diff |${co_pos} & ${co_neg}| are removed."),
    #tag = "Figure 1",
    x = "Difference of count",
    y = ""
  )+
  scale_y_discrete(labels = function(x) str_wrap(str_to_title(x), 
                                                 width = y_lab_wrap_length))

################################################################################
#                                            MF
co_pos <- 4
co_neg <- -2
Top10_df_long <- ExpandGotermDF(Top10_df, "Gene.Ontology..molecular.function.")
Last10_df_long <- ExpandGotermDF(Last10_df, "Gene.Ontology..molecular.function.")
BP_toplast10 <- rbind(Top10_df_long, Last10_df_long)
head(BP_toplast10)


GOmerge <- c(BP_toplast10$GO_Text)
length(GOmerge)
GOmerge <- unique(GOmerge) 
length(GOmerge)

toplast10 <- list()
for(i in seq_along(GOmerge)){
  #i=120
  tmp_df <- BP_toplast10[BP_toplast10$GO_Text == GOmerge[i],]
  toplast10$GOterm[i] <- GOmerge[i]
  toplast10$Text[i] <- tmp_df$GO_Term[1]
  toplast10$Entry[i] <- tmp_df$Entry[1]
  topcount <- sum(Top10_df_long$GO_Text == GOmerge[i])
  toplast10$Top10_Count[i] <- topcount
  lastcount <- sum(Last10_df_long$GO_Text == GOmerge[i])
  toplast10$Last10_Count[i] <- lastcount
  
  if(topcount >lastcount){
    toplast10$Rq50_Av[i] <- median(tmp_df$Rq50max)
  }else toplast10$Rq50_Av[i] <- median(tmp_df$Rq50min)
  #toplast10$Rq50_Av[i] <- round(sum(tmp_df$Rq50)/(dim(tmp_df)[1]), 2)
}

toplast10 <- as.data.frame(toplast10)
toplast10$Diff_Count <- toplast10$Top10_Count - toplast10$Last10_Count

toplast10_BP <- toplast10[order(toplast10$Diff_Count, decreasing = T),]
toplast10_BP <- toplast10_BP[((toplast10_BP$Diff_Count > co_pos) | (toplast10_BP$Diff_Count < co_neg)),]

head(toplast10_BP)
tail(toplast10_BP)


toplast10_BP$Text <- factor(toplast10_BP$Text, levels = toplast10_BP$Text)
toplast10_BP$GOterm <- factor(toplast10_BP$GOterm, levels = toplast10_BP$GOterm)
#col_ryb <- colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(dim(toplast10_BP)[1])
col_list <- c(crp_rw(sum(toplast10_BP$Diff_Count>=0)), rev(crp_bw(sum(toplast10_BP$Diff_Count <0))))

size=toplast10_BP$Rq50_Av
ggplot(toplast10_BP, aes(x=Diff_Count, y=GOterm, colour=Diff_Count)) +
  geom_point(alpha=1, size = size) +
  scale_colour_gradientn(colours = redblue(20), values=seq(-10, 20, length.out=30)/30) +
  scale_colour_gradientn(colours = redblue(20), limits=c(-5, 14))+
  geom_point(aes(color =Diff_Count, size = size))+
  theme_bw(base_size = theme_base_size)+
  labs(
    title = "Molecular Function",
    subtitle = "Rq50_Av is mean of Medians",
    caption = stringr::str_interp("Count diff |${co_pos} & ${co_neg}| are removed."),
    #tag = "Figure 1",
    x = "Difference of count",
    y = ""
  )+
  
  scale_y_discrete(labels = function(x) str_wrap(str_to_title(x), 
                                                 width = y_lab_wrap_length))
################################################################################
#                                            CC
co_pos <- 3
co_neg <- -2
Top10_df_long <- ExpandGotermDF(Top10_df, "Gene.Ontology..cellular.component.")
Last10_df_long <- ExpandGotermDF(Last10_df, "Gene.Ontology..cellular.component.")
BP_toplast10 <- rbind(Top10_df_long, Last10_df_long)
head(BP_toplast10)


GOmerge <- c(BP_toplast10$GO_Text)
length(GOmerge)
GOmerge <- unique(GOmerge) 
length(GOmerge)

toplast10 <- list()
for(i in seq_along(GOmerge)){
  #i=120
  tmp_df <- BP_toplast10[BP_toplast10$GO_Text == GOmerge[i],]
  toplast10$GOterm[i] <- GOmerge[i]
  toplast10$Text[i] <- tmp_df$GO_Term[1]
  toplast10$Entry[i] <- tmp_df$Entry[1]
  topcount <- sum(Top10_df_long$GO_Text == GOmerge[i])
  toplast10$Top10_Count[i] <- topcount
  lastcount <- sum(Last10_df_long$GO_Text == GOmerge[i])
  toplast10$Last10_Count[i] <- lastcount
  
  if(topcount >lastcount){
    toplast10$Rq50_Av[i] <- median(tmp_df$Rq50max)
  }else toplast10$Rq50_Av[i] <- median(tmp_df$Rq50min)
  #toplast10$Rq50_Av[i] <- round(sum(tmp_df$Rq50)/(dim(tmp_df)[1]), 2)
}

toplast10 <- as.data.frame(toplast10)
toplast10$Diff_Count <- toplast10$Top10_Count - toplast10$Last10_Count

toplast10_BP <- toplast10[order(toplast10$Diff_Count, decreasing = T),]
toplast10_BP <- toplast10_BP[((toplast10_BP$Diff_Count > co_pos) | (toplast10_BP$Diff_Count < co_neg)),]

head(toplast10_BP)
tail(toplast10_BP)


toplast10_BP$Text <- factor(toplast10_BP$Text, levels = toplast10_BP$Text)
toplast10_BP$GOterm <- factor(toplast10_BP$GOterm, levels = toplast10_BP$GOterm)
#col_ryb <- colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(dim(toplast10_BP)[1])
col_list <- c(crp_rw(sum(toplast10_BP$Diff_Count>=0)), rev(crp_bw(sum(toplast10_BP$Diff_Count <0))))

#ggplot(toplast10_BP, aes(x = Diff_Count, y = Text, size = Rq50_Av, color=col_ryb))+
size=toplast10_BP$Rq50_Av
ggplot(toplast10_BP, aes(x=Diff_Count, y=GOterm, colour=Diff_Count)) +
  geom_point(alpha=1, size = size) +
  scale_colour_gradientn(colours = redblue(20), values=seq(-10, 20, length.out=30)/30) +
  scale_colour_gradientn(colours = redblue(20), limits=c(-5, 14))+
  theme_bw(base_size = theme_base_size)+
  labs(
    title = "CC",
    subtitle = "Rq50_Av is mean of Medians",
    caption = stringr::str_interp("Count diff |${co_pos} & ${co_neg}| are removed."),
    #tag = "Figure 1",
    x = "Difference of count",
    y = ""
  )+
  scale_y_discrete(labels = function(x) str_wrap(str_to_title(x), 
                                                 width = y_lab_wrap_length))


if(SAVEIMAGE) dev.off()





library(ggplot2)

# Define a color palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
redblue<-colorRampPalette(rev(c("red","orange","blue")))
size=toplast10_BP$Rq50_Av
# Create a plot with the defined color scale
p1 <- ggplot(toplast10_BP, aes(x=Diff_Count, y=GOterm, colour=Diff_Count)) +
  geom_point(alpha=1, size = size) +
  scale_colour_gradientn(colours = redblue(30), values=seq(-10, 20, length.out=30)/30) +
  ggtitle("Z: 0 - 100")

# Create another plot with the same color scale
p2 <- ggplot(d2, aes(x=Diff_Count, y=GOterm, colour=Diff_Count)) +
  geom_point(alpha=.5, size = 6) +
  scale_colour_gradientn(colours = redblue(30), values=seq(-10, 20, length.out=30)/30) +
  ggtitle("Z: 50 - 150")

# Add the color scale to both plots
sc <- scale_colour_gradientn(colours = redblue(20), limits=c(-5, 14))
p1 + sc
p2 + sc


