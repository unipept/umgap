setwd("/Users/Aranka/Documents/Unief/Thesis/Thesis/Kmer2consensus/OptiesVGL")

mat <- read.csv("HiSeq_PR_SE.csv")
mat.mrl <- read.csv("HiSeq_PRe_SE.MRL.csv")
matM <- read.csv("MiSeq_PR_SE.csv")
matM.mrl <- read.csv("MiSeq_PRe_SE.MRL.csv")
data_frame <- as.data.frame(mat)
data_frame.mrl <- as.data.frame(mat.mrl)
data_frameM <- as.data.frame(matM)
data_frameM.mrl <- as.data.frame(matM.mrl)
data_frame$length_type <- paste(data_frame$seed_length,data_frame$type, sep=";")
data_frame.mrl$length_type <- paste(data_frame.mrl$seed_length,data_frame$type, sep=";")
data_frameM$length_type <- paste(data_frameM$seed_length,data_frameM$type, sep=";")
data_frameM.mrl$length_type <- paste(data_frameM.mrl$seed_length,data_frameM$type, sep=";")

matM08.mrl <- read.csv("MiSeq_PRe_SE.BothTree08.csv")
data_frameM08.mrl <- as.data.frame(matM08.mrl)
data_frameM08.mrl$length_type <- paste(data_frameM08.mrl$seed_length,data_frameM08.mrl$type, sep=";")
mat08.mrl <- read.csv("HiSeq_PRe_SE.BothTree08.csv")
data_frame08.mrl <- as.data.frame(mat08.mrl)
data_frame08.mrl$length_type <- paste(data_frame08.mrl$seed_length,data_frame08.mrl$type, sep=";")
matM02.mrl <- read.csv("MiSeq_PRe_SE.BothTree02.csv")
data_frameM02.mrl <- as.data.frame(matM02.mrl)
data_frameM02.mrl$length_type <- paste(data_frameM02.mrl$seed_length,data_frameM02.mrl$type, sep=";")
mat02.mrl <- read.csv("HiSeq_PRe_SE.BothTree02.csv")
data_frame02.mrl <- as.data.frame(mat02.mrl)
data_frame02.mrl$length_type <- paste(data_frame02.mrl$seed_length,data_frame02.mrl$type, sep=";")

matM05.mrl <- read.csv("MiSeq_PRe_SE.BothTree05.csv")
data_frameM05.mrl <- as.data.frame(matM05.mrl)
data_frameM05.mrl$length_type <- paste(data_frameM05.mrl$seed_length,data_frameM05.mrl$type, sep=";")
mat05.mrl <- read.csv("HiSeq_PRe_SE.BothTree05.csv")
data_frame05.mrl <- as.data.frame(mat05.mrl)
data_frame05.mrl$length_type <- paste(data_frame05.mrl$seed_length,data_frame05.mrl$type, sep=";")

matM04.mrl <- read.csv("MiSeq_PRe_SE.BothTree04.csv")
data_frameM04.mrl <- as.data.frame(matM04.mrl)
data_frameM04.mrl$length_type <- paste(data_frameM04.mrl$seed_length,data_frameM04.mrl$type, sep=";")
mat04.mrl <- read.csv("HiSeq_PRe_SE.BothTree04.csv")
data_frame04.mrl <- as.data.frame(mat04.mrl)
data_frame04.mrl$length_type <- paste(data_frame04.mrl$seed_length,data_frame04.mrl$type, sep=";")

matM06.mrl <- read.csv("MiSeq_PRe_SE.BothTree06.csv")
data_frameM06.mrl <- as.data.frame(matM06.mrl)
data_frameM06.mrl$length_type <- paste(data_frameM06.mrl$seed_length,data_frameM06.mrl$type, sep=";")
mat06.mrl <- read.csv("HiSeq_PRe_SE.BothTree06.csv")
data_frame06.mrl <- as.data.frame(mat06.mrl)
data_frame06.mrl$length_type <- paste(data_frame06.mrl$seed_length,data_frame06.mrl$type, sep=";")

matM07.mrl <- read.csv("MiSeq_PRe_SE.BothTree07.csv")
data_frameM07.mrl <- as.data.frame(matM07.mrl)
data_frameM07.mrl$length_type <- paste(data_frameM07.mrl$seed_length,data_frameM07.mrl$type, sep=";")
mat07.mrl <- read.csv("HiSeq_PRe_SE.BothTree07.csv")
data_frame07.mrl <- as.data.frame(mat07.mrl)
data_frame07.mrl$length_type <- paste(data_frame07.mrl$seed_length,data_frame07.mrl$type, sep=";")


gap_pen28 <- c(which(mat[,"gap_pen"]==0.2),which(mat[,"gap_pen"]==0.8))
mat_reduced <- mat[-gap_pen28,]
data_frame_reduced <- as.data.frame(mat_reduced)
data_frame_reduced$length_type <- paste(data_frame_reduced$seed_length,data_frame_reduced$type, sep=";")

gap_pen28 <- c(which(mat.mrl[,"gap_pen"]==0.2),which(mat.mrl[,"gap_pen"]==0.8))
mat.mrl_reduced <- mat.mrl[-gap_pen28,]
data_frame.mrl_reduced <- as.data.frame(mat.mrl_reduced)
data_frame.mrl_reduced$length_type <- paste(data_frame.mrl_reduced$seed_length,data_frame.mrl_reduced$type, sep=";")

gap_pen28M <- c(which(matM[,"gap_pen"]==0.2),which(matM[,"gap_pen"]==0.8))
matM_reduced <- matM[-gap_pen28M,]
data_frameM_reduced <- as.data.frame(matM_reduced)
data_frameM_reduced$length_type <- paste(data_frameM_reduced$seed_length,data_frameM_reduced$type, sep=";")

gap_pen28M <- c(which(matM.mrl[,"gap_pen"]==0.2),which(matM.mrl[,"gap_pen"]==0.8))
matM.mrl_reduced <- matM.mrl[-gap_pen28M,]
data_frameM.mrl_reduced <- as.data.frame(matM.mrl_reduced)
data_frameM.mrl_reduced$length_type <- paste(data_frameM.mrl_reduced$seed_length,data_frameM.mrl_reduced$type, sep=";")

library(tikzDevice)
library(car)


tikz("gap_pen_PR_Mmrl.tex",width=4.5,height=3.5)
scatterplot(mean_prec~gap_pen|length_type,data=data_frameM.mrl,smoother=FALSE,main="MiSeq,MRL",xlab = "gap penalty",ylab="mean precision",legend.title = "seedlengte;type")
dev.off()
tikz("gap_pen_SE_Mmrl.tex",width=4.5,height=3.5)
scatterplot(mean_sens~gap_pen|length_type,data=data_frameM.mrl,smoother=FALSE,main="MiSeq,MRL",xlab = "gap penalty",ylab="mean sensitivity",legend.title = "seedlengte;type")
dev.off()
tikz("gap_pen_PR_Hmrl.tex",width=4.5,height=3.5)
scatterplot(mean_prec~gap_pen|length_type,data=data_frame.mrl,smoother=FALSE,main="HiSeq,MRL",xlab = "gap penalty",ylab="mean precision",legend.title = "seedlengte;type")
dev.off()
tikz("gap_pen_SE_Hmrl.tex",width=4.5,height=3.5)
scatterplot(mean_sens~gap_pen|length_type,data=data_frame.mrl,smoother=FALSE,main="HiSeq,MRL",xlab = "gap penalty",ylab="mean sensitivity",legend.title = "seedlengte;type")
dev.off()
tikz("gap_pen_PR_M.tex",width=4.5,height=3.5)
scatterplot(mean_prec~gap_pen|length_type,data=data_frameM,smoother=FALSE,main="MiSeq,LCA*",xlab = "gap penalty",ylab="mean precision",legend.title = "seedlengte;type")
dev.off()
tikz("gap_pen_SE_M.tex",width=4.5,height=3.5)
scatterplot(mean_sens~gap_pen|length_type,data=data_frameM,smoother=FALSE,main="MiSeq,LCA*",xlab = "gap penalty",ylab="mean sensitivity",legend.title = "seedlengte;type")
dev.off()
tikz("gap_pen_PR_H.tex",width=4.5,height=3.5)
scatterplot(mean_prec~gap_pen|length_type,data=data_frame,smoother=FALSE,main="HiSeq,LCA*",xlab = "gap penalty",ylab="mean precision",legend.title = "seedlengte;type")
dev.off()
tikz("gap_pen_SE_H.tex",width=4.5,height=3.5)
scatterplot(mean_sens~gap_pen|length_type,data=data_frame,smoother=FALSE,main="HiSeq,LCA*",xlab = "gap penalty",ylab="mean sensitivity",legend.title = "seedlengte;type")
dev.off()

tikz("H_LCA_overview.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.5,0.85),ylim=c(0.8,0.99),
            labels = mat_reduced[,"gap_size"],id.n = nrow(data_frame_reduced),id.cex = 0.6,
            main="HiSeq,LCA*",xlab="genus sensitivity", ylab="genus precision",legend.title = "seedlengte;type")
dev.off()
tikz("H_MRL_overview.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame.mrl_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.5,0.85),ylim=c(0.8,0.99),
            main="HiSeq,MRL",xlab="genus sensitivity", ylab="genus precision",legend.title = "seedlengte;type")
dev.off()
tikz("M_LCA_overview.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.30,0.70),ylim=c(0.6,0.96),
            labels = matM_reduced[,"gap_size"],id.n = nrow(data_frameM_reduced),id.cex=0.6,
            main="MiSeq,LCA*",xlab="genus sensitivity", ylab="genus precision",legend.title = "seedlengte;type")
dev.off()
tikz("M_MRL_overview.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM.mrl_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.30,0.70),ylim=c(0.6,0.96),
            main="MiSeq,MRL",xlab="genus sensitivity", ylab="genus precision",legend.title = "seedlengte;type")
dev.off()

tikz("H_ID_vs_length.tex", width=4.5,height=3)
scatterplot(IDs~seed_length, data = data_frame, smoother=FALSE,main="HiSeq",xlab="seedlengte", ylab="aantal identificaties")
dev.off()
tikz("M_ID_vs_length.tex", width=4.5,height=3)
scatterplot(IDs~seed_length, data = data_frameM, smoother=FALSE,main="MiSeq",xlab="seedlengte", ylab="aantal identificaties")
dev.off


tikz("H_LCA_species.tex",width=4.5,height=3.5)
scatterplot(species_prec ~ species_sens | length_type,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.25,0.45),ylim=c(0.65,0.96),
            labels = mat_reduced[,"gap_size"],id.n = nrow(data_frame_reduced),id.cex = 0.6,
            main="HiSeq,LCA*",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()
tikz("H_MRL_species.tex",width=4.5,height=3.5)
scatterplot(species_prec ~ species_sens | length_type,data=data_frame.mrl_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.25,0.45),ylim=c(0.65,0.96),
            main="HiSeq,MRL",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()
tikz("M_LCA_species.tex",width=4.5,height=3.5)
scatterplot(species_prec ~ species_sens | length_type,data=data_frameM_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.15,0.45),ylim=c(0.45,0.8),
            labels = matM_reduced[,"gap_size"],id.n = nrow(data_frameM_reduced),id.cex=0.6,
            main="MiSeq,LCA*",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()
tikz("M_MRL_species.tex",width=4.5,height=3.5)
scatterplot(species_prec ~ species_sens | length_type,data=data_frameM.mrl_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.15,0.45),ylim=c(0.45,0.8),
            main="MiSeq,MRL",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()


tikz("H_LCA_species.tex",width=4.5,height=3.5)
scatterplot(phylum_prec ~ phylum_sens | length_type,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,
            labels = mat_reduced[,"gap_size"],id.n = nrow(data_frame_reduced),id.cex = 0.6,
            main="HiSeq,LCA*",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()
tikz("H_MRL_species.tex",width=4.5,height=3.5)
scatterplot(species_prec ~ species_sens | length_type,data=data_frame.mrl_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.25,0.45),ylim=c(0.65,0.96),
            main="HiSeq,MRL",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()
tikz("M_LCA_species.tex",width=4.5,height=3.5)
scatterplot(species_prec ~ species_sens | length_type,data=data_frameM_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.15,0.45),ylim=c(0.45,0.8),
            labels = matM_reduced[,"gap_size"],id.n = nrow(data_frameM_reduced),id.cex=0.6,
            main="MiSeq,LCA*",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()
tikz("M_MRL_species.tex",width=4.5,height=3.5)
scatterplot(species_prec ~ species_sens | length_type,data=data_frameM.mrl_reduced, smoother=FALSE,reg.line=FALSE, 
            legend.coords="bottomleft",legend.columns=3,xlim=c(0.15,0.45),ylim=c(0.45,0.8),
            main="MiSeq,MRL",xlab="species sensitivity", ylab="species precision",legend.title = "seedlengte;type")
dev.off()


tikz("M_SP_PR_Hyb08.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM08.mrl, smoother=FALSE,reg.line=FALSE,
            xlim=c(0.45,0.67),ylim=c(0.85,0.96),legend.coords="bottomleft",legend.columns=3,
            main="MiSeq,hybrid (0.8)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("M_SP_PR_Hyb02.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM02.mrl, smoother=FALSE,reg.line=FALSE
            ,xlim=c(0.45,0.67),ylim=c(0.85,0.96), legend.coords="bottomleft",legend.columns=3,
            main="MiSeq,hybrid (0.2)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("M_SP_PR_Hyb05.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM05.mrl, smoother=FALSE,reg.line=FALSE,
            xlim=c(0.45,0.67),ylim=c(0.85,0.96), legend.coords="bottomleft",legend.columns=3,
            main="MiSeq,hybrid (0.5)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("M_SP_PR_Hyb04.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM04.mrl, smoother=FALSE,reg.line=FALSE,
            xlim=c(0.45,0.67),ylim=c(0.85,0.96), legend.coords="bottomleft",legend.columns=3,
            main="MiSeq,hybrid (0.4)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("M_SP_PR_Hyb06.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM06.mrl, smoother=FALSE,reg.line=FALSE,
            xlim=c(0.45,0.67),ylim=c(0.85,0.96), legend.coords="bottomleft",legend.columns=3,
            main="MiSeq,hybrid (0.6)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("M_SP_PR_Hyb07.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frameM07.mrl, smoother=FALSE,reg.line=FALSE,
            xlim=c(0.45,0.67),ylim=c(0.85,0.96), legend.coords="bottomleft",legend.columns=3,
            main="MiSeq,hybrid (0.7)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()


tikz("H_SP_PR_Hyb08.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame08.mrl, smoother=FALSE,reg.line=FALSE,
            xlim=c(0.62,0.81),ylim=c(0.92,0.99), legend.coords="bottomleft",legend.columns=3,
            main="HiSeq,hybrid (0.8)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("H_SP_PR_Hyb02.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame02.mrl, smoother=FALSE,reg.line=FALSE,
            xlim=c(0.62,0.81),ylim=c(0.92,0.99), legend.coords="bottomleft",legend.columns=3,
            main="HiSeq,hybrid (0.2)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("H_SP_PR_Hyb05.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame05.mrl, smoother=FALSE,reg.line=FALSE, 
            xlim=c(0.62,0.81),ylim=c(0.92,0.99), legend.coords="bottomleft",legend.columns=3,
            main="HiSeq,hybrid (0.5)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("H_SP_PR_Hyb04.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame04.mrl, smoother=FALSE,reg.line=FALSE, 
            xlim=c(0.62,0.81),ylim=c(0.92,0.99), legend.coords="bottomleft",legend.columns=3,
            main="HiSeq,hybrid (0.4)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()
tikz("H_SP_PR_Hyb06.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame06.mrl, smoother=FALSE,reg.line=FALSE, 
            xlim=c(0.62,0.81),ylim=c(0.92,0.99), legend.coords="bottomleft",legend.columns=3,
            main="HiSeq,hybrid (0.6)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()

tikz("H_SP_PR_Hyb07.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame07.mrl, smoother=FALSE,reg.line=FALSE, 
            xlim=c(0.62,0.81),ylim=c(0.92,0.99), legend.coords="bottomleft",legend.columns=3,
            main="HiSeq,hybrid (0.7)",xlab="genus sensitiviteit", ylab="genus precisie",legend.title = "seedlengte;type")
dev.off()









CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")
tot.reads = 10000
grades <- c("superkingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name")
header <- c("superkingdom","phylum","class","order","family","genus","species")

all.sensitivity <- numeric()
all.specificity <- numeric()
identifications <- numeric()

k <- 18
filename <- paste0("Kraken_HiSeq.11.Opt",k,".csv")
reads <- read.csv(filename, header=TRUE)
true.lineages <- read.csv("../Kraken_hiseq_taxonomy.csv",fill = TRUE, header=TRUE)
rownames(true.lineages)<- true.lineages[,"header"]
reads[,"taxon_rank"] <- ordered(reads[,"taxon_rank"], levels = CLASSES)

correctness <- matrix(FALSE,nrow = nrow(reads), ncol = 7)
colnames(correctness) = grades
identifications <- c(identifications,nrow(reads))

for(i in 1:nrow(reads)){
  for(j in 1:length(grades)){
    correctness[i,j] <- as.character(reads[i,grades[j]]) == as.character(true.lineages[reads[i,"correct"],grades[j]])
  }
}

correct.rank <- numeric(0)
nr.rank <- numeric(0)

correct.rank <- c(correct.rank, length(intersect(which(reads[,"taxon_rank"]>="superkingdom"),which(correctness[,"superkingdom_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="superkingdom",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="phylum"),which(correctness[,"phylum_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="phylum",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="class"),which(correctness[,"class_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="class",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="order"),which(correctness[,"order_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="order",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="family"),which(correctness[,"family_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="family",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="genus"),which(correctness[,"genus_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="genus",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="species",na.rm=TRUE))


specificity <- c(correct.rank/nr.rank,mean(correct.rank/nr.rank))
sensitivity <- c(correct.rank/(correct.rank+(tot.reads-nr.rank)),mean(correct.rank/(correct.rank+(tot.reads-nr.rank))))
all.sensitivity <- rbind(all.sensitivity,sensitivity)
all.specificity <- rbind(all.specificity,specificity)

corr_kingdom <- NULL
corr_phylum <- NULL
corr_class <- NULL
corr_order <- NULL
corr_family <- NULL
corr_genus <- NULL
corr_species <- NULL
for(x in levels(reads[,"correct"])){
  corr_kingdom <- c(corr_kingdom,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="superkingdom"),which(correctness[,"superkingdom_name"])))))
  corr_phylum <- c(corr_phylum,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="phylum"),which(correctness[,"phylum_name"])))))
  corr_class <- c(corr_class,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="class"),which(correctness[,"class_name"])))))
  corr_order <- c(corr_order,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="order"),which(correctness[,"order_name"])))))
  corr_family <- c(corr_family,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="family"),which(correctness[,"family_name"])))))
  corr_genus <- c(corr_genus,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="genus"),which(correctness[,"genus_name"])))))
  corr_species <- c(corr_species,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"])))))
}
correct_mat <- rbind(corr_kingdom,corr_phylum,corr_class,corr_order,corr_family,corr_genus,corr_species)
colnames(correct_mat) <- levels(reads[,"correct"])
correct_dataframe  <- as.data.frame(correct_mat)

# #tikz("H_ids_perRang.tex",width=4.5,height=3.5)
# plot(correct_dataframe$A_hydrophila,type='n',ylim=c(0,1000),xaxt='n',main="HiSeq",
#      ylab = "aantal correcte identificaties", xlab = "taxonrang")
# for(i in 1:10){
#   lines(correct_dataframe[i], col=rainbow(10)[i])
# }
# axis(1, at=1:7, labels=header) 
# legend("bottomleft",legend=c("A. hydrophyla","B. cereus","B. fragilis","M. abscessus","P. fermentans","R. sphaeriodes","S. aureus","S. pneumoniaer","V.cholerae","X. axonopodis"),
#        col=rainbow(10),fill=rainbow(10),ncol=2,cex=0.5)
# #dev.off()


nr_kingdom <- NULL
nr_phylum <- NULL
nr_class <- NULL
nr_order <- NULL
nr_family <- NULL
nr_genus <- NULL
nr_species <- NULL
for(x in levels(reads[,"correct"])){
  nr_kingdom <- c(nr_kingdom,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="superkingdom"))))
  nr_phylum <- c(nr_phylum,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="phylum"))))
  nr_class <- c(nr_class,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="class"))))
  nr_order <- c(nr_order,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="order"))))
  nr_family <- c(nr_family,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="family"))))
  nr_genus <- c(nr_genus,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="genus"))))
  nr_species <- c(nr_species,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="species"))))
}
nr_mat <- rbind(nr_kingdom,nr_phylum,nr_class,nr_order,nr_family,nr_genus,nr_species)
colnames(nr_mat) <- levels(reads[,"correct"])
nr_dataframe  <- as.data.frame(nr_mat)

tikz("H_LCA_ids_perRang.tex",width=4.5,height=3.5)
plot(nr_dataframe$A_hydrophila,type='n',ylim=c(0,1000),xaxt='n',main="HiSeq,LCA*",
     ylab = "aantal correcte identificaties", xlab = "taxonrang")
for(i in 1:10){
  lines(nr_dataframe[i], col=rainbow(10)[i])
  lines(correct_dataframe[i], col = rainbow(10)[i],lty=3)
}
axis(1, at=1:7, labels=header) 
legend("bottomleft",legend=c("A. hydrophyla","B. cereus","B. fragilis","M. abscessus","P. fermentans","R. sphaeriodes","S. aureus","S. pneumoniaer","V.cholerae","X. axonopodis"),
       col=rainbow(10),fill=rainbow(10),ncol=2,cex=0.5)
dev.off()




all.sensitivity <- numeric()
all.specificity <- numeric()
identifications <- numeric()

k <- 18
filename <- paste0("Kraken_MiSeq.11.Opt",k,".csv")
reads <- read.csv(filename, header=TRUE)
true.lineages <- read.csv("../Kraken_miseq_taxonomy.csv",fill = TRUE, header=TRUE)
rownames(true.lineages)<- true.lineages[,"header"]
reads[,"taxon_rank"] <- ordered(reads[,"taxon_rank"], levels = CLASSES)

correctness <- matrix(FALSE,nrow = nrow(reads), ncol = 7)
colnames(correctness) = grades
identifications <- c(identifications,nrow(reads))

for(i in 1:nrow(reads)){
  for(j in 1:length(grades)){
    correctness[i,j] <- as.character(reads[i,grades[j]]) == as.character(true.lineages[reads[i,"correct"],grades[j]])
  }
}

correct.rank <- numeric(0)
nr.rank <- numeric(0)

correct.rank <- c(correct.rank, length(intersect(which(reads[,"taxon_rank"]>="superkingdom"),which(correctness[,"superkingdom_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="superkingdom",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="phylum"),which(correctness[,"phylum_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="phylum",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="class"),which(correctness[,"class_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="class",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="order"),which(correctness[,"order_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="order",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="family"),which(correctness[,"family_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="family",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="genus"),which(correctness[,"genus_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="genus",na.rm=TRUE))

correct.rank <- c(correct.rank,length(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]))))
nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="species",na.rm=TRUE))


specificity <- c(correct.rank/nr.rank,mean(correct.rank/nr.rank))
sensitivity <- c(correct.rank/(correct.rank+(tot.reads-nr.rank)),mean(correct.rank/(correct.rank+(tot.reads-nr.rank))))
all.sensitivity <- rbind(all.sensitivity,sensitivity)
all.specificity <- rbind(all.specificity,specificity)

sum(reads[,"species_name"]=="Proteus vulgaris")

corr_kingdom <- NULL
corr_phylum <- NULL
corr_class <- NULL
corr_order <- NULL
corr_family <- NULL
corr_genus <- NULL
corr_species <- NULL
for(x in levels(reads[,"correct"])){
  corr_kingdom <- c(corr_kingdom,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="superkingdom"),which(correctness[,"superkingdom_name"])))))
  corr_phylum <- c(corr_phylum,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="phylum"),which(correctness[,"phylum_name"])))))
  corr_class <- c(corr_class,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="class"),which(correctness[,"class_name"])))))
  corr_order <- c(corr_order,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="order"),which(correctness[,"order_name"])))))
  corr_family <- c(corr_family,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="family"),which(correctness[,"family_name"])))))
  corr_genus <- c(corr_genus,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="genus"),which(correctness[,"genus_name"])))))
  corr_species <- c(corr_species,length(intersect(which(reads[,"correct"]==x),intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"])))))
}
correct_mat <- rbind(corr_kingdom,corr_phylum,corr_class,corr_order,corr_family,corr_genus,corr_species)
colnames(correct_mat) <- levels(reads[,"correct"])
correct_dataframe  <- as.data.frame(correct_mat)

# tikz("M_ids_perRang.tex",width=4.5,height=3.5)
# plot(correct_dataframe$V_cholerae,type='n',ylim=c(0,1000),xaxt='n',main="MiSeq",
#      ylab = "aantal correcte identificaties",xlab = "taxonrang")
# for(i in 1:10){
#   lines(correct_dataframe[i], col=rainbow(10)[i])
# }
# axis(1, at=1:7, labels=header) 
# legend("bottomleft",legend=c("B. cereus","C. freundii","E. cloacae","K. pneumoniae","M. abscessus","P. vulgaris","R. sphaeroides","S. aureus","S. enterica","V. cholerae"),
#        col=rainbow(10),fill=rainbow(10),ncol=2,cex=0.5)
# dev.off()




nr_kingdom <- NULL
nr_phylum <- NULL
nr_class <- NULL
nr_order <- NULL
nr_family <- NULL
nr_genus <- NULL
nr_species <- NULL
for(x in levels(reads[,"correct"])){
  nr_kingdom <- c(nr_kingdom,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="superkingdom"))))
  nr_phylum <- c(nr_phylum,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="phylum"))))
  nr_class <- c(nr_class,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="class"))))
  nr_order <- c(nr_order,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="order"))))
  nr_family <- c(nr_family,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="family"))))
  nr_genus <- c(nr_genus,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="genus"))))
  nr_species <- c(nr_species,length(intersect(which(reads[,"correct"]==x),which(reads[,"taxon_rank"]>="species"))))
}
nr_mat <- rbind(nr_kingdom,nr_phylum,nr_class,nr_order,nr_family,nr_genus,nr_species)
colnames(nr_mat) <- levels(reads[,"correct"])
nr_dataframe  <- as.data.frame(nr_mat)

tikz("M_LCA_ids_perRang.tex",width=4.5,height=3.5)
plot(nr_dataframe$V_cholerae,type='n',ylim=c(0,1000),xaxt='n',main="MiSeq,LCA*",
     ylab = "aantal identificaties",xlab = "taxonrang")
for(i in 1:10){
  lines(nr_dataframe[i], col=rainbow(10)[i])
  lines(correct_dataframe[i], col = rainbow(10)[i], lty=3)
}
axis(1, at=1:7, labels=header) 
legend("bottomleft",legend=c("B. cereus","C. freundii","E. cloacae","K. pneumoniae","M. abscessus","P. vulgaris","R. sphaeroides","S. aureus","S. enterica","V. cholerae"),
       col=rainbow(10),fill=rainbow(10),ncol=2,cex=0.5)
dev.off()

