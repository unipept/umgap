library(tikzDevice)
library(car)

setwd("/Users/Aranka/Documents/Unief/Thesis/Thesis/Kmer2consensus/metabenchmark")

mat <- read.csv("Metabench_analyse.LCA.csv")
mat.mrl <- read.csv("Metabench_analyse.MRL.csv")
mat.hyb <- read.csv("Metabench_analyse.Hybrid06.csv")
mat.hyb2 <- read.csv("Metabench_analyse.Hybrid04.csv")
data_frame <- as.data.frame(mat)
data_frame.mrl <- as.data.frame(mat.mrl)
data_frame.hyb <- as.data.frame(mat.hyb)
data_frame.hyb2 <- as.data.frame(mat.hyb2)
data_frame$length_type <- paste(data_frame$seed_length,data_frame$score, sep=";")
data_frame.mrl$length_type <- paste(data_frame.mrl$seed_length,data_frame.mrl$score, sep=";")
data_frame.hyb$length_type <- paste(data_frame.hyb$seed_length,data_frame.hyb$score, sep=";")
data_frame.hyb2$length_type <- paste(data_frame.hyb2$seed_length,data_frame.hyb2$score, sep=";")


### Specificiteit en NPV
tikz("metabench_LCA_SP_NPV.tex",width=4.5,height=4)
scatterplot(genus_spec ~ genus_npv | length_type,data=data_frame, smoother=FALSE,reg.line=FALSE,
            labels = mat[,"gap_size"],id.n = nrow(data_frame),id.cex=0.8,legend.coords="bottomright",
            legend.columns=2,main="LCA*",xlab = "genus NPV",ylab="genus specificiteit",legend.title = "seedlengte;score",
            xlim=c(0.15,0.55), ylim=c(0.42,0.92))
dev.off()
tikz("metabench_MRL_SP_NPV.tex",width=4.5,height=4)
scatterplot(genus_spec ~ genus_npv | length_type,data=data_frame.mrl, smoother=FALSE,reg.line=FALSE,
            labels = mat[,"gap_size"],id.n = nrow(data_frame),id.cex=0.8,legend.coords="bottomleft",
            legend.columns=2,main="MRL",xlab = "genus NPV",ylab="genus specificiteit",legend.title = "seedlengte;score",
            xlim=c(0.15,0.55), ylim=c(0.42,0.92))
dev.off()
tikz("metabench_Hyb06_SP_NPV.tex",width=4.5,height=4)
scatterplot(genus_spec ~ genus_npv | length_type,data=data_frame.hyb, smoother=FALSE,reg.line=FALSE,
            labels = mat[,"gap_size"],id.n = nrow(data_frame),id.cex=0.8,legend.coords="topleft",
            legend.columns=2,main="Hybride 0.6",xlab = "genus NPV",ylab="genus specificiteit",legend.title = "seedlengte;score",
            xlim=c(0.15,0.55), ylim=c(0.42,0.92))
dev.off()
tikz("metabench_Hyb04_SP_NPV.tex",width=4.5,height=4)
scatterplot(genus_spec ~ genus_npv | length_type,data=data_frame.hyb2, smoother=FALSE,reg.line=FALSE,
            labels = mat[,"gap_size"],id.n = nrow(data_frame),id.cex=0.8,legend.coords="topleft",
            legend.columns=2,main="Hybride 0.4",xlab = "genus NPV",ylab="genus specificiteit",legend.title = "seedlengte;score",
            xlim=c(0.15,0.55), ylim=c(0.42,0.92))
dev.off()

###### Precisie en sensitiviteit

tikz("metabench_LCA_PR_SE_abstr.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame, smoother=FALSE,reg.line=FALSE,
            labels = mat[,"gap_size"],id.n = nrow(data_frame),id.cex = 0.8, legend.coords="topleft",
            legend.columns=2,main="LCA*",xlab = "genus sensitivity",ylab="genus precision",
            legend.title = "seedlength;score", xlim=c(0.2,0.87),ylim=c(0.75,0.97))
dev.off()
tikz("metabench_MRL_PR_SE_abstr.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame.mrl, smoother=FALSE,reg.line=FALSE,
            labels = mat.mrl[,"gap_size"],id.n = nrow(data_frame.mrl),id.cex = 0.8, legend.coords="topleft",
            legend.columns=2,main="MRL",xlab = "genus sensitivity",ylab="genus precision",
            legend.title = "seedlength;score", xlim=c(0.2,0.87),ylim=c(0.75,0.97))
dev.off()
tikz("metabench_Hyb06_PR_SE.tex",width=4.5,height=4)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame.hyb, smoother=FALSE,reg.line=FALSE,
            labels = mat[,"gap_size"],id.n = nrow(data_frame),id.cex = 0.8, legend.coords="topleft",
            legend.columns=2,main="Hybride 0.6",xlab = "genus sensitiviteit",ylab="genus precisie",
            legend.title = "seedlengte;score", xlim=c(0.2,0.87),ylim=c(0.75,0.97))
dev.off()
tikz("metabench_Hyb04_PR_SE_abstr.tex",width=4.5,height=3.5)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame.hyb2, smoother=FALSE,reg.line=FALSE,
            labels = mat[,"gap_size"],id.n = nrow(data_frame),id.cex = 0.8, legend.coords="topleft",
            legend.columns=2,main="Hybrid 0.4",xlab = "genus sensitivity",ylab="genus precision",
            legend.title = "seedlength;score", xlim=c(0.2,0.87),ylim=c(0.75,0.97))
dev.off()


s#### barplots

data_genus <- t(as.matrix(mat.mrl[c(2,6,11),c("genus_sens","genus_spec","genus_prec","genus_npv","genus_mcc")]))
data_phylum <- t(as.matrix(mat.mrl[c(2,6,11),c("phylum_sens","phylum_spec","phylum_prec","phylum_npv","phylum_mcc")]))

colnames(data_genus) <- c("1","2","3")
colnames(data_phylum)<- c("1","2","3")
barplot(data_genus,beside=TRUE,horiz = TRUE,xlim=c(0,1))
barplot(data_phylum,beside=TRUE,horiz = TRUE,xlim=c(0,1))

