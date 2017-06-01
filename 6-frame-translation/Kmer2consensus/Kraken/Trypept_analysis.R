setwd("/Users/Aranka/Documents/Unief/Thesis/Thesis/Kmer2consensus/Kraken")

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

  filename <- paste0("Kraken_HiSeq.11.pept.MRL.csv")
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
  sensitivity <- c(correct.rank/tot.reads,mean(correct.rank/tot.reads))

tot_perSP <- 1000
IDs_perSP <- summary(reads[,"correct"])
perc_IDs_perSP <- IDs_perSP/tot_perSP


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
correct_mat <- rbind(NA,corr_kingdom,corr_phylum,corr_class,corr_order,corr_family,corr_genus,corr_species)
colnames(correct_mat) <- levels(reads[,"correct"])
correct_dataframe  <- as.data.frame(correct_mat)

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
nr_mat <- rbind(IDs_perSP,nr_kingdom,nr_phylum,nr_class,nr_order,nr_family,nr_genus,nr_species)
colnames(nr_mat) <- levels(reads[,"correct"])
nr_dataframe  <- as.data.frame(nr_mat)

library(tikzDevice)
tikz("HiSeq_trypept_ids_MRL.tex",width=5,height=3.5)
plot(nr_dataframe$A_hydrophila,type='n',ylim=c(-100,1000),xaxt='n',main="HiSeq,tryptische peptiden",
     ylab = "aantal correcte identificaties", xlab = "taxonrang")
for(i in 1:10){
  lines(nr_dataframe[i], col=rainbow(10)[i],lwd=1.5)
  lines(correct_dataframe[i], col = rainbow(10)[i],lty=3,lwd=1.5)
}
axis(1, at=1:8, labels=c("ids",header),cex.axis=0.1)
legend("bottomleft",legend=c("A. hydrophyla","B. cereus","B. fragilis","M. abscessus","P. fermentans","R. sphaeriodes","S. aureus","S. pneumoniae","V.cholerae","X. axonopodis"),
       col=rainbow(10),fill=rainbow(10),ncol=2,cex=0.63)
dev.off()

tikz("MiSeq_trypept_ids_MRL.tex",width=4.5,height=3.5)
plot(nr_dataframe$V_cholerae,type='n',ylim=c(-100,1000),xaxt='n',main="MiSeq,LCA*",
     ylab = "aantal identificaties",xlab = "taxonrang")
for(i in 1:10){
  lines(nr_dataframe[i], col=rainbow(10)[i],lwd=1.5)
  lines(correct_dataframe[i], col = rainbow(10)[i], lty=3,lwd=1.5)
}
axis(1, at=1:7, labels=c("ids",header) 
legend("bottomleft",legend=c("B. cereus","C. freundii","E. cloacae","K. pneumoniae","M. abscessus","P. vulgaris","R. sphaeroides","S. aureus","S. enterica","V. cholerae"),
       col=rainbow(10),fill=rainbow(10),ncol=2,cex=0.63)
dev.off()

