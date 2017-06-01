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

for(k in c(1,2,3,8,9,10,11,12,13,18,19,20,25,26,27,28,29,30,35,36,37,42,43,44,45,46,47)){
  filename <- paste0("Kraken_MiSeq.11.Opt",k,".bothTree07.csv")
  reads <- read.csv(filename, header=TRUE)
  true.lineages <- read.csv("../Kraken_miseq_taxonomy2.csv",fill = TRUE, header=TRUE)
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
  #sensitivity <- c(correct.rank/tot.reads,mean(correct.rank/tot.reads))
  sensitivity <- c(correct.rank/(correct.rank+(tot.reads-nr.rank)),mean(correct.rank/(correct.rank+(tot.reads-nr.rank))))
  all.sensitivity <- rbind(all.sensitivity,sensitivity)
  all.specificity <- rbind(all.specificity,specificity)
  print(k)
}

tot_perSP <- 1000
IDs_perSP <- summary(reads[,"correct"])
perc_IDs_perSP <- IDs_perSP/tot_perSP


model_params <- matrix(c(3,2,1,0.5,
                         3,2,2,NA,
                         3,2,3,0.5,
                         3,0,1,NA,
                         3,0,2,NA,
                         3,0,3,NA,
                         3,4,1,0.5,
                         3,4,2,NA,
                         3,4,3,0.5,
                         2,2,1,0.5,
                         2,2,2,NA,
                         2,2,3,0.5,
                         2,0,1,NA,
                         2,0,2,NA,
                         2,0,3,NA,
                         2,4,1,0.5,
                         2,4,2,NA,
                         2,4,3,0.5,
                         4,2,1,0.5,
                         4,2,2,NA,
                         4,2,3,0.5,
                         4,0,1,NA,
                         4,0,2,NA,
                         4,0,3,NA,
                         4,4,1,0.5,
                         4,4,2,NA,
                         4,4,3,0.5), ncol=4,byrow=TRUE)

mat <- cbind(model_params,identifications/tot.reads, all.sensitivity, all.specificity)
colnames(mat) <- c("seed_length","gap_size","type","gap_pen","IDs","superkingdom_sens","phylum_sens","class_sens","order_sens","family_sens","genus_sens","species_sens","mean_sens","superkingdom_prec","phylum_prec","class_prec","order_prec","family_prec","genus_prec","species_prec","mean_prec")
rownames(mat) <- c(1,2,3,8,9,10,11,12,13,18,19,20,25,26,27,28,29,30,35,36,37,42,43,44,45,46,47)

write.csv(mat,file="MiSeq_PRe_SE.bothTree07.csv",sep=",")

library(car) 
data_frame <- as.data.frame(mat)
summary(data_frame)
data_frame$length_type <- paste(data_frame$seed_length,data_frame$type, sep=";")

scatterplot(genus_prec ~ genus_sens | seed_length,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(genus_prec ~ genus_sens | type,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(IDs ~ seed_length | length_type, data=data_frame, smoother=FALSE,reg.line=FALSE)

gap_pen28 <- c(which(mat[,"gap_pen"]==0.2),which(mat[,"gap_pen"]==0.8))
mat_reduced <- mat[-gap_pen28,]
data_frame_reduced <- as.data.frame(mat_reduced)
data_frame_reduced$length_type <- paste(data_frame_reduced$seed_length,data_frame_reduced$type, sep=";")
