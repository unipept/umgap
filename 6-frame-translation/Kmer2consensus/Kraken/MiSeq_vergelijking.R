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

for(k in 1:51){
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
  print(k)
}

tot_perSP <- 1000
IDs_perSP <- summary(reads[,"correct"])
perc_IDs_perSP <- IDs_perSP/tot_perSP

#### Zowel genus als species fout
beiden_fout <- reads[intersect(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)),which(correctness[,"genus_name"]==F)),c("correct","species_name","genus_name")]
summary(beiden_fout[,"correct"])/IDs_perSP
summary(beiden_fout[beiden_fout[,"correct"]=="C_freundii",])
summary(beiden_fout[beiden_fout[,"correct"]=="B_cereus",])
summary(beiden_fout[beiden_fout[,"correct"]=="E_cloacae",])
summary(beiden_fout[beiden_fout[,"correct"]=="K_pneumoniae",])
summary(beiden_fout[beiden_fout[,"correct"]=="P_vulgaris",])
summary(beiden_fout[beiden_fout[,"correct"]=="S_aureus",])
summary(beiden_fout[beiden_fout[,"correct"]=="S_enterica",])
#### Genus juist, species fout
genus_juist <- reads[intersect(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)),which(correctness[,"genus_name"])),c("correct","species_name","genus_name")]
summary(genus_juist[,"correct"])/IDs_perSP
summary(genus_juist[genus_juist[,"correct"]=="B_cereus",])
summary(genus_juist[genus_juist[,"correct"]=="E_cloacae",])
summary(genus_juist[genus_juist[,"correct"]=="K_pneumoniae",])
summary(genus_juist[genus_juist[,"correct"]=="M_abscessus",])
summary(genus_juist[genus_juist[,"correct"]=="P_vulgaris",])
summary(genus_juist[genus_juist[,"correct"]=="S_aureus",])

#### Percentage van foute species met correct genus
length(intersect(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)),which(correctness[,"genus_name"])))/length(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)))


model_params <- matrix(c(3,2,1,0.5,
                         3,2,2,NA,
                         3,2,3,0.5,
                         3,2,1,0.2,
                         3,2,3,0.2,
                         3,2,1,0.8,
                         3,2,3,0.8,
                         3,0,1,NA,
                         3,0,2,NA,
                         3,0,3,NA,
                         3,4,1,0.5,
                         3,4,2,NA,
                         3,4,3,0.5,
                         3,4,1,0.2,
                         3,4,3,0.2,
                         3,4,1,0.8,
                         3,4,3,0.8,
                         2,2,1,0.5,
                         2,2,2,NA,
                         2,2,3,0.5,
                         2,2,1,0.2,
                         2,2,3,0.2,
                         2,2,1,0.8,
                         2,2,3,0.8,
                         2,0,1,NA,
                         2,0,2,NA,
                         2,0,3,NA,
                         2,4,1,0.5,
                         2,4,2,NA,
                         2,4,3,0.5,
                         2,4,1,0.2,
                         2,4,3,0.2,
                         2,4,1,0.8,
                         2,4,3,0.8,
                         4,2,1,0.5,
                         4,2,2,NA,
                         4,2,3,0.5,
                         4,2,1,0.2,
                         4,2,3,0.2,
                         4,2,1,0.8,
                         4,2,3,0.8,
                         4,0,1,NA,
                         4,0,2,NA,
                         4,0,3,NA,
                         4,4,1,0.5,
                         4,4,2,NA,
                         4,4,3,0.5,
                         4,4,1,0.2,
                         4,4,3,0.2,
                         4,4,1,0.8,
                         4,4,3,0.8), ncol=4,byrow=TRUE)

mat <- cbind(model_params,identifications/tot.reads, all.sensitivity, all.specificity)
colnames(mat) <- c("seed_length","gap_size","type","gap_pen","IDs","superkingdom_sens","phylum_sens","class_sens","order_sens","family_sens","genus_sens","species_sens","mean_sens","superkingdom_prec","phylum_prec","class_prec","order_prec","family_prec","genus_prec","species_prec","mean_prec")
rownames(mat) <- c(1:51)

write.csv(mat,file="MiSeq_PR_SE.csv",sep=",")

library(car) 

data_frame <- as.data.frame(mat)
data_frame$length_type <- paste(data_frame$seed_length,data_frame$type, sep=";")

scatterplot(mean_prec ~ mean_sens | seed_length,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | type,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | gap_size,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | gap_pen,data=data_frame, smoother=FALSE,reg.line=FALSE)

scatterplot(mean_prec ~ mean_sens | length_type,data=data_frame, smoother=FALSE,reg.line=FALSE)

##############################
## SUBSET ####################
##############################

gap_pen28 <- c(which(mat[,"gap_pen"]==0.2),which(mat[,"gap_pen"]==0.8))
mat_reduced <- mat[-gap_pen28,]
data_frame_reduced <- as.data.frame(mat_reduced)
data_frame_reduced$length_type <- paste(data_frame_reduced$seed_length,data_frame_reduced$type, sep=";")

scatterplot(mean_prec ~ mean_sens | type,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | seed_length,data=data_frame_reduced, 
            smoother=FALSE,reg.line=FALSE,labels = data_frame_reduced[,"type"],id.n = nrow(data_frame_reduced),
            id.cex = 0.7,xlim=c(0.59,0.76),ylim=c(0.9,0.99))

scatterplot(seed_length ~ IDs | type, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(type ~ IDs | seed_length, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)

scatterplot(genus_prec ~ genus_sens | seed_length,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE,labels = rownames(mat_reduced),id.n = nrow(data_frame_reduced),
            id.cex = 0.7,xlim=c(0.45,0.6),ylim=c(0.87,0.96))
scatterplot(species_prec ~ species_sens | seed_length,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE,labels = rownames(mat_reduced),id.n = nrow(data_frame_reduced),
            id.cex = 0.7,xlim=c(0.24,0.34),ylim=c(0.68,0.78))

scatterplot(mean_prec ~IDs | seed_length, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_sens ~IDs | seed_length, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)

scatterplot(mean_prec ~ mean_sens | length_type,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | length_type,data=data_frame_reduced,smoother=FALSE,reg.line=FALSE,xlim=c(0.58,0.76),ylim=c(0.90,0.96))

scatterplot(genus_prec ~ genus_sens | length_type,data=data_frame_reduced, 
            smoother=FALSE,reg.line=FALSE,xlim=c(0.45,0.6),ylim=c(0.87,0.96))
scatterplot(species_prec ~ species_sens |length_type,data=data_frame_reduced,
            smoother=FALSE,reg.line=FALSE,xlim=c(0.24,0.34),ylim=c(0.68,0.78))
