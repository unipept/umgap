setwd("/Users/Aranka/Documents/Unief/Thesis/Thesis/Kmer2consensus/metabenchmark")

CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")
tot.reads = 27658749
tot.random = 4970853
grades <- c("superkingdom_id","phylum_id","class_id","order_id","family_id","genus_id","species_id")

all.sensitivity <- numeric()
all.npv <- numeric()
all.precision <- numeric()
all.specificity <- numeric()
all.mcc <- numeric()
identifications <- numeric()

for(k in c(1,8,11,18,19,25,26,28,29,36,43,46)){
  filename <- paste0("setA1.",k,".Hybrid04.counts")
  reads <- read.csv(filename, header=FALSE)
  colnames(reads) <- c("count","phylum", "embl","taxon_rank","taxon_id","superkingdom_id","kingdom_id","subkingdom_id","superphylum_id","phylum_id","subphylum_id","superclass_id","class_id","subclass_id","infraclass_id","superorder_id","order_id","suborder_id","infraorder_id","parvorder_id","superfamily_id","family_id","subfamily_id","tribe_id","subtribe_id","genus_id","subgenus_id","species_group_id","species_subgroup_id","species_id","subspecies_id","varietas_id","forma_id")
  true.lineages <- read.csv("present_lineages.csv",fill = TRUE, header=FALSE)
  colnames(true.lineages) <- c("phylum", "embl","taxon_id","superkingdom_id","kingdom_id","subkingdom_id","superphylum_id","phylum_id","subphylum_id","superclass_id","class_id","subclass_id","infraclass_id","superorder_id","order_id","suborder_id","infraorder_id","parvorder_id","superfamily_id","family_id","subfamily_id","tribe_id","subtribe_id","genus_id","subgenus_id","species_group_id","species_subgroup_id","species_id","subspecies_id","varietas_id","forma_id")
  #rownames(true.lineages)<- true.lineages[,"embl"]
  reads[,"taxon_rank"] <- ordered(reads[,"taxon_rank"], levels = CLASSES)
  
  random.mapped <- c(which(as.character(reads[,"embl"])=="SHUFFLED"),which(as.character(reads[,"phylum"])=="Random"))
  reads.reduced <- reads[-random.mapped,]
  
  correctness <- matrix(FALSE,nrow = nrow(reads.reduced), ncol = 7)
  colnames(correctness) = grades
  identifications <- c(identifications,sum(reads[,"count"]))
  
  for(i in 1:nrow(reads.reduced)){
    for(j in length(grades):1){
      juist <- which(as.character(true.lineages[,"embl"])==as.character(reads.reduced[i,"embl"]))[1]
      corr <- reads.reduced[i,grades[j]] == true.lineages[juist,grades[j]]
      if(is.na(corr)){
        correctness[i,j] <- correctness[i,j%%length(grades)+1]
      }else{
        correctness[i,j] <- corr
      }
    }
  }
  
  #random.mapped <- sum(reads[which(as.character(reads[,"embl"])=="SHUFFLED"),"count"]) + sum(reads[which(as.character(reads[,"phylum"])=="Random"),"count"])
  
  correct.rank <- numeric(0)
  nr.rank <- numeric(0)
  random.mapped.rank <- numeric(0)
  
  for(rank in header){
    random.mapped.rank <- c(random.mapped.rank,sum(reads[intersect(which(reads[,"taxon_rank"]>=rank),random.mapped),"count"]))
    correct.rank <- c(correct.rank,sum(reads.reduced[intersect(which(reads.reduced[,"taxon_rank"]>=rank),which(correctness[,paste0(rank,"_id")])),"count"]))
    nr.rank <- c(nr.rank,sum(reads.reduced[which(reads.reduced[,"taxon_rank"]>=rank),"count"],na.rm=TRUE))
  }

  TN <- tot.random - random.mapped.rank
  FP <- nr.rank - correct.rank + random.mapped.rank
  TP <- correct.rank
  FN <- tot.reads - TN - FP - TP
  precision <- c(TP/(TP+FP),mean(TP/(TP+FP)))
  specificity <- c(TN/(TN+FP),mean(TN/TN+FP))
  sensitivity <- c(TP/(TP+FN),mean(TP/(TP+FN)))
  NPV <- c(TN/(TN+FN),mean(TN/TN+FN))
  MCC <- c((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),mean((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  all.precision <- rbind(all.precision,precision)
  all.sensitivity <- rbind(all.sensitivity,sensitivity)
  all.specificity <- rbind(all.specificity,specificity)
  all.npv <- rbind(all.npv,NPV)
  all.mcc <- rbind(all.mcc,MCC)
  print(k)
}



model_params <- matrix(c(3,2,2,0.5,
                         3,0,2,NA,
                         3,4,2,0.5,
                         2,2,2,0.5,
                         2,2,1,NA,
                         2,0,2,NA,
                         2,0,1,NA,
                         2,4,2,0.5,
                         2,4,1,NA,
                         4,2,1,NA,
                         4,0,1,NA,
                         4,4,1,NA), ncol=4,byrow=TRUE)

mat <- cbind(model_params,identifications/tot.reads, all.sensitivity, all.precision,all.specificity,all.npv,all.mcc)
colnames(mat) <- c("seed_length","gap_size","score","gap_pen","IDs",
                   "superkingdom_sens","phylum_sens","class_sens","order_sens","family_sens","genus_sens","species_sens","mean_sens",
                   "superkingdom_prec","phylum_prec","class_prec","order_prec","family_prec","genus_prec","species_prec","mean_prec",
                   "superkingdom_spec","phylum_spec","class_spec","order_spec","family_spec","genus_spec","species_spec","mean_spec",
                   "superkingdom_npv","phylum_npv","class_npv","order_npv","family_npv","genus_npv","species_npv","mean_npv",
                   "superkingdom_mcc","phylum_mcc","class_mcc","order_mcc","family_mcc","genus_mcc","species_mcc","mean_mcc")
rownames(mat) <- c(1,8,11,18,19,25,26,28,29,36,43,46)

write.csv(mat,file="Metabench_analyse.Hybrid04.csv",sep=",")

