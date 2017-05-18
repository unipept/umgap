setwd("/Users/Aranka/Documents/Unief/Thesis/Thesis/Kmer2consensus/metabenchmark")

CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")
tot.reads = 10000
grades <- c("superkingdom_id","phylum_id","class_id","order_id","family_id","genus_id","species_id")
header <- c("superkingdom","phylum","class","order","family","genus","species")

all.sensitivity <- numeric()
all.specificity <- numeric()
identifications <- numeric()

for(k in 1:51){
  k=8
  filename <- paste0("setA1.",k,".LCA.counts")
  reads <- read.csv(filename, header=FALSE)
  colnames(reads) <- c("count","phylum", "embl","taxon_rank","taxon_id","superkingdom_id","kingdom_id","subkingdom_id","superphylum_id","phylum_id","subphylum_id","superclass_id","class_id","subclass_id","infraclass_id","superorder_id","order_id","suborder_id","infraorder_id","parvorder_id","superfamily_id","family_id","subfamily_id","tribe_id","subtribe_id","genus_id","subgenus_id","species_group_id","species_subgroup_id","species_id","subspecies_id","varietas_id","forma_id")
  true.lineages <- read.csv("present_lineage.csv",fill = TRUE, header=FALSE)
  colnames(true.lineages) <- c("phylum", "embl","taxon_id","superkingdom_id","kingdom_id","subkingdom_id","superphylum_id","phylum_id","subphylum_id","superclass_id","class_id","subclass_id","infraclass_id","superorder_id","order_id","suborder_id","infraorder_id","parvorder_id","superfamily_id","family_id","subfamily_id","tribe_id","subtribe_id","genus_id","subgenus_id","species_group_id","species_subgroup_id","species_id","subspecies_id","varietas_id","forma_id")
  rownames(true.lineages)<- true.lineages[,"embl"]
  #reads[,"taxon_rank"] <- ordered(reads[,"taxon_rank"], levels = CLASSES)
  
  correctness <- matrix(FALSE,nrow = nrow(reads), ncol = 7)
  colnames(correctness) = grades
  identifications <- c(identifications,sum(reads[,"count"]))
  
  for(i in 1:nrow(reads)){
    for(j in 1:length(grades)){
      correctness[i,j] <- as.character(reads[i,grades[j]]) == as.character(true.lineages[reads[i,"embl"],grades[j]])
    }
  }
  
  sum(na.omit(correctness[,"genus_id"]))
  
  correct.rank <- numeric(0)
  nr.rank <- numeric(0)
  
  correct.rank <- c(correct.rank,sum(reads[1,intersect(which(reads[,"taxon_rank"]>="superkingdom"),which(correctness[,"superkingdom_id"]))]))
  nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="superkingdom",na.rm=TRUE))
  
  correct.rank <- c(correct.rank,sum(reads[1,intersect(which(reads[,"taxon_rank"]>="phylum"),which(correctness[,"phylum_id"]))]))
  nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="phylum",na.rm=TRUE))
  
  correct.rank <- c(correct.rank,sum(reads[1,intersect(which(reads[,"taxon_rank"]>="class"),which(correctness[,"class_id"]))]))
  nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="class",na.rm=TRUE))
  
  correct.rank <- c(correct.rank,sum(reads[1,intersect(which(reads[,"taxon_rank"]>="order"),which(correctness[,"order_id"]))]))
  nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="order",na.rm=TRUE))
  
  correct.rank <- c(correct.rank,sum(reads[1,intersect(which(reads[,"taxon_rank"]>="family"),which(correctness[,"family_id"]))]))
  nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="family",na.rm=TRUE))
  
  correct.rank <- c(correct.rank,sum(reads[1,intersect(which(reads[,"taxon_rank"]>="genus"),which(correctness[,"genus_id"]))]))
  nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="genus",na.rm=TRUE))
  
  correct.rank <- c(correct.rank,sum(reads[1,intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_id"]))]))
  nr.rank <- c(nr.rank,sum(reads[,"taxon_rank"]>="species",na.rm=TRUE))
  
  
  specificity <- c(correct.rank/nr.rank,mean(correct.rank/nr.rank))
  sensitivity <- c(correct.rank/(correct.rank+(tot.reads-nr.rank)),mean(correct.rank/(correct.rank+(tot.reads-nr.rank))))
  all.sensitivity <- rbind(all.sensitivity,sensitivity)
  all.specificity <- rbind(all.specificity,specificity)
  print(k)
}