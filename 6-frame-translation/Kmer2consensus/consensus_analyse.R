tot.reads <- 10000
present.taxons <- read.csv("Kraken_HiSeq_presentTaxons", header=FALSE)
summary(present.taxons)
levels(present.taxons[,3])
true.lineages <- read.csv("Kraken_HiSeq_Lineages",fill = TRUE, header=FALSE)
reads.consensus <- read.csv("Kraken_HiSeq_readConsensus.csv", header=FALSE)
reads.consensus[,2] <- as.factor(reads.consensus[,2])
reads.consensus[,5] <- as.factor(reads.consensus[,5])
reads.consensus[,ncol(reads.consensus)+1] <- reads.consensus[,4]-reads.consensus[,3]
reads.consensus[,ncol(reads.consensus)+1] <- "null"
reads.consensus[,ncol(reads.consensus)+1] <- "null"
reads.consensus[,ncol(reads.consensus)+1] <- FALSE
for(i in 1:nrow(reads.consensus)){
  reads.consensus[i,7] <- as.character(present.taxons[which(present.taxons[,2]==reads.consensus[i,5]),3])
  reads.consensus[i,8] <- as.character(present.taxons[which(present.taxons[,2]==reads.consensus[i,5]),4])
  if(any(as.character(reads.consensus[i,7])==true.lineages)||as.character(reads.consensus[i,7])=="root"){
    reads.consensus[i,9] <- TRUE
  }else{
    reads.consensus[i,9] <- FALSE
  }
}
reads.consensus[,7] <- as.factor(reads.consensus[,7])
reads.consensus[,8] <- as.factor(reads.consensus[,8])
colnames(reads.consensus) <- c("Header","Frame","start","end","taxonID","n_kmers","taxonName","taxonRank","correct")
summary(reads.consensus)
CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")
reads.consensus[,"taxonRank"] <- ordered(reads.consensus[,"taxonRank"], levels = CLASSES)

count(reads.consensus[,"correct"])


Voorspellingen <- nrow(reads.consensus)

Genera <- sum(reads.consensus[,"taxonRank"]>="genus")
cor.genera <- sum(reads.consensus[reads.consensus[,"taxonRank"]>="genus","correct"]==TRUE)
incor.voorspelling <- sum(reads.consensus[,"correct"]==FALSE)
precision <- cor.genera/Genera
sensitivity

summary(reads.consensus)
CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")
reads.consensus[,"taxonRank"] <- ordered(reads.consensus[,"taxonRank"], levels = CLASSES)
barplot(table(reads.consensus[,"correct"],droplevels(reads.consensus[,"taxonRank"])))
unique(reads.consensus[,5])
library(plyr)
count(reads.consensus, 'taxonID')
system(command = "unipept taxonomy 365349")
