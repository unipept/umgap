
taxonomy <- read.csv("Kraken_hiseq_taxonomy.csv")

reads.consensus <- read.csv("Kraken_HiSeq_consensus.csv", header=FALSE)
true.lineages <- read.csv("Kraken_HiSeq_Lineages",fill = TRUE, header=FALSE)
reads.consensus[,ncol(reads.consensus)+1] <- FALSE
for(i in 1:nrow(reads.consensus)){
  if(any(as.character(reads.consensus[i,3])==true.lineages)||as.character(reads.consensus[i,3])=="root"){
    reads.consensus[i,6] <- TRUE
  }else{
    if(reads.consensus[i,5]==0){
      reads.consensus[i,6] <- TRUE
    }else{
      reads.consensus[i,6] <- FALSE
    }
  }
}

colnames(reads.consensus) <- c("Header","taxonID","taxonName","taxonRank","distance","correct")

CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")
reads.consensus[,"taxonRank"] <- ordered(reads.consensus[,"taxonRank"], levels = CLASSES)

sum(reads.consensus[,6])

Voorspellingen <- nrow(reads.consensus)
Genera <- sum(reads.consensus[,"taxonRank"]>="genus")
cor.genera <- sum(reads.consensus[reads.consensus[,"taxonRank"]>="genus","correct"]==TRUE)
incor.voorspelling <- sum(reads.consensus[,"correct"]==FALSE)
precision <- cor.genera/Genera
sensitivity <- cor.genera/Voorspellingen


barplot(table(reads.consensus[,"correct"],droplevels(reads.consensus[,"taxonRank"])))

reads.consensus[reads.consensus[,"taxonRank"]>"species",]

reads.consensus[reads.consensus[,"taxonRank"]=="species",]
Foute_reads <- reads.consensus[reads.consensus[,"correct"]==FALSE,]
barplot(table(droplevels(Foute_reads[,"taxonRank"])))

summary(Foute_reads[,"distance"])
barplot(table(Foute_reads[,"distance"]))

species <- reads.consensus[reads.consensus[,"taxonRank"]=="species",]
species[species[,"correct"]==FALSE,]
