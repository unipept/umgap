foundTaxons <- read.csv("setA1_1.1.taxonListWithLineage",header = FALSE)
metabench_genus <- read.csv("metabench_genus.csv",sep=";")
present_genera <- metabench_genus[,1]
metabench_phylum <- read.csv("metabench_phylum.csv",sep=";")
present_phyla <- metabench_phylum[2:18,1]

CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")

header <- c("amount","taxon_ID","classification","depth",CLASSES[-1])
colnames(foundTaxons) <- header

summary(foundTaxons[,"phylum"])
tot_identification <- sum(foundTaxons[,1])
total <- 28912773

tot_identification/total
