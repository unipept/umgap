tabel <- read.csv("acinetobacter_reads_analyseer.csv")
juist <- tabel[tabel[,7]==1,]
fout <- tabel[tabel[,7]==0,]
juist <- juist[,-1]
fout <- fout[,-1]
juist <- juist[-which(is.na(juist$taxon_rank)),]
fout <- fout[-which(is.na(fout$taxon_rank)),]
juist.n <- as.numeric(paste(juist[,5]))
fout.n <- as.numeric(paste(fout[,5]))
sum(fout.n >=16)/length(fout.n)
sum(juist.n==5)/length(juist.n)
summary(as.numeric(paste(juist[,5])))
summary(as.numeric(paste(fout[,5])))
hist(as.numeric(paste(juist[,5])))
histo <-hist(as.numeric(paste(fout[,5])))
histo

sum(tabel$taxon_id==1)/nrow(tabel)

library(plyr)
CLASSES = c("no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum",
            "subphylum", "superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", 
            "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus",
            "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma")

juist$taxon_rank <- factor(juist$taxon_rank, levels = CLASSES)
juist <- juist[order(juist$taxon_rank),]
counts <- count(juist$taxon_rank)
counts.m <- as.matrix(counts)
as.numeric(counts.m[,2])
counts.m
barplot(as.numeric(counts.m[,2]),names.arg = counts.m[,1])

fout <- fout[-which(fout$taxon_id==1),]

fout$taxon_rank <- factor(fout$taxon_rank, levels = CLASSES)
fout <- fout[order(fout$taxon_rank),]
counts.f <- count(fout$taxon_rank)
counts.m <- as.matrix(counts.f)
as.numeric(counts.m[,2])
counts.m
barplot(as.numeric(counts.m[,2]),names.arg = counts.m[,1])

sum(fout$taxon_rank=="species")
sum(juist$taxon_rank=="species")

#A$animal <- factor(A$animal, levels = c("dog", "elephant","cat"))