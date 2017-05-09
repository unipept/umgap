setwd("/Users/Aranka/Documents/Unief/Thesis/Thesis/Kmer2consensus/OptiesVGL")

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
  filename <- paste0("Kraken_HiSeq.11.Opt",k,".csv")
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
  #sensitivity <- c(correct.rank/tot.reads,mean(correct.rank/tot.reads))
  sensitivity <- c(correct.rank/(correct.rank+(tot.reads-nr.rank)),mean(correct.rank/(correct.rank+(tot.reads-nr.rank))))
  all.sensitivity <- rbind(all.sensitivity,sensitivity)
  all.specificity <- rbind(all.specificity,specificity)
  print(k)
}

tot_perSP <- 1000
IDs_perSP <- summary(reads[,"correct"])
perc_IDs_perSP <- IDs_perSP/tot_perSP

#### Zowel genus als species fout
beide_fout <- reads[intersect(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)),which(correctness[,"genus_name"]==F)),c("correct","species_name","genus_name")]
summary(beide_fout[,"correct"])/IDs_perSP
summary(beide_fout[beide_fout[,"correct"]=="A_hydrophila",])
summary(beide_fout[beide_fout[,"correct"]=="B_cereus",])
summary(beide_fout[beide_fout[,"correct"]=="B_fragilis",])
summary(beide_fout[beide_fout[,"correct"]=="M_abscessus",])
summary(beide_fout[beide_fout[,"correct"]=="P_fermentans",])
summary(beide_fout[beide_fout[,"correct"]=="R_sphaeriodes",])
summary(beide_fout[beide_fout[,"correct"]=="S_aureus",])
summary(beide_fout[beide_fout[,"correct"]=="S_pneumoniae",])
summary(beide_fout[beide_fout[,"correct"]=="V_cholerae",])
summary(beide_fout[beide_fout[,"correct"]=="X_axanopodis",])
#### Genus juist, species fout
genus_juist <- reads[intersect(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)),which(correctness[,"genus_name"])),c("correct","species_name","genus_name")]
summary(genus_juist[,"correct"])/IDs_perSP
summary(genus_juist[genus_juist[,"correct"]=="A_hydrophila",])
summary(genus_juist[genus_juist[,"correct"]=="B_cereus",])
summary(genus_juist[genus_juist[,"correct"]=="B_fragilis",])
summary(genus_juist[genus_juist[,"correct"]=="M_abscessus",])
summary(genus_juist[genus_juist[,"correct"]=="P_fermentans",])
summary(genus_juist[genus_juist[,"correct"]=="R_sphaeroides",])
summary(genus_juist[genus_juist[,"correct"]=="S_aureus",])
summary(genus_juist[genus_juist[,"correct"]=="S_pneumoniae",])
summary(genus_juist[genus_juist[,"correct"]=="V_cholerae",])
summary(genus_juist[genus_juist[,"correct"]=="X_axanopodis",])

#### Percentage van foute species met correct genus
length(intersect(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)),which(correctness[,"genus_name"])))/length(intersect(which(reads[,"taxon_rank"]>="species"),which(correctness[,"species_name"]==FALSE)))

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
correct_mat <- rbind(corr_kingdom,corr_phylum,corr_class,corr_order,corr_family,corr_genus,corr_species)
colnames(correct_mat) <- levels(reads[,"correct"])
correct_dataframe  <- as.data.frame(correct_mat)

plot(correct_dataframe$A_hydrophila,type='n',ylim=c(0,1000))
for(i in 1:10){
  lines(correct_dataframe[i], col=rainbow(10)[i])
}
legend("bottomleft",legend=levels(reads[,"correct"]), col=rainbow(10),fill=rainbow(10),ncol=2,cex=0.5)
correct_dataframe[1]

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

write.csv(mat,file="HiSeq_PR_SE.csv",sep=",")

data_frame <- as.data.frame(mat)
summary(data_frame)

plot(mat[,"superkingdom_sens"],mat[,"superkingdom_prec"])

library(car) 
?scatterplot
scatterplot(mean_prec ~ mean_sens | seed_length,data=data_frame, smoother=FALSE,reg.line=FALSE,
            labels = data_frame[,"type"],id.n = nrow(data_frame),id.cex = 0.7)
scatterplot(mean_prec ~ mean_sens | seed_length,data=data_frame, 
            smoother=FALSE,reg.line=FALSE,labels = data_frame[,"type"],id.n = nrow(data_frame),
            id.cex = 0.7,xlim=c(0.66,0.76),ylim=c(0.95,0.99))
scatterplot(mean_prec ~ mean_sens | type,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | gap_size,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | gap_pen,data=data_frame, smoother=FALSE,reg.line=FALSE)

scatterplot(genus_prec ~ genus_sens | seed_length,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(genus_prec ~ genus_sens | type,data=data_frame, smoother=FALSE,reg.line=FALSE)

scatterplot(phylum_prec ~ phylum_sens | seed_length,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(genus_prec~genus_sens| seed_length,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(species_prec~species_sens| seed_length,data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(phylum_prec ~ phylum_sens | type,data=data_frame, smoother=FALSE,reg.line=FALSE)

scatterplot(mean_prec ~IDs | type, data=data_frame, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_sens ~IDs | type, data=data_frame, smoother=FALSE,reg.line=FALSE)

scatterplot(IDs~seed_length, data = data_frame, smoother=FALSE)
scatterplot(IDs~type|seed_length,data = data_frame, smoother=FALSE)

scatterplot(seed_length~mean_prec | type, data=data_frame, smoother=FALSE)
scatterplot(seed_length~mean_sens | type, data=data_frame, smoother=FALSE)

scatterplot(type~mean_prec | seed_length, data=data_frame, smoother=FALSE)
scatterplot(type~mean_sens | seed_length, data=data_frame, smoother=FALSE)

scatterplot(seed_length~phylum_prec | type, data=data_frame, smoother=FALSE)
scatterplot(seed_length~phylum_sens | type, data=data_frame, smoother=FALSE)

scatterplot(type~phylum_prec | seed_length, data=data_frame, smoother=FALSE)
scatterplot(type~phylum_sens | seed_length, data=data_frame, smoother=FALSE)

data_frame$length_type <- paste(data_frame$seed_length,data_frame$type, sep=";")

scatterplot(mean_prec~gap_size|length_type,data=data_frame,smoother=FALSE)
scatterplot(mean_prec~gap_size|type,data=data_frame,smoother=FALSE)
scatterplot(mean_sens~gap_size|length_type,data=data_frame,smoother=FALSE)
scatterplot(mean_sens~gap_size|type,data=data_frame,smoother=FALSE)

scatterplot(mean_prec ~ mean_sens | length_type,data=data_frame, smoother=FALSE,reg.line=FALSE)

summary(lm(mean_sens~length_type + gap_pen,data=data_frame))
summary(lm(mean_prec~length_type + gap_pen,data=data_frame))


library(MASS)
parcoord(data_frame[,6:12], col=rainbow(51), var.label=TRUE)



##############################
## SUBSET ####################
##############################

gap2 <- which(mat[,"gap_size"]==2)
gap_pen28 <- c(which(mat[,"gap_pen"]==0.2),which(mat[,"gap_pen"]==0.8))
mat_reduced <- mat[-gap_pen28,]
data_frame_reduced <- as.data.frame(mat_reduced)
data_frame_reduced$length_type <- paste(data_frame_reduced$seed_length,data_frame_reduced$type, sep=";")

scatterplot(mean_spec ~ mean_sens | seed_length,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE,labels = rownames(mat_reduced),id.n = nrow(data_frame_reduced),
            id.cex = 0.7)
scatterplot(mean_spec ~ mean_sens | seed_length,data=data_frame_reduced, 
            smoother=FALSE,reg.line=FALSE,labels = data_frame_reduced[,"type"],id.n = nrow(data_frame_reduced),
            id.cex = 0.7,xlim=c(0.66,0.76),ylim=c(0.95,0.99))
scatterplot(mean_spec ~ mean_sens | type,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(seed_length ~ IDs | type, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(type ~ IDs | seed_length, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(genus_prec ~ genus_sens | seed_length,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE,labels = rownames(mat_reduced),id.n = nrow(data_frame_reduced),
            id.cex = 0.7,xlim=c(0.66,0.76),ylim=c(0.95,0.99))
scatterplot(species_prec ~ species_sens | seed_length,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE,labels = rownames(mat_reduced),id.n = nrow(data_frame_reduced),
            id.cex = 0.7,xlim=c(0.31,0.36),ylim=c(0.85,0.96))

scatterplot(mean_prec ~IDs | seed_length, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_sens ~IDs | seed_length, data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)

scatterplot(mean_prec ~ mean_sens | length_type,data=data_frame_reduced, smoother=FALSE,reg.line=FALSE)
scatterplot(mean_prec ~ mean_sens | length_type,data=data_frame_reduced,   smoother=FALSE,reg.line=FALSE,xlim=c(0.65,0.76),ylim=c(0.95,0.99))
###########################
## Sorteeroefeningen ######
###########################
# "superkingdom","phylum","class","order","family","genus","species"

more_info <- cbind(mat_reduced, mat_reduced[,"mean_sens"]*mat_reduced[,"mean_prec"],
                   mat_reduced[,"superkingdom_sens"]*mat_reduced[,"superkingdom_prec"],
                   mat_reduced[,"phylum_sens"]*mat_reduced[,"phylum_prec"], 
                   mat_reduced[,"class_sens"]*mat_reduced[,"class_prec"],
                   mat_reduced[,"order_sens"]*mat_reduced[,"order_prec"],
                   mat_reduced[,"family_sens"]*mat_reduced[,"family_prec"],
                   mat_reduced[,"genus_sens"]*mat_reduced[,"genus_prec"],
                   mat_reduced[,"species_sens"]*mat_reduced[,"species_prec"])
rownames(more_info) <- rownames(mat_reduced)
per_mean <- sort(more_info[,23],decreasing = TRUE,index.return = TRUE)$x
names(per_mean) <- rownames(more_info)[sort(more_info[,23],decreasing = TRUE,index.return = TRUE)$ix]
per_kingdom <- sort(more_info[,24],decreasing=TRUE,index.return = TRUE)$x
names(per_kingdom) <- rownames(more_info)[sort(more_info[,24],decreasing=TRUE,index.return = TRUE)$ix]
per_phylum <- sort(more_info[,25],decreasing = TRUE,index.return = TRUE)$x
names(per_phylum) <- rownames(more_info)[sort(more_info[,25],decreasing = TRUE,index.return = TRUE)$ix]
per_class <- sort(more_info[,26],decreasing = TRUE,index.return = TRUE)$x
names(per_class) <- rownames(more_info)[sort(more_info[,26],decreasing = TRUE,index.return = TRUE)$ix]
per_order <- sort(more_info[,27],decreasing = TRUE,index.return = TRUE)$x
names(per_order) <- rownames(more_info)[sort(more_info[,27],decreasing = TRUE,index.return = TRUE)$ix]
per_family <- sort(more_info[,28],decreasing = TRUE,index.return = TRUE)$x
names(per_family) <- rownames(more_info)[sort(more_info[,28],decreasing = TRUE,index.return = TRUE)$ix]
per_genus <- sort(more_info[,29],decreasing = TRUE,index.return = TRUE)$x
names(per_genus) <- rownames(more_info)[sort(more_info[,29],decreasing = TRUE,index.return = TRUE)$ix]
per_species <- sort(more_info[,30],decreasing = TRUE,index.return = TRUE)$x
names(per_species) <- rownames(more_info)[sort(more_info[,30],decreasing = TRUE,index.return = TRUE)$ix]

intersect(intersect(intersect(intersect(intersect(intersect(intersect(names(per_mean[1:15]),
                                                            names(per_kingdom[1:15])),
                                                  names(per_phylum[1:15])),
                                        names(per_class[1:15])),
                              names(per_order[1:15])),
                    names(per_family[1:15])),
          names(per_genus[1:15])),
  names(per_species[1:15]))

sens_kingdom <- sort(mat_reduced[,"superkingdom_sens"],decreasing = TRUE,index.return=TRUE)$x
names(sens_kingdom) <- rownames(mat_reduced)[sort(mat_reduced[,"superkingdom_sens"],decreasing = TRUE,index.return=TRUE)$ix]
sens_phylum <- sort(mat_reduced[,"phylum_sens"],decreasing = TRUE,index.return=TRUE)$x
names(sens_phylum) <- rownames(mat_reduced)[sort(mat_reduced[,"phylum_sens"],decreasing = TRUE,index.return=TRUE)$ix]
sens_class <- sort(mat_reduced[,"class_sens"],decreasing = TRUE,index.return=TRUE)$x
names(sens_class) <- rownames(mat_reduced)[sort(mat_reduced[,"class_sens"],decreasing = TRUE,index.return=TRUE)$ix]
sens_order <- sort(mat_reduced[,"order_sens"],decreasing = TRUE,index.return=TRUE)$x
names(sens_order) <- rownames(mat_reduced)[sort(mat_reduced[,"order_sens"],decreasing = TRUE,index.return=TRUE)$ix]
sens_family <- sort(mat_reduced[,"family_sens"],decreasing = TRUE,index.return=TRUE)$x
names(sens_family) <- rownames(mat_reduced)[sort(mat_reduced[,"family_sens"],decreasing = TRUE,index.return=TRUE)$ix]
sens_genus <- sort(mat_reduced[,"genus_sens"],decreasing = TRUE,index.return=TRUE)$x
names(sens_genus) <- rownames(mat_reduced)[sort(mat_reduced[,"genus_sens"],decreasing = TRUE,index.return=TRUE)$ix]
sens_species <- sort(mat_reduced[,"species_sens"],decreasing = TRUE,index.return=TRUE)$x
names(sens_species) <- rownames(mat_reduced)[sort(mat_reduced[,"species_sens"],decreasing = TRUE,index.return=TRUE)$ix]

intersect(intersect(intersect(intersect(intersect(intersect(names(sens_kingdom[1:15]),
                                                            names(sens_phylum[1:15])),
                                                  names(sens_class[1:15])),
                                        names(sens_order[1:15])),
                              names(sens_family[1:15])),
                    names(sens_genus[1:15])),
          names(sens_species[1:15]))

prec_kingdom <- sort(mat_reduced[,"superkingdom_prec"],decreasing = TRUE,index.return=TRUE)$x
names(prec_kingdom) <- rownames(mat_reduced)[sort(mat_reduced[,"superkingdom_prec"],decreasing = TRUE,index.return=TRUE)$ix]
prec_phylum <- sort(mat_reduced[,"phylum_prec"],decreasing = TRUE,index.return=TRUE)$x
names(prec_phylum) <- rownames(mat_reduced)[sort(mat_reduced[,"phylum_prec"],decreasing = TRUE,index.return=TRUE)$ix]
prec_class <- sort(mat_reduced[,"class_prec"],decreasing = TRUE,index.return=TRUE)$x
names(prec_class) <- rownames(mat_reduced)[sort(mat_reduced[,"class_prec"],decreasing = TRUE,index.return=TRUE)$ix]
prec_order <- sort(mat_reduced[,"order_prec"],decreasing = TRUE,index.return=TRUE)$x
names(prec_order) <- rownames(mat_reduced)[sort(mat_reduced[,"order_prec"],decreasing = TRUE,index.return=TRUE)$ix]
prec_family <- sort(mat_reduced[,"family_prec"],decreasing = TRUE,index.return=TRUE)$x
names(prec_family) <- rownames(mat_reduced)[sort(mat_reduced[,"family_prec"],decreasing = TRUE,index.return=TRUE)$ix]
prec_genus <- sort(mat_reduced[,"genus_prec"],decreasing = TRUE,index.return=TRUE)$x
names(prec_genus) <- rownames(mat_reduced)[sort(mat_reduced[,"genus_prec"],decreasing = TRUE,index.return=TRUE)$ix]
prec_species <- sort(mat_reduced[,"species_prec"],decreasing = TRUE,index.return=TRUE)$x
names(prec_species) <- rownames(mat_reduced)[sort(mat_reduced[,"species_prec"],decreasing = TRUE,index.return=TRUE)$ix]

intersect(intersect(intersect(intersect(intersect(intersect(names(prec_kingdom[1:15]),
                                                            names(prec_phylum[1:15])),
                                                  names(prec_class[1:15])),
                                        names(prec_order[1:15])),
                              names(prec_family[1:15])),
                    names(prec_genus[1:15])),
          names(prec_species[1:15]))
