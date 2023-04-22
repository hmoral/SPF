#!/usr/bin/Rscript
library(dplyr)
library(tidyverse)
library(ggplot2)
library(stringr)
library(readr)

indexf <- "/shared/volume/hologenomics/data/gfemer/Ref_genomes/Terpsiphone_corvina_B10K.fasta.fai"

# 1 - get the ids of the largest scaffolds
index <- read.table(indexf, header =F, sep = "\t")[,c(1,2)]
largest_scaffolds <- index %>% arrange(desc(V2)) %>% top_n(5) %>% select(V1)

# 2 - get the the average persite depth per scaffold (from all individuals)
globalposf <- "/groups/hologenomics/gfemer/data/Real_Reference_results/Global_Depth/Global_depthpersite_modern.pos.gz"

n_all <- 19 # Number of individuals
datalist = list()
for(i in 1:length(largest_scaffolds$V1)){
    scaf <- largest_scaffolds$V1[i]
    dat <- read.table(text = system(paste('zcat',globalposf,'| grep',scaf), intern = T),header= F)[,c(1,3)]
    datalist[[i]] <- dat
}
allscafs <- do.call(rbind, datalist)
normfactor <- mean(allscafs$V3/n_all) # Mean coverage over the 5 largest scaffolds (divided by N because V3 is total counts not average)

# 3 - get the average depth per scaffold per group
get_average <- function(depth,n,normfactor) {
    normdepth <- tryCatch(
        mean( depth$totDepth / n ) / normfactor,
        error = function(err) {
            print(paste0("No depth in ",scaffold," assuming 0"))
            normdepth <- 0
            return(normdepth)
        }
    )
    print(normdepth)
    return(normdepth)
}

min_length <- 500
results_dir <- "/groups/hologenomics/gfemer/data/Real_Reference_results/Depth_per_scaffold"
names <- vector()
depthlist <- list()
# for (dir in list.dirs(results_dir, recursive =F)) {
for (dir in c("Males","Females")) {
    depth_vector <- vector()
    dir <- paste0(results_dir,"/",dir)
    names <- c(names,basename(dir)) 
    cmd <- paste0("wc -l ",paste0(dir,"/",list.files(path=dir,pattern="*.depthSample",recursive = F)[1]))
    n <- as.numeric(unlist(strsplit(system(cmd, intern = T), " "))[1])
    for (scaf in list.files(path=dir,pattern="*.pos.gz",recursive = F)) {
      scaffold <- strsplit(basename(path=scaf),split="\\.")[[1]][1]
      print(scaffold)
      if(index$V2[index$V1 == scaffold] < min_length){
          print(paste0("skipped ",scaffold))
          next
      }
      scaft <- read_tsv(paste0(dir,"/",scaf),col_names=T)[,3]
      scaffold <- strsplit(scaffold,"_")[[1]][2]
      normdepth <- get_average(scaft,n,normfactor)
      depth_vector <- c(depth_vector, normdepth)
      names(depth_vector)[length(depth_vector)] <- scaffold
   }
   depthlist[[length(depthlist)+1]] <- depth_vector 
}
depth <- do.call(cbind, depthlist)
colnames(depth) <- names

write.table(depth,sep="\t", quote = F, row.names =F , file=paste0("Normalized_Depth_all_scaffolds_fl",min_length,".txt"))

# Calculating the differences of coverage between males and females
png(paste0("depth_persex_fl",min_length,".global.r2.png"),width=3000)
depth %>% as.data.frame(row.names = rownames(depth)) %>% 
    rownames_to_column('scaffold') %>% gather("Sex","Normcov", -scaffold) %>% 
    ggplot(.,aes(scaffold,Normcov, fill=Sex)) + geom_bar(stat="identity", position=position_dodge())
dev.off()

png(paste0("depth_persex_fl",min_length,".global.altos.png"),width=1000)
depth %>% as.data.frame(row.names = rownames(depth)) %>% 
    rownames_to_column('scaffold') %>% gather("Sex","Normcov", -scaffold) %>% filter(Normcov > 10) %>%
    ggplot(.,aes(scaffold,Normcov, fill=Sex,group=Sex)) + geom_bar(stat="identity", position=position_dodge())
dev.off()

differences <- depth %>% as.data.frame(row.names = rownames(depth)) %>% 
    rownames_to_column('scaffold') %>% mutate(difference = Males - Females)
write.table(differences, sep="\t", quote = F, row.names = F, file="Depth_differences_MF_500.txt")

# Selecting only scaffolds that have a difference of 0.4
png(paste0("depth_persex_fl",min_length,".diff0.4.png"),width = 1000)
differences %>% filter(difference >= 0.4) %>% select(-difference) %>% gather("Sex","Normcov", -scaffold) %>% 
    ggplot(.,aes(scaffold,Normcov, fill=Sex,group=Sex)) + geom_bar(stat="identity", position=position_dodge()) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 1))
dev.off()

png(paste0("depth_persex_fl",min_length,".diff0.4.nooutliers.png"),width = 1000)
differences %>% filter(abs(difference) >= 0.4) %>% filter(Females <= 1.5) %>% filter(Males <= 1.5) %>% select(-difference) %>% gather("Sex","Normcov", -scaffold) %>% 
    ggplot(.,aes(scaffold,Normcov, fill=Sex,group=Sex)) + 
    geom_bar(stat="identity", position=position_dodge()) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 1))
dev.off()

png(paste0("depth_persex_fl",min_length,".diff0.4.outliers_putative_repetitive.png"),width = 1000)
differences %>% filter(abs(difference) >= 0.4) %>% filter(Females > 1.5) %>% filter(Males > 1.5) %>% select(-difference) %>% gather("Sex","Normcov", -scaffold) %>% 
    ggplot(.,aes(scaffold,Normcov, fill=Sex,group=Sex)) + 
    geom_bar(stat="identity", position=position_dodge()) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 1))
dev.off()

png(paste0("sexchrs_by_difference",min_length,".diff0.4.nooutliers.png"),width = 1000)
differences %>% filter(abs(difference) >= 0.4) %>% filter(Females <= 1.5) %>% filter(Males <= 1.5) %>% select(scaffold, difference) %>% 
mutate(Chromosome = ifelse(difference > 0, "Z","W")) %>%
    ggplot(.,aes(scaffold,difference,group=Chromosome, fill=Chromosome)) + 
    geom_bar(stat="identity", position=position_dodge()) + ylim(-0.7,0.8) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 1))
dev.off()

differences %>% filter(abs(difference) >= 0.4) %>%
write.table(.,sep="\t", quote = F, row.names = F, file=paste0("Sexual_Scaffolds_Depth_differences_MF_",min_length,".txt"))

# Tamaño de chrs sexuales
chrZ <- differences %>% filter(abs(difference) >= 0.4) %>% filter(Females <= 1.5) %>% filter(Males <= 1.5) %>% select(scaffold, difference) %>% mutate(Chromosome = ifelse(difference > 0, "Z","W")) %>% filter(Chromosome == "Z")
chrZ <- chrZ$scaffold
chrZ <- paste0("scaffold_",chrZ)
sum(index$V2[index$V1 %in% chrZ])

chrW <- differences %>% filter(abs(difference) >= 0.4) %>% filter(Females <= 1.5) %>% filter(Males <= 1.5) %>% select(scaffold, difference) %>% mutate(Chromosome = ifelse(difference > 0, "Z","W")) %>% filter(Chromosome == "W")
chrW <- chrW$scaffold
chrW <- paste0("scaffold_",chrW)
sum(index$V2[index$V1 %in% chrW])

Pot_sex <- paste0("scaffold_",differences$scaffold)
write.table(Pot_sex, quote = F, row.names = F, file="Potential_sex_scaffolds.txt")
Regions <- index$V1[!index$V1 %in Pot_sex]
write.table(Regions, quote =F, row.names =F , file="Regions_noSexChr.txt" )

# cut -f1 Sexual_Scaffolds_Depth_differences_MF_500.txt | grep -v "scaffold" | awk '{print "scaffold_"$1}' > Potential_sex_scaffolds.txt
# cut -f1 /shared/volume/hologenomics/data/gfemer/Ref_genomes/Terpsiphone_corvina_B10K.fasta.fai | grep -vwf Potential_sex_scaffolds.txt > Regions_noSexChr.txt