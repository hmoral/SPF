library(data.table)
library(tidyverse)
library(patchwork)

# sample metadata 
popInfo=fread("1_SPF_popInfo.txt")
# time    indv    sex       age   pop                          loc   year
# 1:     Modern   SPF02   male sub-adult laD_M Behind MarieLiz grandmothers 2007-8
# 2:     Modern   SPF03   male     adult laD_M      LUnion old peoples home 2007-8
# 3:     Modern   SPF04 female     adult laD_M                Allee Kersley 2007-8
# 4: Historical SPF1086   <NA>      <NA> laD_H                   La Digue H   1880
# 5: Historical  SPF261   <NA>      <NA> laD_H                   La Digue H   1880
# 6: Historical  SPF262   <NA>      <NA> laD_H                   La Digue H   1880
# 7: Historical   SPF2H   <NA>      <NA> laD_H                   La Digue H   1888

# snpeff results
master_refANC_eff=fread("1_SPF_snpEff_output.txt.gz")
# CHROM    POS REF ALT  DP                                                                                                                      GEN.GT EFF.IMPACT
# 1:     scaffold_1   8562   G   T 194 0/0,0/1,0/0,0/0,0/0,1/1,./.,./.,0/0,0/0,0/0,0/0,0/1,0/0,0/0,0/0,./.,0/1,0/0,0/1,1/1,1/1,0/0,./.,./.,./.,0/1,./.,1/1,1/1,./.        LOW
# 2:     scaffold_1  11717   C   A 221 ./.,0/0,0/0,0/0,0/0,0/0,0/0,./.,0/0,0/0,0/0,0/0,0/0,0/0,0/0,0/0,0/0,./.,0/0,./.,./.,0/0,0/0,0/0,0/1,0/0,0/0,1/1,0/0,0/1,./.        LOW
# 3:     scaffold_1  18629   G   C 169 0/0,0/1,0/0,0/0,0/0,1/1,0/0,./.,0/0,0/0,0/0,./.,1/1,./.,./.,0/0,0/0,0/1,0/0,0/0,0/0,0/1,./.,./.,0/0,./.,0/0,0/0,0/1,./.,0/0   MODERATE
# 4:     scaffold_1  18831   T   G 166 0/0,0/1,0/0,0/0,0/0,1/1,0/0,1/1,0/0,0/0,0/0,./.,1/1,0/0,./.,0/0,0/0,0/1,0/0,./.,./.,1/1,./.,./.,./.,./.,./.,0/0,1/1,./.,./.   MODERATE
# 5:     scaffold_1 394608   C   G 244 0/0,0/0,0/1,0/0,0/1,0/0,0/1,./.,0/0,0/0,0/0,0/0,0/0,0/0,0/0,0/1,./.,0/0,0/0,./.,./.,0/0,0/0,./.,./.,0/0,0/0,1/1,0/1,./.,./.        LOW

boot_alleles=data.table()
boot_rxy=data.table()
for(BOOT in 1:50){
  print("ITERATION")
  print(BOOT)
  #
  refANC_eff=master_refANC_eff
  ### metadata
  refANC_indv=fread("tmp_snpeff_results_vcfMergeBCFtools_refSPF_individuals.txt",header = F)
  names(refANC_indv)="indv"
  # list of bam files that gives the order of individuals for the genotype calls
  # 1: bams/SPF02.merged.dedup.realigned.bam
  # 2: bams/SPF03.merged.dedup.realigned.bam
  # 3: bams/SPF04.merged.dedup.realigned.bam
  # 4: bams/SPF05.merged.dedup.realigned.bam
  # 5: bams/SPF06.merged.dedup.realigned.bam
  # 6: bams/SPF07.merged.dedup.realigned.bam
  # remove the paths from the names to retain the ID only. Add an order identifier and merge it with the sample metadata
  refANC_indv$indv=gsub("/groups/hologenomics/gfemer/data/Data/Historical/Bamfiles/Realigned/|bams/|.merged.dedup.realigned.bam","",refANC_indv$indv)
  refANC_indv$order=1:nrow(refANC_indv)
  refANC_indv=merge(refANC_indv,popInfo,by = "indv")
  refANC_indv=refANC_indv[indv!="SPF02"] # remove the highly related individual

  # equalize modern and historical samples size
  GENO=data.table(str_split_fixed(refANC_eff$GEN.GT,",",Inf))
  names(GENO)=paste(refANC_indv$indv,refANC_indv$time,sep="_")
  refANC_indv=refANC_indv %>% 
    group_by(time) %>%
    slice_sample(n = 13) %>% data.table()
  refANC_indv=refANC_indv[order(order)]
  GENO=GENO[,paste(refANC_indv$indv,refANC_indv$time,sep="_"),with=F]
  GENO=apply(GENO, 1, paste, collapse=",")
  refANC_eff$GEN.GT=GENO

  # ANC alleles - consensus alleles from the sister species assemblies
  refANC_alleles=fread("summary_out_counts_angsd_Tcinnamomea_dumpCounts4.merged.gz",header = T)
  refANC_alleles=refANC_alleles[,names(refANC_alleles)[c(1,2,8,13,14)],with=F]
  names(refANC_alleles)=c("CHROM","POS","anc_alleles","anc_all1","anc_all2")
  # CHROM     POS anc_alleles anc_all1 anc_all2
  #   1:  scaffold_1    3743           1        T     <NA>
  #   2:  scaffold_1    3821           1        C     <NA>
  #   3:  scaffold_1    4590           1        A     <NA>
  #   4:  scaffold_1    4978           1        C     <NA>
  #   5:  scaffold_1    5768           1        A     <NA>
  # ---                                                  
  #   57420: scaffold_83 2026404           2        T        A
  # 57421: scaffold_83 2026429           2        A        T
  # 57422: scaffold_83 2028494           2        C        T
  # 57423: scaffold_83 2028760           2        G        C
  # 57424: scaffold_83 2028774           2        C        T
  POS_ANC=paste(refANC_alleles$CHROM,refANC_alleles$POS,sep="_")
  POS_REF=paste(refANC_eff$CHROM,refANC_eff$POS,sep="_")
  sum(POS_REF%in%POS_ANC)/length(POS_REF) # 92% of called alleles were found in the ancestral set
  sum(!POS_REF%in%POS_ANC) # 66982 were not
  table(refANC_eff[!POS_REF%in%POS_ANC]$EFF.IMPACT) # they occur across all categories
  hist(refANC_eff[!POS_REF%in%POS_ANC]$DP) # they have all kinds of depth

  # make final ancestral ref by selecting randomly one allele
  oneAllele=refANC_alleles[anc_alleles==1]
  nrow(oneAllele) # [1] 737734 - 93%
  oneAllele$anc_ALLELE=oneAllele$anc_all1
  twoAllele=refANC_alleles[anc_alleles==2]
  nrow(twoAllele) # [1] 57424 - 7%
  twoAllele=data.frame(twoAllele)
  twoAllele=twoAllele %>% 
    mutate(pmap_dfr(across(anc_all1:anc_all2), ~ `[<-`(c(...), !seq_along(c(...)) %in% sample(which(!is.na(c(...))), 1), NA))) %>% data.table()
  twoAllele=melt(twoAllele,id.vars=names(twoAllele)[1:3])
  twoAllele=twoAllele[!is.na(value)]
  twoAllele$anc_ALLELE=twoAllele$value
  refANC_alleles=rbind(oneAllele,twoAllele,fill=T)
  refANC_alleles=refANC_alleles[order(CHROM,POS)]
  
  # merge called alleles with ancestral alles, rataining ALL snp eff snps
  refANC_eff=merge(refANC_eff,refANC_alleles[,c("CHROM","POS","anc_ALLELE"),with=F],by=c("CHROM","POS"),all.x=T)

  ####
  # Assume ancestral for all matches
  ####
  # retain sites with ancestral info and without
  refANC_eff_WITHinfo=refANC_eff[!is.na(anc_ALLELE)] # 88%
  # split those that the reference match the ancestral and does that do not match
  refANC_match_REF=refANC_eff_WITHinfo[refANC_eff_WITHinfo$anc_ALLELE==refANC_eff_WITHinfo$REF]
  # find the ones in which the match to the ancestral is the alternative and swap the order
  refANC_match_ALT=refANC_eff_WITHinfo[refANC_eff_WITHinfo$anc_ALLELE==refANC_eff_WITHinfo$ALT]
  GENO=data.table(str_split_fixed(refANC_match_ALT$GEN.GT,",",Inf))
  GENO[GENO=="0/0"]="ALT/ALT"
  GENO[GENO=="0/1"]="ALT/REF"
  GENO[GENO=="1/0"]="REF/ALT"
  GENO[GENO=="1/1"]="REF/REF"
  GENO[GENO=="ALT/ALT"]="1/1"
  GENO[GENO=="ALT/REF"]="1/0"
  GENO[GENO=="REF/ALT"]="0/1"
  GENO[GENO=="REF/REF"]="0/0"
  GENO=apply(GENO, 1, paste, collapse=",")
  refANC_match_ALT$GEN.GT=GENO
  ### 
  refANC_match=rbind(refANC_match_REF,refANC_match_ALT)
  refANC_eff=refANC_match
  
  ####
  # add per site stats
  ####
  ##
  GENO=data.table(str_split_fixed(refANC_eff$GEN.GT,",",Inf))
  names(GENO)=paste(refANC_indv$indv,refANC_indv$time,sep="_")
  #
  all_miss=rowSums(GENO=="./.")*2
  all_ref=rowSums(GENO=="0/0")*2
  all_het=rowSums(GENO=="0/1" | GENO=="1/0")
  all_alt=rowSums(GENO=="1/1")*2
  all_all=ncol(GENO)*2
  (all_all) == (all_miss[1]+all_ref[1]+all_het[1]*2+all_alt[1])
  all_miss_p=all_miss/all_all
  all_freq_ref=((all_ref)+(all_het))/(all_all-all_miss)
  all_freq_alt=((all_alt)+(all_het))/(all_all-all_miss)
  
  mod_miss=rowSums(GENO[,grep("Modern",names(GENO)),with=F]=="./.")*2
  mod_ref=rowSums(GENO[,grep("Modern",names(GENO)),with=F]=="0/0")*2
  mod_het=rowSums(GENO[,grep("Modern",names(GENO)),with=F]=="0/1" | GENO[,grep("Modern",names(GENO)),with=F]=="1/0")
  mod_alt=rowSums(GENO[,grep("Modern",names(GENO)),with=F]=="1/1")*2
  mod_all=ncol(GENO[,grep("Modern",names(GENO)),with=F])*2
  mod_all == (mod_miss[1]+mod_ref[1]+mod_het[1]*2+mod_alt[1])
  mod_miss_p=mod_miss/mod_all
  mod_freq_ref=((mod_ref)+(mod_het))/(mod_all-mod_miss)
  mod_freq_alt=((mod_alt)+(mod_het))/(mod_all-mod_miss)
  hist(mod_freq_alt)
  
  hist_miss=rowSums(GENO[,grep("Historical",names(GENO)),with=F]=="./.")*2
  hist_ref=rowSums(GENO[,grep("Historical",names(GENO)),with=F]=="0/0")*2
  hist_het=rowSums(GENO[,grep("Historical",names(GENO)),with=F]=="0/1" | GENO[,grep("Historical",names(GENO)),with=F]=="1/0")
  hist_alt=rowSums(GENO[,grep("Historical",names(GENO)),with=F]=="1/1")*2
  hist_all=ncol(GENO[,grep("Historical",names(GENO)),with=F])*2
  hist_all == (hist_miss[1]+hist_ref[1]+hist_het[1]+hist_het[1]+hist_alt[1])
  hist_miss_p=hist_miss/hist_all
  hist_freq_ref=((hist_ref)+(hist_het))/(hist_all-hist_miss)
  hist_freq_alt=((hist_alt)+(hist_het))/(hist_all-hist_miss)
  hist(hist_freq_alt)
  
  all_GENO=data.table(all_all,all_miss,all_ref,all_het,all_alt,all_miss_p,all_freq_ref,all_freq_alt,
                      mod_all,mod_miss,mod_ref,mod_het,mod_alt,mod_miss_p,mod_freq_ref,mod_freq_alt,
                      hist_all,hist_miss,hist_ref,hist_het,hist_alt,hist_miss_p,hist_freq_ref,hist_freq_alt)

  ####
  # merge everything
  ####
  refANC_eff=cbind(refANC_eff,all_GENO)

  ####
  # calculate rxy and count alleles
  ####
  X=refANC_eff
  x_rxy=X%>%group_by(EFF.IMPACT)%>%summarise(freq_hist=sum(hist_freq_alt*(1-mod_freq_alt)),freq_mod=sum(mod_freq_alt*(1-hist_freq_alt)))%>%data.table()
  x_rxy$rxy=x_rxy$freq_hist/x_rxy$freq_mod
  x_rxy$boot=BOOT
  boot_rxy=rbind(boot_rxy,x_rxy)

  GENO=data.table(str_split_fixed(X$GEN.GT,",",Inf))
  names(GENO)=paste(refANC_indv$indv,refANC_indv$time,sep="_")
  GENO=cbind(GENO,X)
  GENO_melt=melt(GENO,id.vars = names(GENO)[-c(grep("_Historical|_Modern",names(GENO)))])
  GENO_melt$variable=as.character(GENO_melt$variable)
  GENO_melt$indv=str_split_fixed(GENO_melt$variable,"_",2)[,1]
  GENO_melt$time=str_split_fixed(GENO_melt$variable,"_",2)[,2]
  GENO_melt=GENO_melt%>%group_by(indv,time)%>%mutate(totsyn=sum(EFF.IMPACT=="LOW"&value!="./.",na.rm = T)*2)%>%data.table()
  GENO_melt=merge(GENO_melt,popInfo[,c("indv","pop"),with=],by="indv")
  GENO_melt_indv=GENO_melt %>% 
    group_by(indv,pop,time,EFF.IMPACT,totsyn) %>% 
    summarise(totAlelles=n()*2,
              NAs=sum(value=="./.",na.rm = T)*2,
              hom_ref=sum(value=="0/0",na.rm = T)*2/totsyn,
              hom_alt=sum(value=="1/1",na.rm = T)*2/totsyn,
              het_der=sum(value=="0/1"|value=="1/0",na.rm = T)/totsyn,
              tot_der=((sum(value=="1/1",na.rm = T)*2)+(sum(value=="0/1"|value=="1/0",na.rm = T)))/totsyn) %>% data.table()
  GENO_melt_indv=GENO_melt_indv[!duplicated(paste(indv,pop,time,EFF.IMPACT))]
  x_alleles=GENO_melt_indv
  x_alleles$boot=BOOT
  x_alleles$N_SNPS=nrow(X)
  boot_alleles=rbind(boot_alleles,x_alleles)
}

