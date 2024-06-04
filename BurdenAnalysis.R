## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
require(dplyr)
require(survminer)
require(stringr)
require(survival)
require(readxl)
require(readr)
require(ggplot2)
require(cluster)
require(ggdendro)
require(MesKit)
require(ggcorrplot)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
is_zero <- function(n){
  val=ifelse(n==0, FALSE, TRUE)
  return(val)
}

inject.dots <- function(df) {names(df) <- sub(" ", ".", names(df));df}


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
Caris = read_xlsx('CARIS Data_10-27-23 with Clinical Info_deidentified.xlsx', sheet = 1)
clinical = read_xlsx('CARIS Data_10-27-23 with Clinical Info_deidentified.xlsx', sheet = 2)

names(Caris)=make.names(names(Caris))
names(clinical)=make.names(names(clinical))

# Drop and rename duplicate columns
Caris <- Caris %>% select(-c(`MN1.Fusion...341`,`TP53.Mutation...564`)) %>% rename(`TP53.Mutation`=`TP53.Mutation...565`) %>% rename(`MN1.Fusion`=`MN1.Fusion...342`)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
clinical$Seizure.presentation.ANY = clinical$Seizure.presentation..0.No..1.Yes. | clinical$Seizure.presentation.ONLY.after.surgery


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
id_cols = c("Deidentified.code","TMB","MSI","LOH","HLA.A","HLA.B","HLA.C")

allcounts <- colSums((Caris != 0))
varcounts_genes <- colSums((Caris %>% select(-c(id_cols)) != 0))
unique_values_per_column <- lapply(Caris, unique)
varcounts_pp <- apply(Caris %>% select(-id_cols), 1, function(row) sum(row != 0))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
Caris_red <- Caris[allcounts>0]

#Recode LOH factors
Caris_red$LOH=ifelse(Caris_red$LOH %in% c('Quality Not Sufficient','QNS','NA'), 
                     NA, Caris_red$LOH)
Caris_red$MSI=ifelse(Caris_red$MSI %in% c('Quality Not Sufficient','QNS','NA'), 
                     NA, Caris_red$MSI)
# Recode TMB
Caris_red$TMB_numeric = as.numeric(substr(Caris_red$TMB,1,2))

# Specific columns and counts
cn_columns <- grep("CN", names(Caris_red), value = TRUE)
varcounts_cn <-  colSums(Caris %>% select(contains("CNA")) != 0)
  
fusion_columns <- grep("Fusion", names(Caris_red), value = TRUE)
varcounts_fusion <-  colSums(Caris %>% select(contains("Fusion")) != 0)
fusion_pp <- rowSums(Caris %>% select(contains("Fusion")) != 0)

mut_columns <-  grep("Mutation", names(Caris_red), value = TRUE)
varcounts_mut <-  colSums(Caris %>% select(contains("Mutation")) != 0)
mutation_pp <- rowSums(Caris %>% select(contains("Mutation")) != 0)

egfr_columns <-  grep("EGFR", names(Caris_red), value = TRUE)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Add row sums 
pp_df <- data.frame(Deidentified.code=Caris_red$Deidentified.code, fusions = rowSums(Caris_red %>% select(contains("Fusion")) != 0), CNAs = rowSums(Caris_red %>% select(contains("CNA")) != 0), muts = rowSums(Caris_red %>% select(contains("Mutation")) != 0))

# Add gene information
pp_df<- left_join(pp_df, Caris_red)

# Add clinical
pp_df <- left_join(pp_df, clinical %>% select(Deidentified.code, age, sex, kps_at_diagnosis,Seizure.presentation..0.No..1.Yes., Seizure.presentation.ONLY.after.surgery, Seizure.presentation.ANY,Seizure.freq.AVG.over.1st.6.mo, Keppra, Vimpat, Zonegran, Lamictal, race,ki_67, idh1_mutated,mgmt_methylated, p53_percent_reactivity, loc_other, loc_frontal, loc_parietal, loc_occipital, loc_temporal))

# Full
Caris_clinical <- left_join(Caris_red %>% na.omit(), clinical) 


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
length(varcounts_fusion!=0)
length(varcounts_mut!=0)
length(varcounts_cn!=0)
sum(varcounts_fusion!=0)
sum(varcounts_mut!=0)
sum(varcounts_cn!=0)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
cn_top = names((sort(varcounts_cn, decreasing=TRUE)[1:5]))
mut_top = sort(varcounts_mut, decreasing=TRUE)[1:6]
fus_top = sort(varcounts_fusion, decreasing=TRUE)[1:5]


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------

ggplot(pp_df, aes(y=CNAs, x=as.factor(is.na(p53_percent_reactivity)))) +
  geom_boxplot(fill = "skyblue", color = "black", alpha = 0.7) +
  labs(x='', y = "Frequency") + 
  theme_minimal() 
ggsave(file="fusions.png", width=2, height=4, dpi=300)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create a data frame for ggplot
#pp_df 

# Medians 
summary(pp_df$muts)
summary(pp_df$fusions)
summary(pp_df$CNAs)

# Convert categorical seizure to multi
pp_df$Sz_Pre_Diagnosis = pp_df$Seizure.presentation..0.No..1.Yes.
pp_df <- pp_df |> mutate(CodeSeizure=Sz_Pre_Diagnosis*1 + Seizure.presentation.ONLY.after.surgery*2) |> mutate(CodeSeizure=as.factor(CodeSeizure))


# Wide to Long
pp_df_long <- tidyr::pivot_longer(pp_df, cols=c('muts', 'fusions', 'CNAs'), names_to = "type", values_to = "count")

# Plot histogram of counts within MSI classes
ggplot(pp_df_long, aes(x=MSI, y=count, color=type)) +
  geom_boxplot(alpha = 0.7, notch=TRUE) +
  labs(x='', y = "Count per Sample") 
  #theme_minimal() 


# Plot histogram of counts per sample
ggplot(pp_df_long %>% filter(!is.na(CodeSeizure)), aes(x=type, y=count, col=CodeSeizure)) +
  geom_boxplot(alpha = 0.7, notch=FALSE) +
  labs(x='', y = "Count per Sample")+
  theme_minimal()  + scale_color_discrete(labels=c("No Seizure","Pre-Diagnosis","Post-Diagnosis"))+scale_x_discrete(labels=c("Copy Number","Fusions","Mutations"))+labs(colour = "Seizure Incidence", x=c("Variation Type"))

ggsave(file="allCN_mut_fus.png", width=6, height=4, dpi=300)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(pp_df, aes(x=fusions, y=muts)) + geom_point(color = "black", alpha = 0.7) + labs(x='', y = "Count per Sample")+ theme_minimal()

cor.test(pp_df$CNAs, pp_df$fusions, method = 'spearman')

cor_matrix = cor(pp_df %>% select(c(TMB_numeric,kps_at_diagnosis,Seizure.freq.AVG.over.1st.6.mo, 'muts', 'fusions','CNAs','age')), use = 'pairwise.complete.obs', method = "spearman")
ggcorrplot(cor_matrix, type = 'lower',method = 'square',lab = TRUE,hc.order = TRUE)

ggsave("correlation_GBM.png")

cor_test =cor(pp_df %>% select(c(TMB_numeric,'muts', 'fusions','CNAs','age', 'ki_67')) %>% dplyr::mutate(ki_67=as.numeric(ki_67)),use = 'pairwise.complete.obs',method = "spearman")




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
library(gpairs)
library(GGally)
cont_df = pp_df %>% select(c(TMB_numeric,kps_at_diagnosis,Seizure.freq.AVG.over.1st.6.mo, 'muts', 'fusions','CNAs','age', 'ki_67', CodeSeizure, mgmt_methylated, sex)) %>% rename(TMB=TMB_numeric, KPS=kps_at_diagnosis, Seizure.freq = Seizure.freq.AVG.over.1st.6.mo, Mutations=muts, Fusions=fusions, Age=age, Sex=sex) %>% mutate(ki_67=as.numeric(ki_67), Sex=as.factor(Sex)) %>% filter(!is.na(TMB))
#chart.Correlation(cont_df, method = "s")

ggpairs(cont_df %>% filter(mgmt_methylated!="Indeterminate", !is.na(CodeSeizure)) %>% filter(Sex==0) %>% select(-c("Sex", "mgmt_methylated")), ggplot2::aes(colour=as.factor(CodeSeizure)),upper = list(continuous = wrap("cor", size = 2.5, method="spearman")))
ggsave("CorrelationContinuous_GBM_Seizure_Sex0.png", dpi=300, width=10, height=10)

ggpairs(cont_df %>% filter(mgmt_methylated!="Indeterminate", !is.na(CodeSeizure)) %>% filter(Sex==1) %>% select(-c("Sex", "mgmt_methylated")), ggplot2::aes(colour=as.factor(CodeSeizure)),upper = list(continuous = wrap("cor", size = 2.5, method="spearman")))
ggsave("CorrelationContinuous_GBM_Seizure_Sex1.png", dpi=300, width=10, height=10)
#gpairs(cont_df)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
ggpairs(cont_df %>% filter(mgmt_methylated=="0", !is.na(CodeSeizure)) %>% select(-c("Sex", "mgmt_methylated")), ggplot2::aes(colour=as.factor(CodeSeizure)),upper = list(continuous = wrap("cor", size = 2.5, method="spearman")))
ggsave("CorrelationContinuous_GBM_Seizure_mgmt0.png", dpi=300, width=10, height=10)

ggpairs(cont_df %>% filter(mgmt_methylated=="1", !is.na(CodeSeizure)) %>% select(-c("Sex", "mgmt_methylated")), ggplot2::aes(colour=as.factor(CodeSeizure)),upper = list(continuous = wrap("cor", size = 2.5, method="spearman")))
ggsave("CorrelationContinuous_GBM_Seizure_mgmt1.png", dpi=300, width=10, height=10)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Add singleton/doubleton

fusion_cols = pp_df %>% select(contains(".Fusion"))
singleton_fusions = names(fusion_cols[colSums(pp_df %>% select(contains(".Fusion")) != 0)==1])
pp_df$singleton_fusion = rowSums(pp_df %>% select(all_of(singleton_fusions))!=0)
doubleton_fusions = names(fusion_cols[colSums(pp_df %>% select(contains(".Fusion")) != 0)==2])
pp_df$doubleton_fusion = rowSums(pp_df %>% select(all_of(doubleton_fusions))!=0)

mutation_cols = pp_df %>% select(contains("Mutation"))
singleton_mutations = names(mutation_cols[colSums(mutation_cols != 0)==1])
pp_df$singleton_mutation = rowSums(pp_df %>% select(all_of(singleton_mutations))!=0)
doubleton_mutations = names(mutation_cols[colSums(mutation_cols != 0)==2])
pp_df$doubleton_mutation = rowSums(pp_df %>% select(all_of(doubleton_mutations))!=0)

cna_cols = pp_df %>% select(contains("CNA"))
singleton_cnas = names(cna_cols[colSums(cna_cols != 0)==1])
pp_df$singleton_cnas = rowSums(pp_df %>% select(all_of(singleton_cnas))!=0)
doubleton_cnas = names(cna_cols[colSums(cna_cols != 0)==2])
pp_df$doubleton_cnas = rowSums(pp_df %>% select(all_of(doubleton_cnas))!=0)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
library(PerformanceAnalytics)
library(gpairs)
singletondoubleton_df = pp_df %>% select(c('muts',doubleton_mutation, singleton_mutation, 'fusions',doubleton_fusion,singleton_fusion, 'CNAs', doubleton_cnas, singleton_cnas, kps_at_diagnosis, CodeSeizure, Seizure.freq.AVG.over.1st.6.mo, sex, age, mgmt_methylated)) %>% rename(KPS=kps_at_diagnosis, Sex=sex, Age=age, Mutations=muts, Fusions=fusions, Seizure.freq=Seizure.freq.AVG.over.1st.6.mo) %>% mutate(rare_mut=singleton_mutation+doubleton_mutation, rare_fus = doubleton_fusion+singleton_fusion, rare_cna=doubleton_cnas+singleton_cnas) %>% mutate(rare_fus=rare_fus/sum(rare_fus), rare_mut=rare_mut/sum(rare_mut), rare_cna=rare_cna/sum(rare_cna), rare_all = rare_cna+rare_mut+rare_fus) %>% select(rare_all, rare_mut, rare_fus, rare_cna, Age, Sex, KPS, Seizure.freq, CodeSeizure,mgmt_methylated)



#x=ggpairs(singletondoubleton_df %>% filter(Sex==1, !is.na(CodeSeizure)) %>% select(-c("Sex", "mgmt_methylated")), ggplot2::aes(colour=as.factor(CodeSeizure)),upper = list(continuous = wrap("cor", size = 2.5, method="spearman")))
#ggpairs(mtcars, columns = 2:4, ggplot2::aes(colour=as.character(am)))
#ggsave("rarevars_sex1CorrelationContinuous_GBM_Seizure.png",plot=x, dpi=300, width=10, height=10)




## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
dat = pp_df %>% select(c(LOH, muts)) %>% dplyr::mutate(mut_high=muts>2) %>% na.omit()
table(dat %>% select(-muts))
fisher.test(table(dat %>% select(-muts)))

dat = pp_df %>% select(c(LOH, CNAs)) %>% dplyr::mutate(CNA_any=CNAs>0) %>% na.omit()
table(dat %>% select(-CNAs))
fisher.test(table(dat %>% select(-CNAs)))

dat = pp_df %>% select(c(LOH, fusions)) %>% dplyr::mutate(fus_any=fusions>0) %>% na.omit()
table(dat %>% select(-fusions))
fisher.test(table(dat %>% select(-fusions)))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
dat = pp_df %>% select(c(MSI, muts)) %>% dplyr::mutate(mut_any=muts>2) %>% na.omit()
table(dat %>% select(-muts))
fisher.test(table(dat %>% select(-muts)))

dat = pp_df %>% select(c(MSI, CNAs)) %>% dplyr::mutate(CNA_any=CNAs>0) %>% na.omit()
table(dat %>% select(-CNAs))
fisher.test(table(dat %>% select(-CNAs)))

dat = pp_df %>% select(c(MSI, fusions)) %>% dplyr::mutate(fus_any=fusions>0) %>% na.omit()
table(dat %>% select(-fusions))
fisher.test(table(dat %>% select(-fusions)))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
dat = pp_df %>% select(c(Seizure.presentation.ANY, muts)) %>% dplyr::mutate(mut_any=muts>0) %>% na.omit()
table(dat %>% select(-muts))
fisher.test(table(dat %>% select(-muts)))

dat = pp_df %>% select(c(Seizure.presentation.ANY, CNAs)) %>% dplyr::mutate(CNA_any=CNAs>0) %>% na.omit()
table(dat %>% select(-CNAs))
fisher.test(table(dat %>% select(-CNAs)))

dat = pp_df %>% select(c(Seizure.presentation.ANY, fusions)) %>% dplyr::mutate(fus_any=fusions>0) %>% na.omit()
table(dat %>% select(-fusions))
fisher.test(table(dat %>% select(-fusions)))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, EGFR.CNA.CND,EGFR.Fusion )) %>% dplyr::mutate(egfr_any=(EGFR.CNA.CND!=0)|(EGFR.Fusion!=0)) %>% na.omit()
dat_table = table(dat %>% select(-c(EGFR.CNA.CND, EGFR.Fusion)))
fisher.test(dat_table)

dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, EGFR.Fusion )) %>% dplyr::mutate(egfr_any=EGFR.Fusion!=0) %>% na.omit()
table(dat %>% select(-EGFR.Fusion))
fisher.test(table(dat %>% select(-EGFR.Fusion)))

dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, EGFR.Mutation )) %>% dplyr::mutate(egfr_any=EGFR.Mutation!=0) %>% na.omit()
table(dat %>% select(-EGFR.Mutation))
fisher.test(table(dat %>% select(-EGFR.Mutation)))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
fus_singletons = colSums(Caris %>% select(contains("Fusion")) != 0)==1
mut_singletons = colSums(Caris %>% select(contains("Mutation")) != 0)==1
CN_singletons = colSums(Caris %>% select(contains("CNA")) != 0)==1

fus_doubletons = colSums(Caris %>% select(contains("Fusion")) != 0)==2
mut_doubletons = colSums(Caris %>% select(contains("Mutation")) != 0)==2
CN_doubletons = colSums(Caris %>% select(contains("CNA")) != 0)==2

  
sum(fus_doubletons)
sum(mut_doubletons)
sum(CN_doubletons)

sum(fus_singletons)
sum(mut_singletons)
sum(CN_singletons)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
Caris_clinical$TMB_numeric = as.numeric(substr(Caris_clinical$TMB,1,2))

fusion_cols = Caris_clinical %>% select(contains(".Fusion"))
singleton_fusions = names(fusion_cols[colSums(Caris_clinical %>% select(contains(".Fusion")) != 0)==1])
Caris_clinical$singleton_fusion = rowSums(Caris_clinical %>% select(all_of(singleton_fusions))!=0)
doubleton_fusions = names(fusion_cols[colSums(Caris_clinical %>% select(contains(".Fusion")) != 0)==2])
Caris_clinical$doubleton_fusion = rowSums(Caris_clinical %>% select(all_of(doubleton_fusions))!=0)

mutation_cols = Caris_clinical %>% select(contains("Mutation"))
singleton_mutations = names(mutation_cols[colSums(mutation_cols != 0)==1])
Caris_clinical$singleton_mutation = rowSums(Caris_clinical %>% select(all_of(singleton_mutations))!=0)
doubleton_mutations = names(mutation_cols[colSums(mutation_cols != 0)==2])
Caris_clinical$doubleton_mutation = rowSums(Caris_clinical %>% select(all_of(doubleton_mutations))!=0)

cna_cols = Caris_clinical %>% select(contains("CNA"))
singleton_cnas = names(cna_cols[colSums(cna_cols != 0)==1])
Caris_clinical$singleton_cnas = rowSums(Caris_clinical %>% select(all_of(singleton_cnas))!=0)
doubleton_cnas = names(cna_cols[colSums(cna_cols != 0)==2])
Caris_clinical$doubleton_cnas = rowSums(Caris_clinical %>% select(all_of(doubleton_cnas))!=0)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rare CNAs
dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, singleton_cnas)) %>% dplyr::mutate(singleton_cnas=singleton_cnas!=0) %>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)

dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, doubleton_cnas)) %>% dplyr::mutate(doubleton_cnas=doubleton_cnas!=0) %>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rare fusions
dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, singleton_fusion)) %>% dplyr::mutate(singleton_fusion=singleton_fusion!=0) %>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)

dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, doubleton_fusion)) %>% dplyr::mutate(doubleton_fusion=doubleton_fusion!=0) %>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Any singletons
dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, singleton_fusion, singleton_mutation, singleton_cnas)) %>% dplyr::mutate(singletons=(singleton_fusion!=0|singleton_mutation!=0|singleton_cnas!=0)) %>% select(Seizure.presentation.ANY,singletons) %>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)

#Any doubletons
dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, doubleton_fusion, doubleton_mutation, doubleton_cnas)) %>% dplyr::mutate(doubletons=(doubleton_fusion!=0|doubleton_mutation!=0|doubleton_cnas!=0)) %>% select(Seizure.presentation.ANY,doubletons)%>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rare mutations
dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, singleton_mutation)) %>% dplyr::mutate(singleton_mutation=singleton_mutation!=0) %>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)

dat = Caris_clinical %>% select(c(Seizure.presentation.ANY, doubleton_mutation)) %>% dplyr::mutate(doubleton_mutation=doubleton_mutation!=0) %>% na.omit()
dat_table = table(dat)
print(dat_table)
fisher.test(dat_table)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
multivariable_lm = glm(Seizure.presentation.ONLY.after.surgery ~ 
                     age + 
                     sex +
                     mgmt_methylated +
                     egfr_amplified + 
                     TMB_numeric + 
                     MSI + 
                     LOH, 
                     data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate"), 
                     family = "binomial")

summary(multivariable_lm)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------

multivariable_lm = glm(Seizure.presentation..0.No..1.Yes. ~ 
                     #age + 
                     sex + 
                     # age:sex +
                     mgmt_methylated +
                     #egfr_amplified +
                     #TMB_numeric + 
                     #MSI + 
                     #LOH + 
                       1
                      
                     , data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate"), 
                     family = "binomial")

summary(multivariable_lm)

multivariable_lm = glm(Seizure.presentation.ANY ~ 
                     #age + 
                     sex +
                     mgmt_methylated +
                     egfr_amplified +
                     TMB_numeric + 
                     #MSI + 
                     #LOH 
                     1
                      
                     , data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate"), 
                     family = "binomial")

summary(multivariable_lm)

multivariable_lm = glm(Seizure.presentation.ONLY.after.surgery ~ 
                     age + 
                     sex +
                     mgmt_methylated +
                     egfr_amplified +
                     #TMB_numeric + 
                     #MSI + 
                     #LOH + 
                     1, data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate"), 
                     family = "binomial")

summary(multivariable_lm)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
multivariable_lm = glm(Seizure.presentation..0.No..1.Yes. ~ 
                     #age + 
                     #doubleton_mutation +
                     #LOH+
                     #MSI +
                     #TMB_numeric +
    
                     #doubleton_cnas + 
                     #doubleton_fusion + 
                     #singleton_fusion + 
                     #singleton_cnas + 
                     #singleton_mutation + 
                     mgmt_methylated +
                     #egfr_amplified +
                     #loc_parietal + 
                     #loc_occipital + 
                     #loc_other
                     
                     singleton_cnas:doubleton_mutation +
                     1, data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate"), family = "binomial")
                    # singleton_mutation + singleton_cnas + singleton_fusion 
                   

summary(multivariable_lm)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
multivariable_lm = glm( Seizure.presentation..0.No..1.Yes.~ 
                     sex+
                     #doubleton_mutation+
                     #doubleton_mutation:doubleton_cnas+ 
                     #doubleton_fusion + 
                     #singleton_fusion + 
                     #singleton_cnas +
                     #singleton_mutation +
                     singleton_cnas:singleton_mutation + 
                     #egfr_amplified +
                     mgmt_methylated #+
                     #loc_parietal + 
                     #loc_occipital 
                     # + loc_frontal + loc_temporal
                   , data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate"), family = "binomial")
                    # singleton_mutation + singleton_cnas + singleton_fusion 
                   

summary(multivariable_lm)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
multivariable_lm = glm(Seizure.presentation.ONLY.after.surgery~ 
                     age
                     + sex 
                    
                    # Pathology 
                    #+ as.numeric(ki_67) 
                    + as.numeric(p53_percent_reactivity)
                       
                    # Rare variants 
                     #+ doubleton_mutation
                     #+ doubleton_cnas
                     #+ doubleton_fusion  
                     #+ singleton_fusion  
                     #+ singleton_mutation
                     #+ singleton_cnas:singleton_mutation  
                    + singleton_cnas
                    
                    #+as.numeric(TERT.Mutation!=0) 
                    #+as.numeric(EGFR.CNA.CND!=0)
                    #+as.numeric(STK11.CNA.CND!=0)
                    #+as.numeric(PTEN.Mutation!=0)
                    #+as.numeric(TP53.Mutation...564!=0)
                    #+as.numeric(EGFR.Mutation!=0)
                    #+as.numeric(EGFR.Fusion!=0)
                    #+as.numeric(TIMM23B.Fusion!=0)
                    #+as.numeric(MEF2B.CNA.CND!=0)
                    
                    # Existing prognostic 
                    + egfr_amplified 
                    #+ mgmt_methylated 
                    
                    # Location
                    #+ loc_parietal  
                    #+ loc_occipital 
                    + loc_frontal
                    #+ loc_temporal
                   , data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate"), family = "binomial")
                    # singleton_mutation + singleton_cnas + singleton_fusion 
                   

summary(multivariable_lm)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
multivariable_lm = glm(Seizure.presentation..0.No..1.Yes.~ 
                     #age
                     + sex 
                    
                    # Pathology 
                    
                    #+ as.numeric(ki_67)
                    #+ as.numeric(p53_percent_reactivity) 
                       
                    # Rare variants 
                     + doubleton_mutation
                     #+ doubleton_cnas
                     #+ doubleton_fusion:doubleton_cnas
                     #+ singleton_fusion
                     + doubleton_mutation:doubleton_cnas
                     #+ singleton_cnas 

                    
                    #+as.numeric(TERT.Mutation!=0) 
                    #+as.numeric(EGFR.CNA.CND!=0)
                    #+as.numeric(STK11.CNA.CND!=0)
                    +as.numeric(PTEN.Mutation!=0)
                    #+as.numeric(TP53.Mutation...564!=0)
                    #+as.numeric(EGFR.Mutation!=0)
                    #+as.numeric(TIMM23B.Fusion!=0)
                    #+as.numeric(MEF2B.CNA.CND!=0)

                    # Existing prognostic 
                    #+ egfr_amplified
                    + mgmt_methylated 
                    
                    # Combined
                    #+ egfr_amplified:as.numeric(ki_67)
                    #+ loc_occipital:as.numeric(ki_67)
                    #+ as.numeric(PTEN.Mutation!=0):singleton_mutation

                     # Location
                    #+ loc_parietal  
                    + loc_occipital
                    #+ loc_frontal
                    #+ loc_temporal
                   , data = Caris_clinical %>% filter(mgmt_methylated!="Indeterminate", !is.na(p53_percent_reactivity)), family = "binomial")
                    # singleton_mutation + singleton_cnas + singleton_fusion 
                   

summary(multivariable_lm)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
library(stringr)
pp_df$HLA.A= ifelse(pp_df$HLA.A=='NA',NA,pp_df$HLA.A)
pp_df$HLA.B= ifelse(pp_df$HLA.B=='NA',NA,pp_df$HLA.B)
pp_df$HLA.C= ifelse(pp_df$HLA.C=='NA',NA,pp_df$HLA.C)

library(tidyr)
pp_df_longHLA <- pp_df |>
separate_wider_delim(HLA.A, delim = ",", names = c("HLA.A1", "HLA.A2"), too_few = "align_end")
pp_df_longHLA <- pp_df_longHLA |>
separate_wider_delim(HLA.B, delim = ",", names = c("HLA.B1", "HLA.B2"), too_few = "align_end")
pp_df_longHLA <- pp_df_longHLA |>
separate_wider_delim(HLA.C, delim = ",", names = c("HLA.C1", "HLA.C2"), too_few = "align_end")

pp_df_longHLA <- pp_df_longHLA |> tidyr::pivot_longer(cols=c(HLA.A1, HLA.A2, HLA.B1, HLA.B2, HLA.C1, HLA.C2), names_to = c('HLA'), values_to = "HLA.alleles")

pp_df_longHLA$HLA.alleles = sub('\'','', pp_df_longHLA$HLA.alleles)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------

library(forcats)
pp_df_longHLA <- pp_df_longHLA |> mutate(sex=as.factor(sex)) |> mutate(sex=lvls_revalue(sex,c('M','F')))
pp_df_longHLA <- pp_df_longHLA |> mutate(mgmt_methylated=as.factor(mgmt_methylated))
pp_df_longHLA <- pp_df_longHLA |> mutate(mgmt_methylated=lvls_revalue(mgmt_methylated,c('Not Methylated','Methylated','Indeterminate')))
pp_df_longHLA <- pp_df_longHLA |> rename(Sz_Pre_Diagnosis=Seizure.presentation..0.No..1.Yes.)


# Rename, only ever needed once per session
#pp_df_longHLA <- pp_df_longHLA |> mutate(Seizure.presentation=lvls_revalue(Seizure.presentation,c('None', 'Seizures')))
#pp_df_longHLA <- pp_df_longHLA |> mutate(Seizure.presentation=as.factor(Seizure.presentation))


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Highest frequency to lowest frequency HLA

for (i in c('A','B','C')){ ggplot(pp_df_longHLA %>% na.omit() %>% 
           filter(grepl(i, 
                        HLA.alleles)) %>% 
           filter(mgmt_methylated!='Indeterminate'), 
         aes(y=fct_infreq(HLA.alleles), fill=Pre_Diagnosis_Seizures)) + 
    geom_bar(position="stack") + 
    facet_grid(mgmt_methylated~sex)+ 
    ylab(paste0('HLA-',i,' Genotype'))
  ggsave(paste0('HLA_Genotype_Dist',i,'.png'), width = 6, height=8)
       }


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#Load package
library(tableone)

#Create a variable list which we want in Table 1
listVars <- c("age", "TMB_numeric", "muts", "fusions", "CNAs")

#Define categorical variables
catVars <- c("sex","Seizure.presentation..0.No..1.Yes.","Seizure.presentation.ONLY.after.surgery", "LOH", "MSI") 

pp_df_tb1 <- pp_df |> select(!contains('Mutation')) |> select(!ends_with('.Fusion')) |> select(!contains('CND')) |> select(!TMB) |> mutate(CodeSeizure=Seizure.presentation..0.No..1.Yes.*1 + Seizure.presentation.ONLY.after.surgery*2) |> select(!starts_with("Deidentif")) |> mutate(across(catVars, as.factor)) |> mutate(across(starts_with("loc"), as.factor)) |> mutate(across(c(CodeSeizure,race, mgmt_methylated, idh1_mutated), as.factor)) |> mutate(across(ki_67,as.numeric)) |> mutate(HLA.A=grepl(x=HLA.A,pattern="A*01:01"), HLA.B=grepl(x=HLA.B,pattern="B*08:01"), HLA.C=grepl(x=HLA.C,pattern="C*07:01")) 

# Create summary table using gtsummary
library(gtsummary)

#table1 <- CreateTableOne(data = pp_df_tb1 %>% select(!starts_with('Sei')) ,addOverall = TRUE,includeNA = TRUE, strata ='CodeSeizure' )
table1 <- tbl_summary(
    pp_df_tb1 %>% select(!starts_with('Sei')),
    by = CodeSeizure, # split table by group 
    ) %>%
  add_n() %>% # add column with total number of non-missing observations
  add_p() # test for a difference between groups

#print(table1)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
#table1 <- pp_df%>% select(catVars, listVars)
# Load the packages
#library(ReporteRs) # Not maintained
library(magrittr)
library(flextable)
library(officer)


doc <- read_docx()
# The script
doc %>% 
     body_add_flextable(table1 %>%
     as_flex_table(.) %>%
               theme_zebra( odd_body = "#DDDDDD", even_body = "#FFFFFF" ) ) %>%
     print(.,target = "table1.docx")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------

pp_df_longHLA <- pp_df_longHLA |> mutate(TMB_group=lvls_revalue(as.factor(ntile(TMB_numeric,2)), c('High','Low')))

# Rename, only ever needed once per session
#pp_df_longHLA <- pp_df_longHLA |> mutate(Seizure.presentation=lvls_revalue(Seizure.presentation,c('None', 'Seizures')))
#pp_df_longHLA <- pp_df_longHLA |> mutate(Seizure.presentation=as.factor(Seizure.presentation))
#pp_df_longHLA <- pp_df_longHLA |> rename(Pre_Diagnosis_Seizures=Seizure.presentation)


# Highest frequency to lowest frequency (fct_infreq())
for (i in c('A','B','C')){ ggplot(pp_df_longHLA %>% na.omit() %>% 
           filter(grepl(i, 
                        HLA.alleles)) %>% 
           filter(mgmt_methylated!='Indeterminate'), 
         aes(y=fct_infreq(HLA.alleles), fill=TMB_group)) + 
    geom_bar(position="stack") + 
    facet_grid(mgmt_methylated~sex)+ 
    ylab(paste0('HLA-',i,' Genotype'))
  ggsave(paste0('HLA_Genotype_TMB_',i,'.png'), width = 6, height=8)
       }



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------


pp_df_longHLA <- pp_df_longHLA |> mutate(CodeSeizure=Sz_Pre_Diagnosis*1 + Seizure.presentation.ONLY.after.surgery*2) |> mutate(CodeSeizure=as.factor(CodeSeizure))

binary_data <- model.matrix(~ . - 1, data = pp_df_longHLA %>% select(-id_cols[1:4]) %>% select(contains(".Fusion")) %>% mutate_all(as.factor))

# Calculate Jaccard distance
jaccard_dist <- daisy(binary_data, metric = "gower")

# Perform hierarchical clustering
hclust_result <- hclust(jaccard_dist, method = "complete")

# Create dendrogram
dend <- as.dendrogram(hclust_result)

# Sample labels and corresponding colors
sample_labels <- as.factor(pp_df_longHLA$CodeSeizure)
colors <- c("red", "blue", "green", "orange")

# Assign colors to labels
label_colors <- colors[as.numeric(sample_labels)][hclust_result$order]
labels_colors(dend) <- label_colors

# Plot dendrogram with colored labels
plot(dend, main = "Hierarchical Clustering Dendrogram (Mutation data)", xlab = "Samples", sub = NULL, horiz = FALSE, )


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Convert categorical data to binary indicators
pp_df_longHLA <- pp_df_longHLA |> mutate(CodeSeizure=Sz_Pre_Diagnosis*1 + Seizure.presentation.ONLY.after.surgery*2) |> mutate(CodeSeizure=as.factor(CodeSeizure))

binary_data <- model.matrix(~ . - 1, data = pp_df_longHLA %>% select(-id_cols[1:4]) %>% select(contains("Mutation")) %>% mutate_all(as.factor))

# Calculate Jaccard distance
jaccard_dist <- daisy(binary_data, metric = "gower")

# Perform hierarchical clustering
hclust_result <- hclust(jaccard_dist, method = "complete")

# Create dendrogram
dend <- as.dendrogram(hclust_result)

# Sample labels and corresponding colors
sample_labels <- as.factor(pp_df_longHLA$CodeSeizure)
colors <- c("red", "blue", "green", "orange")

# Assign colors to labels
label_colors <- colors[as.numeric(sample_labels)][hclust_result$order]
labels_colors(dend) <- label_colors

# Plot dendrogram with colored labels
plot(dend, main = "Hierarchical Clustering Dendrogram (Mutation data)", xlab = "Samples", sub = NULL, horiz = FALSE, )


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Convert categorical data to binary indicators
pp_df_longHLA <- pp_df_longHLA |> mutate(CodeSeizure=Sz_Pre_Diagnosis*1 + Seizure.presentation.ONLY.after.surgery*2) |> mutate(CodeSeizure=as.factor(CodeSeizure))

binary_data <- model.matrix(~ . - 1, data = pp_df_longHLA %>% select(-id_cols[1:4]) %>% select(contains("CND")) %>% mutate_all(as.factor))

# Calculate Jaccard distance
jaccard_dist <- daisy(binary_data, metric = "gower")

# Perform hierarchical clustering
hclust_result <- hclust(jaccard_dist, method = "complete")

# Create dendrogram
dend <- as.dendrogram(hclust_result)

# Sample labels and corresponding colors
sample_labels <- as.factor(pp_df_longHLA$CodeSeizure)
colors <- c("red", "blue", "green", "orange")

# Assign colors to labels
label_colors <- colors[as.numeric(sample_labels)][hclust_result$order]
labels_colors(dend) <- label_colors

# Plot dendrogram with colored labels
plot(dend, main = "Hierarchical Clustering Dendrogram (Mutation data)", xlab = "Samples", sub = NULL, horiz = FALSE, )


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Install and load necessary packages
#install.packages("dendextend")
library(dendextend)

# Convert categorical data to binary indicators
binary_data <- model.matrix(~ . - 1, data = as.data.frame(lapply(pp_df_longHLA %>% select(-id_cols[1:4]) %>% select(contains("Fusion")), as.factor)))

# Calculate Jaccard distance
jaccard_dist <- daisy(binary_data, metric = "gower")

# Perform hierarchical clustering
hclust_result <- hclust(jaccard_dist, method = "complete")

# Create dendrogram
dend <- as.dendrogram(hclust_result)

# Sample labels and corresponding colors
sample_labels <- as.factor(pp_df_longHLA$CodeSeizure)
colors <- c("red", "blue", "green", "orange")

# Assign colors to labels
label_colors <- colors[as.numeric(sample_labels)][hclust_result$order]
labels_colors(dend) <- label_colors

# Plot dendrogram with colored labels
plot(dend, main = "Hierarchical Clustering Dendrogram (Gene fusion data)", xlab = "Samples", sub = NULL, horiz = FALSE, )


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# Install and load necessary packages
#install.packages("randomForest")
library(randomForest)

pp_df <- left_join(pp_df, Caris_clinical %>% select(c(Deidentified.code,contains("ton"), chemo_tx_cytotoxic, chemo_tx_biol_targ)))

rare_cols <- grep("ton", names(pp_df), value = TRUE)
loc_columns <-  grep("loc", names(pp_df), value = TRUE)
mut_columns <-  grep("Mutation", names(pp_df), value = TRUE)
cn_columns <-  grep("CND", names(pp_df), value = TRUE)

pp_df <- pp_df %>% na.omit()

### Classification column 
pp_df$Class <- as.factor(pp_df$CodeSeizure)

equal_0 <- function(x){ return(x==0)}

### Columns for classification
'c(Class, p53_percent_reactivity, age, ki_67, muts, CNAs, fusions, TMB_numeric, sex, loc_columns, cn_top, fus_top, rare_cols)'

### Test/Train split 
split_index = sample(1:nrow(pp_df %>% select(-CodeSeizure)), 0.7*nrow(pp_df %>% select(-CodeSeizure)))

train_data = (pp_df %>%  select(-c(Deidentified.code,CodeSeizure)))[split_index, ]  %>% select(Class, chemo_tx_biol_targ,p53_percent_reactivity, age, ki_67, muts, CNAs, fusions, TMB_numeric, sex, rare_cols, loc_columns, names(fus_top)[1:2], names(cn_top)[1:5], names(mut_top)[1:5]) %>% mutate(across(c(names(fus_top)[1:2], names(cn_top)[1:5], names(mut_top)[1:5]), equal_0))
train_data$ki_67 =as.numeric(train_data$ki_67)

test_data = (pp_df %>% select(-c(Deidentified.code,CodeSeizure)))[-split_index, ]  %>% select(Class, chemo_tx_biol_targ, loc_columns, p53_percent_reactivity, age, ki_67, muts, CNAs, fusions, TMB_numeric, sex, rare_cols, names(fus_top)[1:2], names(cn_top)[1:5], names(mut_top)[1:5]) %>% mutate(across(c(names(fus_top)[1:2], names(cn_top)[1:5], names(mut_top)[1:5]), equal_0))

test_data$ki_67 =as.numeric(test_data$ki_67)


rf_model <- randomForest(Class ~ ., data = train_data, ntree = 10, max_depth=6, na.action=na.omit)
predictions <- predict(rf_model, test_data)
confusion_matrix <- table(predictions, test_data$Class)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(confusion_matrix)
print(paste("Accuracy:", round(accuracy, 2)))
print(importance(rf_model))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(61)
library(glmnet)
library(caret)

# this method controls everything about training
# we will just set up 10-fold cross validation
trctrl <- trainControl(method = "cv",number=10)

col_zeros = names(which(colSums(train_data==0) == nrow(train_data)))
# we will now train elastic net model
# it will try
enetFit <- train(Class~., data = train_data %>%na.omit() %>% select(!col_zeros), 
                 method = "glmnet",  
                 trControl=trctrl,
                 na.action = na.omit,
                 # alpha and lambda paramters to try
                 tuneGrid = data.frame(alpha=0.05,
                                       lambda=seq(0.05,0.9,0.01)))

# best alpha and lambda values by cross-validation accuracy
enetFit$bestTune

#test_data =test_data %>%na.omit() %>% select(!col_zeros)

class.res=predict(enetFit,test_data[,-1])
confusionMatrix(test_data$Class,class.res)$overall[1]

confusionMatrix(class.res, test_data$Class)

plot(varImp(enetFit),top=10)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(pp_df, aes(y=CNAs, x="Copy Number Alterations")) +
  geom_boxplot(fill = "skyblue", color = "black", alpha = 0.7) +
  labs(x='', y = "Frequency") + 
  theme_minimal() 
ggsave(file="CNAs.png", width=2, height=4, dpi=300)

