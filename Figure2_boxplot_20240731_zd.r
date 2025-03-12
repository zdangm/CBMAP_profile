rm(list=ls())
gc()
library(dplyr)
library(data.table)
library(readxl)
library(stringr)
library(ggplot2)
setwd("/data/projects/China_Brain_MultiOmics/profile/RIN_analysis/result/boxplot/")


################################### read  CBMAP data

brain_data = as.data.frame(read.csv("/data/shared_data/China_Brain_MultiOmics/sample_information/process/all_qc_sample_summary_raw_20240731.csv",check.names = F, fileEncoding = "GBK"))


#clean
brain_data[,'bank'] = brain_data[,'bank.x']

brain_data[which(brain_data[,'gender'] %in% c('M','男')),'sex_male'] = 1
brain_data[-which(brain_data[,'gender'] %in% c('M','男')),'sex_male'] = 0


brain_data[which(brain_data[,'age.x'] == '28天'),'age.x'] = 0.1 #28 days
brain_data[which(brain_data[,'age.x'] == '8个月'),'age.x'] = 0.7 #8 months

brain_data[which(brain_data[,'PMD'] == 0),'PMD'] = NA


brain_data[,'age'] = as.numeric(brain_data[,'age.x'])

brain_data[,'year'] = as.numeric(substr(as.character(brain_data[, 'year.x']),1,4))


#only adults
brain_data = brain_data[which(brain_data$age>=18), ]
#only PMD <=24
brain_data = brain_data[which(brain_data$PMD<=24 | is.na(brain_data$PMD)), ]


CBMAP = brain_data[,c('id','bank','sex_male','age','year',
                      'PMD','RIN',
                      # 'Tau_grade','Ab_grade',  #only 湘雅
                      # 'APOE_genotype', #only 浙大
                      "ADRP","ADRP-L","ADRP-M","ADRP-H",
                      "LBD","CVD","CVD_atherosclerosis",
                      "CVD_small_arteriosclerosis",
                      "CVD_cerebral_hemorrhage",
                      "CVD_cerebral_infarction",
                      "PART","LATE","ARTAG")]


#write
# write.csv(brain_df, '~/../ZLab/projects/China_Brain_MultiOmics/projects/study_profile/data/brain_df.csv',quote = F,row.names = F)
write.table(CBMAP, '/data/shared_data/China_Brain_MultiOmics/sample_information/process/brain_df_0731.txt',quote = F,row.names = F,sep = '\t',fileEncoding = "GBK")


#---figure 3---

rin_sample_CBMAP<-data.frame(id=CBMAP$id,RIN=CBMAP$RIN,PMI=CBMAP$PMD,age=CBMAP$age,bank=CBMAP$bank,source="CBMAP")

rin_sample_CBMAP$age<-as.numeric(rin_sample_CBMAP$age)
rin_sample_CBMAP$bank<-toupper(rin_sample_CBMAP$bank)
rin_sample_CBMAP$RIN<-as.numeric(rin_sample_CBMAP$RIN)


##################### process gtex ##############

# Read the expression from gtex
#exp_gtex <- as.data.frame(fread("/data/shared_data/gtex/exp/TPM_log_matrix/Brain_Frontal_Cortex_BA9.txt.gz", header = TRUE, sep = "\t",check.names = F))

# Read the RIN from gtex
rin_gtex <- as.data.frame(fread("/data/shared_data/gtex/phenotypes/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", header = TRUE, sep = "\t"))
rin_gtex <- rin_gtex[rin_gtex$SMTSD=="Brain - Frontal Cortex (BA9)",c("SAMPID","SMRIN")]
rin_gtex <- rin_gtex[!is.na(rin_gtex$SMRIN),]

rin_gtex_SAMPID<-as.data.frame(str_split(rin_gtex$SAMPID,"-",simplify = T))

########## check the frequency of individualID
rin_gtex$SAMPID<-paste0(rin_gtex_SAMPID$V1,"-",rin_gtex_SAMPID$V2)
rin_gtex<-unique(rin_gtex)
freq<-as.data.frame(table(rin_gtex$SAMPID))
freq[freq$Freq>1,]
rin_gtex$freq<-freq[match(rin_gtex$SAMPID,freq$Var1),"Freq"]
############ remove the the individual with frequency >1
rin_gtex<-rin_gtex[rin_gtex$freq==1,]

# Read the PMI and age from gtex
sampleinfo_gtex <- as.data.frame(fread("/data/shared_data/gtex/phenotypes/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt", header = TRUE, sep = "\t"))

############### fill the rin,PMI,age to the sample from gtex expression
rin_sample_gtex<-data.frame(id=rin_gtex$SAMPID)
rin_sample_gtex$RIN<-rin_gtex[match(rin_sample_gtex$id,rin_gtex$SAMPID),"SMRIN"]
rin_sample_gtex$PMI<-sampleinfo_gtex[match(rin_sample_gtex$id,sampleinfo_gtex$SUBJID),"TRISCHD"]/60
rin_sample_gtex$age<-sampleinfo_gtex[match(rin_sample_gtex$id,sampleinfo_gtex$SUBJID),"AGE"]
rin_sample_gtex$source<-"GTEx"

################################### process rosmap ################################################

#exp_rosmap<-as.data.frame(fread("/data/shared_data/ROSMAP/Processed/gene_expression/RNA_seq/bulk_brain/ROSMAP_RNAseq_FPKM_gene.tsv", header = TRUE, sep = "\t",check.names = F))

rin_rosmap<-as.data.frame(fread("/data/shared_data/ROSMAP/Processed/metadata/ROSMAP_assay_rnaSeq_metadata.csv",check.names = F))

sampleinfo_rosmap<-as.data.frame(fread("/data/shared_data/ROSMAP/Processed/metadata/ROSMAP_clinical.csv",check.names = F))

rosmap_id_info<-as.data.frame(fread("/data/shared_data/ROSMAP/Processed/metadata/ROSMAP_biospecimen_metadata.csv",check.names = F))
########## filter sample
rosmap_id_info<-rosmap_id_info[rosmap_id_info$tissue=="dorsolateral prefrontal cortex",]
rosmap_id_info<-rosmap_id_info[!str_starts(rosmap_id_info$specimenID,"RISK"),]

############# select sample
rin_sample_rosmap<-data.frame(specimenID=rosmap_id_info$specimenID)
#rin_sample_rosmap<-data.frame(str_split(rin_sample_rosmap$specimenID,"_",simplify = T))
#rin_sample_rosmap$specimenID<-paste0(rin_sample_rosmap$X1,"_",rin_sample_rosmap$X2)
#freq_rosmap<-as.data.frame(table(rin_sample_rosmap$specimenID)) ############## freq 492_120515 =3

rin_sample_rosmap$id<-rosmap_id_info[match(rin_sample_rosmap$specimenID,rosmap_id_info$specimenID),"individualID"]

############### fill the rin,PMI,age to the sample from rosmap expression
rin_sample_rosmap$RIN<-rin_rosmap[match(rin_sample_rosmap$specimenID,rin_rosmap$specimenID),"RIN"]
rin_sample_rosmap$PMI<-sampleinfo_rosmap[match(rin_sample_rosmap$id,sampleinfo_rosmap$individualID),"pmi"]
rin_sample_rosmap$age<-sampleinfo_rosmap[match(rin_sample_rosmap$id,sampleinfo_rosmap$individualID),"age_death"]
rin_sample_rosmap$age<-ifelse(rin_sample_rosmap$age=="90+","90",rin_sample_rosmap$age)
rin_sample_rosmap$age<-as.numeric(rin_sample_rosmap$age)
rin_sample_rosmap$assay<-rosmap_id_info[match(rin_sample_rosmap$specimenID,rosmap_id_info$specimenID),"assay"]
rin_sample_rosmap<-rin_sample_rosmap[rin_sample_rosmap$id!="",]
rin_sample_rosmap<-rin_sample_rosmap[!is.na(rin_sample_rosmap$RIN),]

########## check the frequency of individualID
rin_sample_rosmap<-rin_sample_rosmap[,-1*which(colnames(rin_sample_rosmap)=="specimenID")]
rin_sample_rosmap<-unique(rin_sample_rosmap)
freq_rosmap<-as.data.frame(table(rin_sample_rosmap$id))
rin_sample_rosmap$freq<-freq_rosmap[match(rin_sample_rosmap$id,freq_rosmap$Var1),"Freq"]

############ remove the the individual with frequency >1
rin_sample_rosmap<-rin_sample_rosmap[rin_sample_rosmap$freq==1,]
rin_sample_rosmap$source<-"ROSMAP"



################################################## plot ##################
data<-rbind(rin_sample_CBMAP[,c("source","RIN","PMI","age")],rin_sample_gtex[,c("source","RIN","PMI","age")],rin_sample_rosmap[,c("source","RIN","PMI","age")])
data$age<-as.numeric(data$age)
data$RIN<-as.numeric(data$RIN)



p1<-ggplot(data, aes(x = source, y = as.numeric(RIN), fill = factor(source))) +
  geom_violin(trim = FALSE , alpha = 0.5,bw=1) +
  geom_boxplot(width=0.1,position=position_dodge(1)) +
  #geom_jitter(width = 0.2, size = 0.01, alpha = 0.6) +
  labs(title = "RIN",
       x="",
       y = "RIN value",
       tag = 'a') +
  # theme_minimal()+
  scale_fill_brewer(palette = "Set2")+ 

  theme_classic() +  # 去掉底纹网格
  theme(legend.position = "none",
    axis.line = element_line(color = "black"))  # 坐标轴用黑色线条



p2<-ggplot(data, aes(x = source, y = PMI, fill = factor(source))) +
  geom_violin(trim = FALSE , alpha = 0.5,bw=1) +
  geom_boxplot(width=0.1,position=position_dodge(1)) +
  #geom_jitter(width = 0.2, size = 0.01, alpha = 0.6) +
  labs(title = "PMI",
       x="",
       y = "PMI value",
       tag = 'b') +
  scale_fill_brewer(palette = "Set2")+ 
  ylim(0, 50)+
  theme_classic() +  # 去掉底纹网格
  theme(legend.position = "none",
        axis.line = element_line(color = "black"))  # 坐标轴用黑色线条

# 设置 y 轴范围


data_cbmap <- data[data$source == "CBMAP", ]
data_gtex <- data[data$source == "GTEx", ]
data_rosmap <- data[data$source == "ROSMAP", ]

p3<-ggplot() +
  geom_violin(data = data_cbmap, aes(x = source, y = age, fill = factor(source)), width = 0.8, trim = FALSE, alpha = 0.5,bw=3) +
  geom_violin(data = data_gtex, aes(x = source, y = age, fill = factor(source)), width = 0.8, trim = FALSE, alpha = 0.5,bw=2) +
  geom_violin(data = data_rosmap, aes(x = source, y = age, fill = factor(source)), width = 0.8, trim = FALSE, alpha = 0.5,bw=2) +
  geom_boxplot(data = data, aes(x = source, y = age, fill = factor(source)),width=0.1,position=position_dodge(1)) +
  #geom_jitter(data = data, aes(x = x, y = y), width = 0.2, size = 1, alpha = 0.6) +
  labs(title = "AGE",
       x="",
       y = "AGE value",
       fill = "Project",
       tag = 'c') +
  scale_fill_brewer(palette = "Set2")+
  theme_classic() +  # 去掉底纹网格
  theme(axis.line = element_line(color = "black"))  # 坐标轴用黑色线条



#p3<-ggplot(data, aes(x = source, y = age, fill = factor(source))) +
#  geom_violin(trim = FALSE , alpha = 0.5,bw=1) +
#  geom_boxplot(width=0.1,position=position_dodge(1)) +
#  #geom_jitter(width = 0.2, size = 0.01, alpha = 0.6) +
#  labs(title = "AGE",
#       x="",
#       y = "AGE value",
#       fill = "Project Source") +
#  theme_minimal()+
#  scale_fill_brewer(palette = "Set2")




p4<-ggplot(rin_sample_CBMAP, aes(x = bank, y = as.numeric(RIN), fill = factor(bank))) +
  geom_violin(trim = FALSE , alpha = 0.5,bw=1) +
  geom_boxplot(width=0.1,position=position_dodge(1)) +
  #geom_jitter(width = 0.2, size = 0.01, alpha = 0.6) +
  labs(title = "RIN",
       x="",
       y = "RIN value",
       tag = 'd') +
  scale_fill_brewer(palette = "Set4")+ 
  theme_classic() +  # 去掉底纹网格
  theme(legend.position = "none",
        axis.line = element_line(color = "black"))  # 坐标轴用黑色线条




p5<-ggplot(rin_sample_CBMAP, aes(x = bank, y = PMI, fill = factor(bank))) +
  geom_violin(trim = FALSE , alpha = 0.5,bw=1) +
  geom_boxplot(width=0.1,position=position_dodge(1)) +
  #geom_jitter(width = 0.2, size = 0.01, alpha = 0.6) +
  labs(title = "PMI",
       x="",
       y = "PMI value",
       tag = 'e') +
  scale_fill_brewer(palette = "Set4")+ 
  ylim(0, 50)+ # 设置 y 轴范围
  theme_classic() +  # 去掉底纹网格
  theme(legend.position = "none",
        axis.line = element_line(color = "black"))  # 坐标轴用黑色线条





p6<-ggplot(rin_sample_CBMAP, aes(x = bank, y = age, fill = factor(bank))) +
  geom_violin(trim = FALSE , alpha = 0.5,bw=4) +
  geom_boxplot(width=0.1,position=position_dodge(1)) +
  #geom_jitter(width = 0.2, size =0.01, alpha = 0.6) +
  labs(title = "AGE",
       x="",
       y = "AGE value",
       fill = "CBMAP Center",
       tag = 'f') +
  scale_fill_brewer(palette = "Set4")+
  theme_classic() +  # 去掉底纹网格
  theme(axis.line = element_line(color = "black"))  # 坐标轴用黑色线条




pdf('/data/projects/China_Brain_MultiOmics/profile/RIN_analysis/result/boxplot/RIN_boxplot_20240731.pdf',height = 8, width = 12)

p1+p2+p3+p4+p5+p6

dev.off()



#write.csv(result,"20240618_zju_firstBatch_summary_phospho.csv", row.names = F)








