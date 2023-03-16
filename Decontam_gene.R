library(momr)
library(data.table)
library(ggfortify)
library(ggsci)
library(dplyr)
library("textshape")
c_m<-read.csv("/home/junwoojo/다운로드/All_data/Colon.tsv",sep = "\t",header = T)
e_m<-read.csv("/home/junwoojo/다운로드/All_data/esophagus.tsv",sep = "\t",header = T)
li_m<-read.csv("/home/junwoojo/다운로드/All_data/Liver.tsv",sep = "\t",header = T)
lu_m<-read.csv("/home/junwoojo/다운로드/All_data/Lung.tsv",sep = "\t",header = T)
s_m<-read.csv("/home/junwoojo/다운로드/All_data/Stomach.tsv",sep = "\t",header = T)
h_m<-read.csv("/home/junwoojo/다운로드/All_data//head_neck.tsv",sep = "\t",header = T)



c_b_m<-read.csv("/home/junwoojo/다운로드/All_blood_data/colon_blood.tsv",sep = "\t",header = T)
e_b_m<-read.csv("/home/junwoojo/다운로드/All_blood_data/esophagus.tsv",sep = "\t",header = T)
li_b_m<-read.csv("/home/junwoojo/다운로드/All_blood_data/liver_blood.tsv",sep = "\t",header = T)
lu_b_m<-read.csv("/home/junwoojo/다운로드/All_blood_data/lung_blood.tsv",sep = "\t",header = T)
s_b_m<-read.csv("/home/junwoojo/다운로드/All_blood_data/stomach_blood.tsv",sep = "\t",header = T)
h_b_m<-read.csv("/home/junwoojo/다운로드/All_blood_data/head_neck.tsv",sep = "\t",header = T)







metadata<-read.csv('/home/junwoojo/all_metadata',sep="\t")


c_m <- c_m[,-1]
e_m <- e_m[,-1]
li_m <- li_m[,-1]
lu_m <- lu_m[,-1]
s_m <- s_m[,-1]
h_m <- h_m[,-1]

c_b_m <- c_b_m[,-1]
e_b_m <- e_b_m[,-1]
li_b_m <- li_b_m[,-1]
lu_b_m <- lu_b_m[,-1]
s_b_m <- s_b_m[,-1]
h_b_m <- h_b_m[,-1]


final_matrix<-do.call(cbind, list(c_m,e_m,li_m,lu_m,s_m,h_m))
final_matrix<-final_matrix[,-c(629)]

final_blood_matrix<-do.call(cbind, list(c_b_m,e_b_m,li_b_m,lu_b_m,s_b_m,h_b_m))


 
final_matrix_row<-final_matrix%>%rownames()
final_matrix_for_msp_miner<-cbind(final_matrix_row,final_matrix)
dim(final_matrix_for_msp_miner)

write_tsv(final_matrix_for_msp_miner,file="~/final_matrix_for_msp_miner",quote = 'none',col_names = T,)
final_matrix_for_msp_miner%>%dim()

#write.table(c_m[colon_tissue_enriched_gene,],file="~/MSPminer/for_msp_decontam_colon",quote = F,row.names = T,sep = "\t",col.names = F)





# final_matrix_row<-final_matrix%>%rownames()
# final_matrix<-cbind(final_matrix_row,final_matrix)
# final_matrix<-cbind(final_matrix_row,final_matrix)
# final_matrix<-final_matrix[enriched_all_gene,]







######################################################################################################################
######################################################################################################################
######################################################################################################################
######################################################################################################################






final_matrix<-cbind(final_matrix,final_blood_matrix)


length_data<-read.csv('/home/junwoojo/다운로드//gene_catalogue_lite_annotation',sep="\t",header=F)
head(length_data)

gene_name<-length_data[,1]%>%as.character()
gene_length<-length_data[,2]
length_data=c(gene_name=gene_length)
names(length_data)=gene_name

length(gene_name)
length(rownames(final_matrix))
rownames(final_matrix)=gene_name




dim(final_matrix)

final_last_matrix<-t(final_matrix)
final_last_matrix[2>final_last_matrix]<-0
final_last_matrix[final_last_matrix>=2]<-1



metadata_name<-gsub('-','.',metadata[,1])
metadata_name<-gsub('_','.',metadata_name)
metadata$X<-paste0(metadata_name,'.d')




last_mat_name<-final_last_matrix%>%rownames()
last_mat_name<-gsub('-','.',last_mat_name)
last_mat_name<-gsub('_','.',last_mat_name)
rownames(final_last_matrix)<-last_mat_name


names(metadata)[1]<-'file_name'

final_norm_mat_rowname<-rownames(final_last_matrix)

zxc<-final_norm_mat_rowname%>%as.data.frame()
colnames(zxc)<-'file_name'

qqq<-gsub('\\.bam','',zxc$file_name)
qqq<-gsub('TCGA\\.','',qqq)
zxc$file_name<-qqq
mzc<-gsub('\\.bam','',metadata$file_name)
mzc<-gsub('TCGA\\.','',mzc)
metadata$file_name<-mzc



right_metadata<-left_join(zxc,metadata,by='file_name')
right_metadata<-right_metadata[,-1]


last_abundance_matrix<-cbind(final_last_matrix,right_metadata)

dim(last_abundance_matrix)
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#################################colon_tissue_dataframe######################################
library(reshape)
last_abundance_matrix$primary_site%>%table()


colon_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Colorectal")
colon_tissue_matrix<-colon_total_matrix%>%filter(sample_type== "Primary Tumor"|sample_type=="Solid Tissue Normal")
#colon_tissue_matrix$subject_id%>%group_by(subject_id)%>%

subject_id<-colon_tissue_matrix$subject_id
colon_tissue_matrix_value<-colon_tissue_matrix[,1:806882]
colon_tissue_matrix<-cbind(colon_tissue_matrix_value,subject_id)

colon_tissue_matrix<-colon_tissue_matrix%>%as.data.table()

gene_name<-colon_tissue_matrix_value%>%colnames()
unique_id<-unique(colon_tissue_matrix$subject_id)
length(colon_tissue_matrix$subject_id)


first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-colon_tissue_matrix_value[colon_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}


default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/colon_tissue.tsv",quote = F,row.names = FALSE)

#####################################################################################


library(reshape)
colon_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Colorectal")

colon_tissue_matrix<-colon_total_matrix%>%filter(sample_type=="Blood Derived Normal" )

subject_id<-colon_tissue_matrix$subject_id
colon_tissue_matrix_value<-colon_tissue_matrix[,1:806882]
colon_tissue_matrix<-cbind(colon_tissue_matrix_value,subject_id)

colon_tissue_matrix<-colon_tissue_matrix%>%as.data.table()
colon_tissue_matrix[,806883]%>%table()

gene_name<-colon_tissue_matrix_value%>%colnames()
unique_id<-unique(colon_tissue_matrix$subject_id)
first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-colon_tissue_matrix_value[colon_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}

default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/colon_blood.tsv",quote = F,row.names = FALSE)

#########################################################################################################################
#########################################################################################################################















#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#################################Esophagus_tissue_dataframe######################################
Esophagus_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Esophagus")
Esophagus_tissue_matrix<-Esophagus_total_matrix%>%filter(sample_type== "Primary Tumor"|sample_type=="Solid Tissue Normal")
#Esophagus_tissue_matrix$subject_id%>%group_by(subject_id)%>%

subject_id<-Esophagus_tissue_matrix$subject_id
Esophagus_tissue_matrix_value<-Esophagus_tissue_matrix[,1:806882]
Esophagus_tissue_matrix<-cbind(Esophagus_tissue_matrix_value,subject_id)

Esophagus_tissue_matrix<-Esophagus_tissue_matrix%>%as.data.table()

gene_name<-Esophagus_tissue_matrix_value%>%colnames()
unique_id<-unique(Esophagus_tissue_matrix$subject_id)
length(Esophagus_tissue_matrix$subject_id)


first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-Esophagus_tissue_matrix_value[Esophagus_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}


default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/Esophagus_tissue.tsv",quote = F,row.names = FALSE)

#####################################################################################
library(reshape)
Esophagus_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Esophagus")

Esophagus_tissue_matrix<-Esophagus_total_matrix%>%filter(sample_type=="Blood Derived Normal" )

subject_id<-Esophagus_tissue_matrix$subject_id
Esophagus_tissue_matrix_value<-Esophagus_tissue_matrix[,1:806882]
Esophagus_tissue_matrix<-cbind(Esophagus_tissue_matrix_value,subject_id)

Esophagus_tissue_matrix<-Esophagus_tissue_matrix%>%as.data.table()
Esophagus_tissue_matrix[,806883]%>%table()

gene_name<-Esophagus_tissue_matrix_value%>%colnames()
unique_id<-unique(Esophagus_tissue_matrix$subject_id)
first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-Esophagus_tissue_matrix_value[Esophagus_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}

default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/Esophagus_blood.tsv",quote = F,row.names = FALSE)
####################################################################################




#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#################################liver_tissue_dataframe######################################
liver_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Liver")
liver_tissue_matrix<-liver_total_matrix%>%filter(sample_type== "Primary Tumor"|sample_type=="Solid Tissue Normal")
#liver_tissue_matrix$subject_id%>%group_by(subject_id)%>%

subject_id<-liver_tissue_matrix$subject_id
liver_tissue_matrix_value<-liver_tissue_matrix[,1:806882]
liver_tissue_matrix<-cbind(liver_tissue_matrix_value,subject_id)

liver_tissue_matrix<-liver_tissue_matrix%>%as.data.table()

gene_name<-liver_tissue_matrix_value%>%colnames()
unique_id<-unique(liver_tissue_matrix$subject_id)
length(liver_tissue_matrix$subject_id)


first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-liver_tissue_matrix_value[liver_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}


default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/liver_tissue.tsv",quote = F,row.names = FALSE)

#####################################################################################
library(reshape)
liver_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Liver")

liver_tissue_matrix<-liver_total_matrix%>%filter(sample_type=="Blood Derived Normal" )

subject_id<-liver_tissue_matrix$subject_id
liver_tissue_matrix_value<-liver_tissue_matrix[,1:806882]
liver_tissue_matrix<-cbind(liver_tissue_matrix_value,subject_id)

liver_tissue_matrix<-liver_tissue_matrix%>%as.data.table()
liver_tissue_matrix[,806883]%>%table()

gene_name<-liver_tissue_matrix_value%>%colnames()
unique_id<-unique(liver_tissue_matrix$subject_id)
first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-liver_tissue_matrix_value[liver_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}

default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/liver_blood.tsv",quote = F,row.names = FALSE)
####################################################################################





#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#################################lung_tissue_dataframe######################################
lung_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Lung")
lung_tissue_matrix<-lung_total_matrix%>%filter(sample_type== "Primary Tumor"|sample_type=="Solid Tissue Normal")
#lung_tissue_matrix$subject_id%>%group_by(subject_id)%>%

subject_id<-lung_tissue_matrix$subject_id
lung_tissue_matrix_value<-lung_tissue_matrix[,1:806882]
lung_tissue_matrix<-cbind(lung_tissue_matrix_value,subject_id)

lung_tissue_matrix<-lung_tissue_matrix%>%as.data.table()

gene_name<-lung_tissue_matrix_value%>%colnames()
unique_id<-unique(lung_tissue_matrix$subject_id)
length(lung_tissue_matrix$subject_id)


first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-lung_tissue_matrix_value[lung_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}


default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/lung_tissue.tsv",quote = F,row.names = FALSE)

#####################################################################################
library(reshape)
lung_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Lung")

lung_tissue_matrix<-lung_total_matrix%>%filter(sample_type=="Blood Derived Normal" )

subject_id<-lung_tissue_matrix$subject_id
lung_tissue_matrix_value<-lung_tissue_matrix[,1:806882]
lung_tissue_matrix<-cbind(lung_tissue_matrix_value,subject_id)

lung_tissue_matrix<-lung_tissue_matrix%>%as.data.table()
lung_tissue_matrix[,806883]%>%table()

gene_name<-lung_tissue_matrix_value%>%colnames()
unique_id<-unique(lung_tissue_matrix$subject_id)
first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-lung_tissue_matrix_value[lung_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}

default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/lung_blood.tsv",quote = F,row.names = FALSE)
####################################################################################








#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#################################head_tissue_dataframe######################################
head_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Head and Neck")
head_tissue_matrix<-head_total_matrix%>%filter(sample_type== "Primary Tumor"|sample_type=="Solid Tissue Normal")
#head_tissue_matrix$subject_id%>%group_by(subject_id)%>%

subject_id<-head_tissue_matrix$subject_id
head_tissue_matrix_value<-head_tissue_matrix[,1:806882]
head_tissue_matrix<-cbind(head_tissue_matrix_value,subject_id)

head_tissue_matrix<-head_tissue_matrix%>%as.data.table()

gene_name<-head_tissue_matrix_value%>%colnames()
unique_id<-unique(head_tissue_matrix$subject_id)
length(head_tissue_matrix$subject_id)


first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-head_tissue_matrix_value[head_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}


default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/head_tissue.tsv",quote = F,row.names = FALSE)

#####################################################################################
library(reshape)
head_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Head and Neck")

head_tissue_matrix<-head_total_matrix%>%filter(sample_type=="Blood Derived Normal" )

subject_id<-head_tissue_matrix$subject_id
head_tissue_matrix_value<-head_tissue_matrix[,1:806882]
head_tissue_matrix<-cbind(head_tissue_matrix_value,subject_id)

head_tissue_matrix<-head_tissue_matrix%>%as.data.table()
head_tissue_matrix[,806883]%>%table()

gene_name<-head_tissue_matrix_value%>%colnames()
unique_id<-unique(head_tissue_matrix$subject_id)
first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-head_tissue_matrix_value[head_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}

default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/head_blood.tsv",quote = F,row.names = FALSE)
####################################################################################








#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#################################stomach_tissue_dataframe######################################
stomach_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Stomach")
stomach_tissue_matrix<-stomach_total_matrix%>%filter(sample_type== "Primary Tumor"|sample_type=="Solid Tissue Normal")
#stomach_tissue_matrix$subject_id%>%group_by(subject_id)%>%

subject_id<-stomach_tissue_matrix$subject_id
stomach_tissue_matrix_value<-stomach_tissue_matrix[,1:806882]
stomach_tissue_matrix<-cbind(stomach_tissue_matrix_value,subject_id)

stomach_tissue_matrix<-stomach_tissue_matrix%>%as.data.table()

gene_name<-stomach_tissue_matrix_value%>%colnames()
unique_id<-unique(stomach_tissue_matrix$subject_id)
length(stomach_tissue_matrix$subject_id)


first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-stomach_tissue_matrix_value[stomach_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}


default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/stomach_tissue.tsv",quote = F,row.names = FALSE)

#####################################################################################
library(reshape)
stomach_total_matrix<-last_abundance_matrix%>%filter(primary_site=="Stomach")

stomach_tissue_matrix<-stomach_total_matrix%>%filter(sample_type=="Blood Derived Normal" )

subject_id<-stomach_tissue_matrix$subject_id
stomach_tissue_matrix_value<-stomach_tissue_matrix[,1:806882]
stomach_tissue_matrix<-cbind(stomach_tissue_matrix_value,subject_id)

stomach_tissue_matrix<-stomach_tissue_matrix%>%as.data.table()
stomach_tissue_matrix[,806883]%>%table()

gene_name<-stomach_tissue_matrix_value%>%colnames()
unique_id<-unique(stomach_tissue_matrix$subject_id)
first_num=0
default_df<-data.frame(gene_name=gene_name)
for (i in unique_id){
  first_num=first_num+1
  print(first_num)
  print(i)
  new_df<-stomach_tissue_matrix_value[stomach_tissue_matrix$subject_id==i,]%>%colSums()%>%as.data.frame()
  default_df<-cbind(default_df,new_df)
}

default_df<-column_to_rownames(default_df,loc=1)
colnames(default_df)<-unique_id
fwrite(default_df,file="/home/junwoojo/다운로드/decomtam_matrix/stomach_blood.tsv",quote = F,row.names = FALSE)



###########################################################################################
