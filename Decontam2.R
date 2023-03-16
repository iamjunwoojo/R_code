
###################################################################################
###################################################################################
#########################Prevalence table_colon####################################
###################################################################################
###################################################################################
library(data.table)
colon_tissue_mat<-fread(file="~/다운로드/decomtam_matrix/colon_tissue.tsv")%>%as.data.frame()
colon_blood_mat<-fread(file="~/다운로드/decomtam_matrix/colon_blood.tsv")%>%as.data.frame()
colon_matched_id<-colnames(colon_tissue_mat)[colnames(colon_tissue_mat)%in%colnames(colon_blood_mat)]

############################################Tissue
colon_tissue_mat
colon_tissue_prevalence <-apply(colon_tissue_mat[,colon_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
colon_tissue_prevalence <- apply(colon_tissue_mat[,colon_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
colon_tissue_prevalence<-t(colon_tissue_prevalence%>%as.data.frame())

rownames(colon_tissue_prevalence)<-c(1:806882)
colnames(colon_tissue_prevalence)<-c('Exist','NotExist')
colon_tissue_prevalence<-colon_tissue_prevalence%>%as.data.frame()%>%mutate(Type='Tissue',gene_name=c(1:806882))

###########################################Blood
colon_blood_mat
colon_blood_prevalence <-apply(colon_blood_mat[,colon_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
colon_blood_prevalence <- apply(colon_blood_mat[,colon_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
colon_blood_prevalence<-t(colon_blood_prevalence%>%as.data.frame())

rownames(colon_blood_prevalence)<-c(1:806882)
colnames(colon_blood_prevalence)<-c('Exist','NotExist')
colon_blood_prevalence<-colon_blood_prevalence%>%as.data.frame()%>%mutate(Type='Blood',gene_name=c(1:806882))

###############################################
total_colon_prevalence<-rbind(colon_tissue_prevalence,colon_blood_prevalence)
total_colon_prevalence<-total_colon_prevalence%>%group_by(gene_name)%>%mutate(sumdata=sum(Exist))%>%as.data.frame()
total_colon_prevalence<-total_colon_prevalence%>%filter(sumdata!=0)%>%select(c(1:4))


###################################################################################
###################################################################################
#########################Prevalence Esophagus####################################
###################################################################################
###################################################################################

Esophagus_tissue_mat<-fread(file="~/다운로드/decomtam_matrix/Esophagus_tissue.tsv")%>%as.data.frame()
Esophagus_blood_mat<-fread(file="~/다운로드/decomtam_matrix/Esophagus_blood.tsv")%>%as.data.frame()



Esophagus_matched_id<-colnames(Esophagus_tissue_mat)[colnames(Esophagus_tissue_mat)%in%colnames(Esophagus_blood_mat)]

############################################Tissue

Esophagus_tissue_prevalence <-apply(Esophagus_tissue_mat[,Esophagus_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
Esophagus_tissue_prevalence <- apply(Esophagus_tissue_mat[,Esophagus_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
Esophagus_tissue_prevalence<-t(Esophagus_tissue_prevalence%>%as.data.frame())

rownames(Esophagus_tissue_prevalence)<-c(1:806882)
colnames(Esophagus_tissue_prevalence)<-c('Exist','NotExist')
Esophagus_tissue_prevalence<-Esophagus_tissue_prevalence%>%as.data.frame()%>%mutate(Type='Tissue',gene_name=c(1:806882))

###########################################Blood

Esophagus_blood_prevalence <-apply(Esophagus_blood_mat[,Esophagus_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
Esophagus_blood_prevalence <- apply(Esophagus_blood_mat[,Esophagus_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
Esophagus_blood_prevalence<-t(Esophagus_blood_prevalence%>%as.data.frame())

rownames(Esophagus_blood_prevalence)<-c(1:806882)
colnames(Esophagus_blood_prevalence)<-c('Exist','NotExist')
Esophagus_blood_prevalence<-Esophagus_blood_prevalence%>%as.data.frame()%>%mutate(Type='Blood',gene_name=c(1:806882))

###############################################
total_Esophagus_prevalence<-rbind(Esophagus_tissue_prevalence,Esophagus_blood_prevalence)
total_Esophagus_prevalence<-total_Esophagus_prevalence%>%group_by(gene_name)%>%mutate(sumdata=sum(Exist))%>%as.data.frame()
total_Esophagus_prevalence<-total_Esophagus_prevalence%>%filter(sumdata!=0)%>%select(c(1:4))

 
###################################################################################
###################################################################################
#########################Prevalence Lung####################################
###################################################################################
###################################################################################

lung_tissue_mat<-fread(file="~/다운로드/decomtam_matrix/lung_tissue.tsv")%>%as.data.frame()
lung_blood_mat<-fread(file="~/다운로드/decomtam_matrix/lung_blood.tsv")%>%as.data.frame()

lung_matched_id<-colnames(lung_tissue_mat)[colnames(lung_tissue_mat)%in%colnames(lung_blood_mat)]

############################################Tissue

lung_tissue_prevalence <-apply(lung_tissue_mat[,lung_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
lung_tissue_prevalence <- apply(lung_tissue_mat[,lung_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
lung_tissue_prevalence<-t(lung_tissue_prevalence%>%as.data.frame())

rownames(lung_tissue_prevalence)<-c(1:806882)
colnames(lung_tissue_prevalence)<-c('Exist','NotExist')
lung_tissue_prevalence<-lung_tissue_prevalence%>%as.data.frame()%>%mutate(Type='Tissue',gene_name=c(1:806882))

###########################################Blood

lung_blood_prevalence <-apply(lung_blood_mat[,lung_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
lung_blood_prevalence <- apply(lung_blood_mat[,lung_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
lung_blood_prevalence<-t(lung_blood_prevalence%>%as.data.frame())

rownames(lung_blood_prevalence)<-c(1:806882)
colnames(lung_blood_prevalence)<-c('Exist','NotExist')
lung_blood_prevalence<-lung_blood_prevalence%>%as.data.frame()%>%mutate(Type='Blood',gene_name=c(1:806882))

###############################################
total_lung_prevalence<-rbind(lung_tissue_prevalence,lung_blood_prevalence)
total_lung_prevalence<-total_lung_prevalence%>%group_by(gene_name)%>%mutate(sumdata=sum(Exist))%>%as.data.frame()
total_lung_prevalence<-total_lung_prevalence%>%filter(sumdata!=0)%>%select(c(1:4))

###################################################################################
###################################################################################
#########################Prevalence Stomach####################################
###################################################################################
###################################################################################
stomach_tissue_mat<-fread(file="~/다운로드/decomtam_matrix/Stomach_tissue.tsv")%>%as.data.frame()
stomach_blood_mat<-fread(file="~/다운로드/decomtam_matrix/Stomach_blood.tsv")%>%as.data.frame()

stomach_matched_id<-colnames(stomach_tissue_mat)[colnames(stomach_tissue_mat)%in%colnames(stomach_blood_mat)]

############################################Tissue

stomach_tissue_prevalence <-apply(stomach_tissue_mat[,stomach_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
stomach_tissue_prevalence <- apply(stomach_tissue_mat[,stomach_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
stomach_tissue_prevalence<-t(stomach_tissue_prevalence%>%as.data.frame())

rownames(stomach_tissue_prevalence)<-c(1:806882)
colnames(stomach_tissue_prevalence)<-c('Exist','NotExist')
stomach_tissue_prevalence<-stomach_tissue_prevalence%>%as.data.frame()%>%mutate(Type='Tissue',gene_name=c(1:806882))

###########################################Blood

stomach_blood_prevalence <-apply(stomach_blood_mat[,stomach_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
stomach_blood_prevalence <- apply(stomach_blood_mat[,stomach_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
stomach_blood_prevalence<-t(stomach_blood_prevalence%>%as.data.frame())

rownames(stomach_blood_prevalence)<-c(1:806882)
colnames(stomach_blood_prevalence)<-c('Exist','NotExist')
stomach_blood_prevalence<-stomach_blood_prevalence%>%as.data.frame()%>%mutate(Type='Blood',gene_name=c(1:806882))

###############################################
total_stomach_prevalence<-rbind(stomach_tissue_prevalence,stomach_blood_prevalence)
total_stomach_prevalence<-total_stomach_prevalence%>%group_by(gene_name)%>%mutate(sumdata=sum(Exist))%>%as.data.frame()
total_stomach_prevalence<-total_stomach_prevalence%>%filter(sumdata!=0)%>%select(c(1:4))


###################################################################################
###################################################################################
#########################Prevalence Head & Neck ####################################
###################################################################################
###################################################################################

head_tissue_mat<-fread(file="~/다운로드/decomtam_matrix/head_tissue.tsv")%>%as.data.frame()
head_blood_mat<-fread(file="~/다운로드/decomtam_matrix/head_blood.tsv")%>%as.data.frame()

head_matched_id<-colnames(head_tissue_mat)[colnames(head_tissue_mat)%in%colnames(head_blood_mat)]

############################################Tissue

head_tissue_prevalence <-apply(head_tissue_mat[,head_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
head_tissue_prevalence <- apply(head_tissue_mat[,head_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
head_tissue_prevalence<-t(head_tissue_prevalence%>%as.data.frame())

rownames(head_tissue_prevalence)<-c(1:806882)
colnames(head_tissue_prevalence)<-c('Exist','NotExist')
head_tissue_prevalence<-head_tissue_prevalence%>%as.data.frame()%>%mutate(Type='Tissue',gene_name=c(1:806882))

###########################################Blood

head_blood_prevalence <-apply(head_blood_mat[,head_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
head_blood_prevalence <- apply(head_blood_mat[,head_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
head_blood_prevalence<-t(head_blood_prevalence%>%as.data.frame())

rownames(head_blood_prevalence)<-c(1:806882)
colnames(head_blood_prevalence)<-c('Exist','NotExist')
head_blood_prevalence<-head_blood_prevalence%>%as.data.frame()%>%mutate(Type='Blood',gene_name=c(1:806882))

###############################################
total_head_prevalence<-rbind(head_tissue_prevalence,head_blood_prevalence)
total_head_prevalence<-total_head_prevalence%>%group_by(gene_name)%>%mutate(sumdata=sum(Exist))%>%as.data.frame()
total_head_prevalence<-total_head_prevalence%>%filter(sumdata!=0)%>%select(c(1:4))



###################################################################################
###################################################################################
#########################Prevalence liver #########################################
###################################################################################
###################################################################################

liver_tissue_mat<-fread(file="~/다운로드/decomtam_matrix/liver_tissue.tsv")%>%as.data.frame()
liver_blood_mat<-fread(file="~/다운로드/decomtam_matrix/liver_blood.tsv")%>%as.data.frame()

liver_matched_id<-colnames(liver_tissue_mat)[colnames(liver_tissue_mat)%in%colnames(liver_blood_mat)]

############################################Tissue

liver_tissue_prevalence <-apply(liver_tissue_mat[,liver_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
liver_tissue_prevalence <- apply(liver_tissue_mat[,liver_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
liver_tissue_prevalence<-t(liver_tissue_prevalence%>%as.data.frame())

rownames(liver_tissue_prevalence)<-c(1:806882)
colnames(liver_tissue_prevalence)<-c('Exist','NotExist')
liver_tissue_prevalence<-liver_tissue_prevalence%>%as.data.frame()%>%mutate(Type='Tissue',gene_name=c(1:806882))

###########################################Blood

liver_blood_prevalence <-apply(liver_blood_mat[,liver_matched_id],1,function(x){sum(x>0)/length(x)})%>%as.data.frame()
liver_blood_prevalence <- apply(liver_blood_mat[,liver_matched_id],1,function(x){c(sum(x>0),length(x)-sum(x>0))})%>%as.data.frame()
liver_blood_prevalence<-t(liver_blood_prevalence%>%as.data.frame())

rownames(liver_blood_prevalence)<-c(1:806882)
colnames(liver_blood_prevalence)<-c('Exist','NotExist')
liver_blood_prevalence<-liver_blood_prevalence%>%as.data.frame()%>%mutate(Type='Blood',gene_name=c(1:806882))

###############################################
total_liver_prevalence<-rbind(liver_tissue_prevalence,liver_blood_prevalence)
total_liver_prevalence<-total_liver_prevalence%>%group_by(gene_name)%>%mutate(sumdata=sum(Exist))%>%as.data.frame()
total_liver_prevalence<-total_liver_prevalence%>%filter(sumdata!=0)%>%select(c(1:4))


######################################################################################
######################################################################################
######################################################################################
###########################   Fisher Exact   #########################################
######################################################################################
######################################################################################
######################################################################################

###################################colon################################################
length(total_colon_prevalence$gene_name%>%unique())
ttt<-c()

k=0
for (i in total_colon_prevalence$gene_name%>%unique()){
  test<-total_colon_prevalence%>%filter(gene_name==i)
  rownames(test)<-test$Type
  test<-test[,1:2]
  fisher_test<-fisher.test(test, alternative = "greater")
  k=k+1
  print(k)
  ttt<-c(ttt,fisher_test$p.value)
}

Tissue<-total_colon_prevalence%>%filter(Type=='Tissue')
Blood<-total_colon_prevalence%>%filter(Type=='Blood')

df_preval<-data.frame(Tissue_prevalence=Tissue[,1], Blood_prevalence=Blood[,1],gene_name=total_colon_prevalence$gene_name%>%unique(),pval=ttt)

ggg<-df_preval%>%filter(Tissue_prevalence>Blood_prevalence & pval<0.05)
colon_sample_number=sum(total_colon_prevalence[1,1:2])
ggg<-ggg%>%filter( !(Blood_prevalence/colon_sample_number)>0.2 )
colon_tissue_enriched_gene<-ggg$gene_name


######################################liver############################################################
length(total_liver_prevalence$gene_name%>%unique())
ttt<-c()

k=0
for (i in total_liver_prevalence$gene_name%>%unique()){
  test<-total_liver_prevalence%>%filter(gene_name==i)
  rownames(test)<-test$Type
  test<-test[,1:2]
  fisher_test<-fisher.test(test, alternative = "greater")
  k=k+1
  print(k)
  ttt<-c(ttt,fisher_test$p.value)
}



Tissue<-total_liver_prevalence%>%filter(Type=='Tissue')
Blood<-total_liver_prevalence%>%filter(Type=='Blood')

df_preval<-data.frame(Tissue_prevalence=Tissue[,1], Blood_prevalence=Blood[,1],gene_name=total_liver_prevalence$gene_name%>%unique(),pval=ttt)

ggg<-df_preval%>%filter(Tissue_prevalence>Blood_prevalence & pval<0.05)
liver_sample_number=sum(total_liver_prevalence[1,1:2])
ggg<-ggg%>%filter( !(Blood_prevalence/liver_sample_number)>0.2 )
liver_tissue_enriched_gene<-ggg$gene_name



##########################head & neck ###############################################


length(total_head_prevalence$gene_name%>%unique())
ttt<-c()

k=0
for (i in total_head_prevalence$gene_name%>%unique()){
  test<-total_head_prevalence%>%filter(gene_name==i)
  rownames(test)<-test$Type
  test<-test[,1:2]
  fisher_test<-fisher.test(test, alternative = "greater")
  k=k+1
  print(k)
  ttt<-c(ttt,fisher_test$p.value)
}




Tissue<-total_head_prevalence%>%filter(Type=='Tissue')
Blood<-total_head_prevalence%>%filter(Type=='Blood')

df_preval<-data.frame(Tissue_prevalence=Tissue[,1], Blood_prevalence=Blood[,1],gene_name=total_head_prevalence$gene_name%>%unique(),pval=ttt)

ggg<-df_preval%>%filter(Tissue_prevalence>Blood_prevalence & pval<0.05)
head_sample_number=sum(total_head_prevalence[1,1:2])
ggg<-ggg%>%filter( !(Blood_prevalence/head_sample_number)>0.2 )
head_tissue_enriched_gene<-ggg$gene_name


##############################################Esophagus ###############################################################################lung liver colon stomach


length(total_Esophagus_prevalence$gene_name%>%unique())
ttt<-c()

k=0
for (i in total_Esophagus_prevalence$gene_name%>%unique()){
  test<-total_Esophagus_prevalence%>%filter(gene_name==i)
  rownames(test)<-test$Type
  test<-test[,1:2]
  fisher_test<-fisher.test(test, alternative = "greater")
  k=k+1
  print(k)
  ttt<-c(ttt,fisher_test$p.value)
}


Tissue<-total_Esophagus_prevalence%>%filter(Type=='Tissue')
Blood<-total_Esophagus_prevalence%>%filter(Type=='Blood')

df_preval<-data.frame(Tissue_prevalence=Tissue[,1], Blood_prevalence=Blood[,1],gene_name=total_Esophagus_prevalence$gene_name%>%unique(),pval=ttt)

ggg<-df_preval%>%filter(Tissue_prevalence>Blood_prevalence & pval<0.05)
Esophagus_sample_number=sum(total_Esophagus_prevalence[1,1:2])
ggg<-ggg%>%filter( !(Blood_prevalence/Esophagus_sample_number)>0.2 )
Esophagus_tissue_enriched_gene<-ggg$gene_name
length(Esophagus_tissue_enriched_gene)


##############################################Stomach ###############################################################################lung liver colon stomach


length(total_stomach_prevalence$gene_name%>%unique())
ttt<-c()

k=0
for (i in total_stomach_prevalence$gene_name%>%unique()){
  test<-total_stomach_prevalence%>%filter(gene_name==i)
  rownames(test)<-test$Type
  test<-test[,1:2]
  fisher_test<-fisher.test(test, alternative = "greater")
  k=k+1
  print(k)
  ttt<-c(ttt,fisher_test$p.value)
}


Tissue<-total_stomach_prevalence%>%filter(Type=='Tissue')
Blood<-total_stomach_prevalence%>%filter(Type=='Blood')

df_preval<-data.frame(Tissue_prevalence=Tissue[,1], Blood_prevalence=Blood[,1],gene_name=total_stomach_prevalence$gene_name%>%unique(),pval=ttt)

ggg<-df_preval%>%filter(Tissue_prevalence>Blood_prevalence & pval<0.05)
stomach_sample_number=sum(total_stomach_prevalence[1,1:2])
ggg<-ggg%>%filter( !(Blood_prevalence/stomach_sample_number)>0.2 )
stomach_tissue_enriched_gene<-ggg$gene_name
length(stomach_tissue_enriched_gene)



##############################################lung ###############################################################################lung liver colon stomach


length(total_lung_prevalence$gene_name%>%unique())
ttt<-c()

k=0
for (i in total_lung_prevalence$gene_name%>%unique()){
  test<-total_lung_prevalence%>%filter(gene_name==i)
  rownames(test)<-test$Type
  test<-test[,1:2]
  fisher_test<-fisher.test(test, alternative = "greater")
  k=k+1
  print(k)
  ttt<-c(ttt,fisher_test$p.value)
}


Tissue<-total_lung_prevalence%>%filter(Type=='Tissue')
Blood<-total_lung_prevalence%>%filter(Type=='Blood')

df_preval<-data.frame(Tissue_prevalence=Tissue[,1], Blood_prevalence=Blood[,1],gene_name=total_lung_prevalence$gene_name%>%unique(),pval=ttt)

ggg<-df_preval%>%filter(Tissue_prevalence>Blood_prevalence & pval<0.05)
lung_sample_number=sum(total_lung_prevalence[1,1:2])
ggg<-ggg%>%filter( !(Blood_prevalence/lung_sample_number)>0.2 )
lung_tissue_enriched_gene<-ggg$gene_name
length(lung_tissue_enriched_gene)


########################################################################################################

















enriched_all_gene<-unique(c(Esophagus_tissue_enriched_gene,colon_tissue_enriched_gene,
    liver_tissue_enriched_gene,lung_tissue_enriched_gene,
    stomach_tissue_enriched_gene,head_tissue_enriched_gene))

for (i in c('stomach_tissue_enriched_gene','colon_tissue_enriched_gene',
            'liver_tissue_enriched_gene','lung_tissue_enriched_gene',
            'stomach_tissue_enriched_gene','head_tissue_enriched_gene')){
print(paste(i,length(get(i))))

  }


length(enriched_all_gene)
