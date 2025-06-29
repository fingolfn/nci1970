####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##       ####      ##     #########   ##               ##   #######   ##########   ########             #########     ########                                               
##       ## ##     ##    ##           ##              ###  ##     ##          ##  ##      ##            ##      ##    ##     ##                    
##       ##  ##    ##   ##            ##               ##  ##     ##         ##   ##      ##            ##       ##   ##     ##                     
##       ##   ##   ##   ##            ##               ##   ########        ##    ##      ##            ##       ##   ########                                   
##       ##    ##  ##   ##            ##               ##         ##       ##     ##      ##            ##       ##   ##     ##                   
##       ##     ## ##    ##           ##               ##         ##      ##      ##      ##            ##      ##    ##     ##                
##       ##      ####     #########   ##               ##         ##     ##        ########             #########     ########                                
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
#DATE:   01-13-2022
#GOAL:   Cleaned up Code for Main Figures & Supporting Information Figures organized in to parallel manuscript structure
#STRUCTURE:
#       FIGURE #1:      Comparison of clinical (NCI1970meta), in vitro chemo-sensitivity data with FDA indications
#       FIGURE #2:      Biomarker definitions and correlation with clinical (an invitro) chemo-sensitivity
#       FIGURE #S7-9    Multi-linear regression of expr-biomarkers demonstrated independent correlations with drug-resposne
#       FIGURE #S10     Literature experimental data on anthracycline intracellular concentration and IC50 changes (to validate chemo-resistance model)
#       FIGURE #3:      2D mechanistic interpretation of drug sensitivity:   total metabolism (x-axis) vs total efflux (y-axis)
#       FIGURE #4:      Clincial correlation of Oncotype/Pam50/Mammaprint signatures (estrogen & proliferation) with chemo-resistance signatures
#       FIGURE #5:      Estimate total chemo-transport and metabolism at human barrier tissues and in immune cells using scRNAseq data
#
#################################################################################################################################################################################################################################################################### 





####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##            ###                                                          
##     ##          ##    ##                          ##   ##           ####                                      
##     ##          ##   ##                      #################        ##                                            
##     ########    ##   ##   ####                    ##   ##             ##                                             
##     ##          ##   ##      ##              #################        ##                                             
##     ##          ##    ##     ##                  ##   ##              ##                                    
##     ##          ##     #######                   ##   ##           ########                                                     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       ####      ##     #########   ##               ##   #######   ##########   ########             #########     ########                                               
##       ## ##     ##    ##           ##              ###  ##     ##          ##  ##      ##            ##      ##    ##     ##                    
##       ##  ##    ##   ##            ##               ##  ##     ##         ##   ##      ##            ##       ##   ##     ##                     
##       ##   ##   ##   ##            ##               ##   ########        ##    ##      ##            ##       ##   ########                                   
##       ##    ##  ##   ##            ##               ##         ##       ##     ##      ##            ##       ##   ##     ##                   
##       ##     ## ##    ##           ##               ##         ##      ##      ##      ##            ##      ##    ##     ##                
##       ##      ####     #########   ##               ##         ##     ##        ########             #########     ########                                
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##                         #########  #######        ##                                ########   ########## ########      #######                          
##                         ##         ##    ##      ####                              ##              ##     ##     ##     ##    ##      
##                         ##         ##     ##    ##  ##                            ##               ##     ##     ##     ##    ##                         
##          ##    ##       ########   ##     ##   ##    ##            ##    ##       ##               ##     ########      #######                                                  
##           ##  ##        ##         ##     ##  ##########            ##  ##        ##               ##     ##    ##      ##                                                    
##            ####         ##         ##    ##   ##      ##             ####          ##              ##     ##     ##     ##      
##             ##          ##         #######    ##      ##              ##            ########       ##     ##      ##    ##                                                             
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################




####################################################################################################################################################################################################################################################################  
#FIGURE S3A-B:   Data distribution (19 cancers x 34 drugs) for NCI-textbook (A) and textbook + 63 metareviews (B)
####################################################################################################################################################################################################################################################################  

#NEEDED LIBRARIES:
library("gplots")
library("RColorBrewer")

#PLOTTING FUNCTION:
hcluster_unsortBW<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=NA,Colv=NA, cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("white", "black"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab = units,symm=F,dendrogram="none")
        
        return(rowMeans(matrix))
}


#LOAD DATA
load("early_db/NCI1970-monotherapy-DB_simple.RData")
load("early_db/NCI1970-META_monotherapyDB.RData")

##########################################
#ORGANIZE DATA BY ROW & COL SUMS
##########################################
#sort drugs & cancers for BW-heatmap
sorted_cancers<-names(sort(colSums(NCI1970_META_numbers),decreasing=T))
sorted_drugs<-names(sort(rowSums(NCI1970_META_numbers),decreasing=T))
drugs_na_num<-sort(rowSums(is.na(data.frame(NCI1970_META_ORrate))))
good_drugs<-names(drugs_na_num)[1:34]

##########################################
#FINAL FORMATTING TO SHOW GAPS:
##########################################
#ADD BLANK COLUMNS TO NCI1970
new_drugs<-good_drugs[-which(good_drugs %in% rownames(NCI1970_simp_numbers))]
new_blank_matrix<-matrix(0,ncol=ncol(NCI1970_simp_numbers),nrow=length(new_drugs))
colnames(new_blank_matrix)<-colnames(NCI1970_simp_numbers)
rownames(new_blank_matrix)<-new_drugs

#FINAL DATA-matrices for paper:
NCI1970_orig_number_FINAL<-rbind(NCI1970_simp_numbers,new_blank_matrix)[good_drugs,sorted_cancers]
NCI1970_update_number_FINAL<-NCI1970_META_numbers[good_drugs,sorted_cancers]
NCI1970_update_ORnum_FINAL<-NCI1970_META_ORnum[good_drugs,sorted_cancers]
NCI1970_update_ORrate_FINAL<-NCI1970_META_ORrate[good_drugs,sorted_cancers]


##########################################
#STATISTICS OF EACH MATRIX:
##########################################

#ORIGINAL Patients:  19,958
sum(NCI1970_simp_numbers)

#ORIGINAL NA's:  58.4%   =   211/361 
orig_nas<-sum(is.na(data.frame(NCI1970_simp_ORrate)))
orig_cells<-ncol(NCI1970_simp_ORrate)*nrow(NCI1970_simp_ORrate)# 19 drugs x 19 cancers

#FINAL Patients:  49,002
sum(NCI1970_update_number_FINAL)

#FINAL NA's: 41.7%   =   270/646  
final_nas<-sum(is.na(data.frame(NCI1970_META_ORrate[good_drugs,])))
final_cells<-ncol(NCI1970_META_ORrate[good_drugs,])*nrow(NCI1970_META_ORrate[good_drugs,])# 34 drugs x 19 cancers


##########################################
#FINAL PLOTS:
##########################################

#Figure S3A
par(las=1)
hcluster_unsortBW(t(log10(NCI1970_orig_number_FINAL+1)),units="log10(patient #)",title="NCI1970 Textbook alone",
                  col_size = 0.7,row_size = 1.2)



#Figure S3B
par(las=1)
hcluster_unsortBW(t(log10(NCI1970_update_number_FINAL+1)),units="log10(patient #)",title="NCI1970 + Meta Reviews",
                  col_size = 0.7,row_size = 1.2)






####################################################################################################################################################################################################################################################################  
#FIGURE S3C:   Total # of patients per cancer-type
####################################################################################################################################################################################################################################################################  

cancer_num<-sort(colSums(NCI1970_update_number_FINAL),decreasing=T)
#brca nsclc  coad lymph    ov  leuk  skcm  hnsc  sarc  cesc  paad  tgct  kirc  stad  blca  esca  ucec   gbm  lihc 
#7070  6136  5555  5230  4982  4668  3004  2892  1422  1093  1065  1061  1056  1030   856   635   425   418   404 

#Figure 1C
par(las=2,mar=c(5, 4, 4, 2)) # make label text perpendicular to axis
barplot(sort(colSums(NCI1970_update_number_FINAL),decreasing=F), main="NCI1970:  # Patients / Cancer", 
        beside=F,cex.axis = 1.5,cex=1.4,cex.main=2,cex.names=1,ylab=c(0,100),col="black",horiz=T,log="x")




####################################################################################################################################################################################################################################################################  
#FIGURE S3D:  Correlation Total data with Cancer Incidence in USA (according to NCI SEER report for 2017)  
####################################################################################################################################################################################################################################################################  

load("clinical_dbs/SEER_cancer_incidence.RData")

par(las=1)
linear_fit_plot<-function(x,y,pt_colors="black"){
        #linear regression:
        lin_reg<-lm(y~x)
        #make fit equation + r-squared:
        pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
        fit_coeff<-round(coef(lin_reg),6)
        r2 <- round(summary(lin_reg)$r.squared, 2)
        rmse <- round(sqrt(mean(resid(lin_reg)^2)), 2)
        eq_r2<-paste("y = ", fit_coeff[2], "x + ", fit_coeff[1]," "," "," ","r^2 = ", r2," "," "," ","rmse = ",rmse)
        #plot data:
        plot(x,y, pch = 16, cex = 1.3, col = pt_colors)
        #add fit-line and fit-equation:
        abline(lin_reg)
        mtext(eq_r2, 3, line=-2,cex=1)
}


seer_vector<-SEER_cancer_incidence[names(cancer_num)]

################
#PLOT:
################
#BEST FIT LINE:  y = 96.44x + 861.98       r^2 = 0.58    rmse = 1436.51
linear_fit_plot(seer_vector,cancer_num)

# Figure 1D
#ACTUAL PLOT:
plot(seer_vector,cancer_num,pch=16,xlim=c(0,80),ylim=c(0,8000),xlab="Cancer Incidence (Millions)",ylab="# patients")
text(seer_vector+2.5,cancer_num,labels=names(cancer_num))
abline(a=0,b=121.18)




######################################################################################################################################################################
#FIGURE 1A-B:  Heatmap of ORR for 34 cancers x 19 cancers (A) and average drug-sensitivity of 19 cancers (B)
######################################################################################################################################################################


load("NCI1970meta_Database.RData")

load("clinical_dbs/FDA_indictations.RData")

load("cell-line_dbs/invitro_AUC-lineage-AVG_dbs.RData")

hcluster_unsort<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=NA,Colv=NA, cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "lightcyan","orangered","red","darkred"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab = units,symm=F,dendrogram='none')
        
        return(rowMeans(matrix,na.rm=T))
}
hcluster_unsortFDA<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=NA,Colv=NA, cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "lightcyan","red"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab = units,symm=F,dendrogram='none')
        
        return(rowMeans(matrix,na.rm=T))
}
hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",key="z-score"){
        #distance matrix:
        dist_columns <- dist(t(matrix))
        dist_rows <- dist(matrix)
        
        #hclustering:
        hclust_columns<-hclust(dist_columns, method="average")
        hclust_rows<-hclust(dist_rows, method="average")
        
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), 
                  cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "lightcyan","orangered","red","darkred"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = key,key.xlab = key)
}
hcluster_corr<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",key="z-score"){
        #distance matrix:
        cor_columns <- cor(matrix,method="spearman")
        cor_rows <- cor(t(matrix),method="spearman")
        
        
        #distance matrix:
        dist_columns <- as.dist(1-cor_columns)
        dist_rows <- as.dist(1-cor_rows)
        
        #hclustering:
        hclust_columns<-hclust(dist_columns, method="average")
        hclust_rows<-hclust(dist_rows, method="average")
        
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), 
                  cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "lightcyan","orangered","red","darkred"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = key,key.xlab = key)
}



#################################
#PARSE OVERLAP
#################################

#FILTER CANCER BY POST-DRUG-FILTERED NA's
cancers_na_num<-sort(colSums(is.na(data.frame(NCI1970meta_rate))))
good_cancers<-names(cancers_na_num)[1:18]

filtered_ORR<-NCI1970meta_rate[,good_cancers]

orr_sorted_cancers<-names(sort(colMeans(filtered_ORR,na.rm=T),decreasing=T))
orr_sorted_drugs<-names(sort(rowMeans(filtered_ORR,na.rm=T),decreasing=T))
drug_overlap<-orr_sorted_drugs[which(orr_sorted_drugs %in% rownames(FDA_indictations))]


#Figure 2A
hcluster_unsort(NCI1970meta_rate[drug_overlap,orr_sorted_cancers],title="NCI1970meta Response Rates",units="% Response",col_size = 1.3,row_size = 0.7)


#Figure 2B
#AVG CANCER-TYPE ORR Barchart:
par(las=2,mar=c(5, 4, 4, 2)) # make label text perpendicular to axis
barplot(sort(colMeans(NCI1970meta_rate[,orr_sorted_cancers],na.rm=T),decreasing=F), main="AVG chemo response rate", beside=F,cex.axis = 1,cex=1,cex.main=2,ylab=c(0,100),col="black",horiz=T)






######################################################################################################################################################################
#FIGURE 1C-D:  Heatmap of FDA-indications (A) correlation of avg-ORR with total # FDA-indications (B) 
######################################################################################################################################################################

linear_fit_plot<-function(x,y,pt_colors="black"){
        #linear regression:
        lin_reg<-lm(y~x)
        #make fit equation + r-squared:
        pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
        fit_coeff<-round(coef(lin_reg),6)
        r2 <- round(summary(lin_reg)$r.squared, 2)
        rmse <- round(sqrt(mean(resid(lin_reg)^2)), 2)
        eq_r2<-paste("y = ", fit_coeff[2], "x + ", fit_coeff[1]," "," "," ","r^2 = ", r2," "," "," ","rmse = ",rmse)
        #plot data:
        plot(x,y, pch = 16, cex = 1.3, col = pt_colors)
        #add fit-line and fit-equation:
        abline(lin_reg)
        mtext(eq_r2, 3, line=-2,cex=1)
}


# Figure 2C

hcluster_unsortFDA(FDA_indictations[drug_overlap,orr_sorted_cancers],title="FDA Indications",units="binary",col_size = 1.3,row_size = 0.7)


####################################
#MULTIPLE LINEAR REGRESSION:   FDA = F( frequency, efficacy)
#####################################


mean_ORR<-colMeans(filtered_ORR[,orr_sorted_cancers],na.rm=T)
sum_FDAindic<-colSums(FDA_indictations[,orr_sorted_cancers],na.rm=T)
sum_NCIpatients<-colSums(NCI1970meta_total[,orr_sorted_cancers])

############################
#FIT:  Efficacy r2 = 0.41
############################
par(las=1)
linear_fit_plot(mean_ORR,sum_FDAindic)

############################
#FIT:  patient-data r2 = 0.59
############################
linear_fit_plot(sum_NCIpatients,sum_FDAindic)

linear_fit_plot(sum_NCIpatients,mean_ORR)



############################
#MULTI-FIT:  patient-data r2 = 0.78
############################

dataframe_multiLinReg<-cbind("FDA"=sum_FDAindic,"patients"=sum_NCIpatients,"efficacy"=mean_ORR)

FDAindications_model<-lm(FDA~patients+efficacy,data=data.frame(dataframe_multiLinReg))

summary(FDAindications_model)
#Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -3.873041   2.509624  -1.543 0.143595    
#patients     0.002384   0.000473   5.040 0.000147 ***
#efficacy    38.830404  10.844881   3.581 0.002733 ** 


#Multiple R-squared:  0.7807,	Adjusted R-squared:  0.7514 

glm_model<-(mean_ORR*38.830404+sum_NCIpatients*0.002384)-3.87

linear_fit_plot(glm_model,sum_FDAindic)

# Figure 2D

#ACTUAL PLOT:
plot(glm_model,sum_FDAindic,pch=16,xlab="Linear Model (ORR + total data)",ylab="# indicationed drugs")
text(glm_model+0.8,sum_FDAindic,labels=names(sum_FDAindic))
abline(a=0,b=1)



######################################################################################################################################################################
#FIGURE 1E-F & S4 & S5:  Heatmap of CTRPv2 (and GDSC) drug-sensitivity (A) and correlation of clinical & in vitro drug-sensitivity (B)
######################################################################################################################################################################

########################################################
#Figure S4A,C & S5A-B:  CTRP ORGANIZE & NORMALIZE DATA:
########################################################

rownames(ctrp_avg_matrix)<-tolower(rownames(ctrp_avg_matrix))

ctrp_cancer_overlap<-orr_sorted_cancers[which(orr_sorted_cancers %in% colnames(ctrp_avg_matrix))]
ctrp_drug_overlap<-tolower(orr_sorted_drugs[which(tolower(orr_sorted_drugs) %in% rownames(ctrp_avg_matrix))])

norm_nci1970_overlap<-(16-ctrp_avg_matrix[ctrp_drug_overlap,ctrp_cancer_overlap])/16
norm_nci1970_overlap_scaled<-t(scale(t(norm_nci1970_overlap),scale=T,center=F))#Divide by Root-mean-square = normalize scale of differential across lineages for each drug

############################
#FINAL FIGURE S4A,C: CTRP CLUSTERED:
############################
hcluster((16-ctrp_avg_matrix[ctrp_drug_overlap,])/16)

hcluster(t(scale(t((16-ctrp_avg_matrix[ctrp_drug_overlap,])/16),scale=T,center=F)))


############################
#Figure 1E & S5A
############################
hcluster_unsort(norm_nci1970_overlap_scaled,title="CTRP:  lineage-average",units="(16-AUC)/16   / drug-rt-mean-sq",col_size = 1.3,row_size = 0.9)


########################################################
#Figure S4B,D & S5C-D:  GDSC  17 drugs x 16 lineages
########################################################

rownames(gdsc1_avg_matrix)<-tolower(rownames(gdsc1_avg_matrix))

gdsc1_cancer_overlap<-orr_sorted_cancers[which(orr_sorted_cancers %in% colnames(gdsc1_avg_matrix))]
gdsc1_drug_overlap<-tolower(orr_sorted_drugs[which(tolower(orr_sorted_drugs) %in% rownames(gdsc1_avg_matrix))])

norm_nci1970_overlap<-(1-gdsc1_avg_matrix[gdsc1_drug_overlap,gdsc1_cancer_overlap])
norm_nci1970_overlap_scaled<-t(scale(t(norm_nci1970_overlap),scale=T,center=F))#Divide by Root-mean-square = normalize scale of differential across lineages for each drug


############################
#FINAL FIGURE S4B,D: GDSC CLUSTERED:
############################
hcluster((1-gdsc1_avg_matrix[gdsc1_drug_overlap,]))

hcluster(t(scale(t((1-gdsc1_avg_matrix[gdsc1_drug_overlap,])),scale=T,center=F)))


############################
#Figure S5C
############################
hcluster_unsort(norm_nci1970_overlap_scaled,title="GDSC:  lineage-average",units="(1-AUC) / drug-rt-mean-sq",col_size = 1.3,row_size = 0.9)


########################################################
#Figure S5B,C:  NORMALZIE CLINICAL DATA TO ASSESS CLINICAL & IN VITRO CORRELATION
########################################################
NCI1970_META_ORrate_lowercase<-NCI1970meta_rate
rownames(NCI1970_META_ORrate_lowercase)<-tolower(rownames(NCI1970meta_rate))

norm_CTRP<-(16-ctrp_avg_matrix[ctrp_drug_overlap,ctrp_cancer_overlap])/16
norm_NCI_ctrp<-NCI1970_META_ORrate_lowercase[ctrp_drug_overlap,ctrp_cancer_overlap]

################
##3.2  Overall-AVG sensitivity trends:
################

mean_drug_sens_CTRP<-colMeans(norm_CTRP,na.rm=T)
mean_drug_sens_NCI<-colMeans(norm_NCI_ctrp,na.rm=T)


#BEST FIT LINE:  y = 1.46x -0.2       r^2 = 0.68    rmse = 0.06
linear_fit_plot(mean_drug_sens_CTRP,mean_drug_sens_NCI)


#ACTUAL PLOT:
plot(mean_drug_sens_CTRP,mean_drug_sens_NCI,pch=16,xlim=c(0.2,0.45),ylim=c(0.05,0.5),xlab="AVG in vitro sensitivity (1-AUC)",ylab="Clinical Response Rate")
text(mean_drug_sens_CTRP+0.008,mean_drug_sens_NCI,labels=names(mean_drug_sens_NCI))
abline(a=-0.2,b=1.46)





####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##           ########                                                            
##     ##          ##    ##                          ##   ##                 ##                                   
##     ##          ##   ##                      #################            ##                                       
##     ########    ##   ##   ####                    ##   ##           ########                                              
##     ##          ##   ##      ##              #################      ##                                              
##     ##          ##    ##     ##                  ##   ##           ##                                      
##     ##          ##     #######                   ##   ##           ########                                                     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########     ##     #######   ###      ###       ##        #######     ##     ##  #######   #######      ######                                
##    ##     ##    ##    ##     ##  ####    ####      ####       ##    ##    ##    ##   ##        ##    ##    ##                                         
##    ##     ##    ##    ##     ##  ## ##  ## ##     ##  ##      ##    ##    ##   ##    ##        ##    ##    ##                                                
##    ########     ##    ##     ##  ##  ####  ##    ##    ##     #######     ######     #######   #######      ######                                                   
##    ##     ##    ##    ##     ##  ##   ##   ##   ##########    ##   ##     ##   ##    ##        ##   ##           ##                                        
##    ##     ##    ##    ##     ##  ##        ##   ##      ##    ##    ##    ##    ##   ##        ##    ##          ##                       
##    ########     ##     #######   ##        ##   ##      ##    ##     ##   ##     ##  #######   ##     ##    ######                                                                        
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 

#################################################################################################################################################################################################################################################################### 
#FIGURE 2A & S6:   define biomarkers
#################################################################################################################################################################################################################################################################### 
#MECHANISM:  LOAD Corerlation Genes:
load("expr_literature/EXPR_causal-biomarkers.RData")

hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",key="z-score"){
        #distance matrix:
        dist_columns <- dist(t(matrix))
        dist_rows <- dist(matrix)
        
        #hclustering:
        hclust_columns<-hclust(dist_columns, method="average")
        hclust_rows<-hclust(dist_rows, method="average")
        
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), 
                  cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("white","black"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = key,key.xlab = key)
}


#Figure S6A:  drug-targets
hcluster(DrugBank_targets_NCI1970)

#Figure S6B:  drug-transporters
hcluster(DrugBank_transport_NCI1970)

#Figure S6C:  drug-metabolism
hcluster(DrugBank_metab_NCI1970)

#Figure S6D:  indirect drug-associations
hcluster(EFD_repair_NCI1970)

#################################################################################################################################################################################################################################################################### 
#FIGURE 2C-D:   biomarker set correlations
#################################################################################################################################################################################################################################################################### 

#Mechanism of Action:   
load("expr_literature/BIOMARKER-signs_matrices.RData")

#LOAD LINEAGE AVG DRUGS:
load("cell-line_dbs/AUC_AVG-zscore.RData")
load("NCI1970meta_Database.RData")
rownames(NCI1970meta_rate)<-tolower(rownames(NCI1970meta_rate))
rm(NCI1970meta_respond)
rm(NCI1970meta_total)

#LOAD LINEAGE AVG RNA:
load("cell-line_dbs/RNAseq_AVG-zscore.RData")
colnames(tcga_avg_matrix)<-gsub("luad","nsclc",gsub("dlbc","lymph",gsub("laml","leuk",colnames(tcga_avg_matrix))))
load("cell-line_dbs/RNAseq_AVG-tpm.RData")
colnames(tcga_avg_tpm)<-gsub("luad","nsclc",gsub("dlbc","lymph",gsub("laml","leuk",colnames(tcga_avg_tpm))))

#LOAD INDIVIDUAL CELL-LINES:
load("cell-line_dbs/AUC_INDIV_zscore.RData")

load("cell-line_dbs/CCLE_RNAseq.RData")
ccle_RNAseq_zscore<-t(scale(t(ccle_RNAseq_filtered),scale=T,center=T))

###################################
#VISUALIZATION FUNCTIONS:
###################################
#plot individual correlations
correlation_drugs_fuction<-function(biomarkers,expression,drugAUC,plot=T){
        gene_overlap<-intersect(colnames(biomarkers),rownames(expression))
        
        MoA_predict_transport<-biomarkers[,gene_overlap] %*% expression[gene_overlap,]
        rownames(MoA_predict_transport)<-tolower(rownames(MoA_predict_transport))
        
        ##################################3
        #RANK DRUG-DRUG CORRELATION
        ##################################
        drugs_overlap<-intersect(rownames(drugAUC),rownames(MoA_predict_transport))
        lineage_overlap<-intersect(colnames(drugAUC),colnames(MoA_predict_transport))
        
        corr_matrix_invitro<-cor(t(MoA_predict_transport[drugs_overlap,lineage_overlap]),
                                 t(drugAUC[drugs_overlap,lineage_overlap]),use="pairwise.complete.obs")
        
        invitro_tissue_corr<-c()
        for (drug in drugs_overlap){
                tmp_corr<-corr_matrix_invitro[drug,drug]
                invitro_tissue_corr<-c(invitro_tissue_corr,tmp_corr)
        }
        
        names(invitro_tissue_corr)<-drugs_overlap
        
        invitro_tissue_corr_sorted<-sort(invitro_tissue_corr,decreasing=F)
        if(plot==T){
                par(las=2,mar=c(4,6,2,2))
                barplot(invitro_tissue_corr_sorted,horiz=T,cex.names = 0.8,xlim=c(-0.3,0.8))
                par(las=1,mar=(c(5, 4, 4, 2) + 0.1))
        }
        
        
        return(invitro_tissue_corr)
}

#visualize total-biomarker correlation colored in proportion to correaltions of biomarker-classes
correlation_MoA_class<-function(expression,drugAUC){
        
        #ALL BIOMARKERS:
        tmp_all<-correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,expression,drugAUC,plot=F)
        
        #INDIVIDUAL BIOMARKER TYPES:
        tmp_trans<-correlation_drugs_fuction(chemo_transport_matrix,expression,drugAUC,plot=F)
        tmp_metab<-correlation_drugs_fuction(chemo_metab_matrix,expression,drugAUC,plot=F)
        tmp_targt<-correlation_drugs_fuction(chemo_target_matrix,expression,drugAUC,plot=F)
        tmp_ddr<-correlation_drugs_fuction(chemo_ddr_matrix,expression,drugAUC,plot=F)
        
        moa_matrix<-abs(cbind(tmp_trans[names(tmp_all)],
                              tmp_metab[names(tmp_all)],
                              tmp_targt[names(tmp_all)],
                              tmp_ddr[names(tmp_all)]))
        colnames(moa_matrix)<-c("trans","metab","targets","indirect")
        
        #make adjustment matrix
        rescale_col<-tmp_all/rowSums(moa_matrix,na.rm=T)
        
        moa_matrix_na<-moa_matrix
        moa_matrix_na[is.na(moa_matrix_na)]<-0
        
        moa_rescaled<-moa_matrix_na*cbind(rescale_col,rescale_col,rescale_col,rescale_col)
        
        
        #check stacked par plot (moa_rescaled) vs actual total corr (tmp_all):
        rowSums(moa_rescaled)
        tmp_all
        
        
        #PLOT:
        final_plot<-t(moa_rescaled[names(sort(rowSums(moa_rescaled))),])
        
        par(las=2,mar=c(4,6,2,2))
        barplot(final_plot,horiz=T,cex.names = 0.8,xlim=c(-0.3,0.9),col=c("darkred","darkred","darkgreen","black"))
        par(las=1,mar=(c(5, 4, 4, 2) + 0.1))
        
        
}


###################################
#Figure 2C-D
###################################
correlation_MoA_class(tcga_avg_tpm,NCI1970meta_rate)

correlation_MoA_class(ccle_avg_tpm,
                      ctrp_avg_matrix)


###############################################################################################################################################################################
#PLOT INDIVIDUAL correlations:
###############################################################################################################################################################################

######################
#INDIVIDUAL LINES
######################

correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_RNAseq_filtered,
                          ctrp_zscore)
correlation_MoA_class(ccle_RNAseq_filtered,ctrp_zscore)




correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_RNAseq_filtered,
                          prism_zscore)
correlation_MoA_class(ccle_RNAseq_filtered,prism_zscore)


correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_RNAseq_filtered,
                          gdsc1_zscore)

correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_RNAseq_filtered,
                          gdsc2_zscore)

######################
#LINEAGE AVERAGES:
######################

#CLINICAL
correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          tcga_avg_tpm,
                          NCI1970meta_rate)

correlation_drugs_fuction(chemo_metab_matrix,
                          tcga_avg_tpm,
                          NCI1970meta_rate)



#IN VITRO
correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_avg_tpm,
                          ctrp_avg_matrix)

correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_avg_tpm,
                          prism_avg_matrix)

correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_avg_tpm,
                          gdsc1_avg_matrix)

correlation_drugs_fuction(chemo_BIOMARKERS_matrix_DDR,
                          ccle_avg_tpm,
                          gdsc2_avg_matrix)









####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##           ########      ########         #######                                                  
##     ##          ##    ##                          ##   ##          ##                   ##        ##      ##          
##     ##          ##   ##                      #################     ##                  ##         ##      ##              
##     ########    ##   ##   ####                    ##   ##           ########          ##   ######  #########                              
##     ##          ##   ##      ##              #################             ##        ##                   ##              
##     ##          ##    ##     ##                  ##   ##                   ##       ##                   ##       
##     ##          ##     #######                   ##   ##            ########       ##                   ##                           
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########     ##     #######   ###      ###       ##        #######     ##     ##  #######   #######                                   
##    ##     ##    ##    ##     ##  ####    ####      ####       ##    ##    ##    ##   ##        ##    ##                                            
##    ##     ##    ##    ##     ##  ## ##  ## ##     ##  ##      ##    ##    ##   ##    ##        ##    ##                                                  
##    ########     ##    ##     ##  ##  ####  ##    ##    ##     #######     ######     #######   #######                                                       
##    ##     ##    ##    ##     ##  ##   ##   ##   ##########    ##   ##     ##   ##    ##        ##   ##                                                  
##    ##     ##    ##    ##     ##  ##        ##   ##      ##    ##    ##    ##    ##   ##        ##    ##                               
##    ########     ##     #######   ##        ##   ##      ##    ##     ##   ##     ##  #######   ##     ##                                                                          
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########    #######     ########      #######    ##              ##      ########### ########    ########   ####    ##                                                                            
##     ##          ##     ##    ##     ##     ##    ##   ##             ####         ##         ##      ##      ##  ## ##   ##                
##    ##           ##     ##    ##     ##     ##    ##   ##            ##  ##        ##         ##      ##      ##  ##  ##  ##                         
##    ##           ##     ##    ########      #######    ##           ##    ##       ##         ##      ##      ##  ##   ## ##                 
##    ##           ##     ##    ##    ##      ##   ##    ##          ##########      ##         ##      ##      ##  ##    ####                    
##     ##          ##     ##    ##     ##     ##    ##   ##         ##        ##     ##         ##      ##      ##  ##     ###                  
##      ########    #######     ##      ##    ##     ##  ########  ##          ##    ##      ########    ########   ##      ##                             
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########   #########   ##########    #######                                  
##    ##            ##           ##       ##          
##    ##            ##           ##       ##                                
##    #######       ##           ##        ######    
##    ##            ##           ##             ##    
##    ##            ##           ##             ##                  
##    ##         ########        ##       #######
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
#STRUCTURE:
#       #0 LOAD DATA:  drug sensitivity, biomarker-definitions, biomarker expression
#       FIGURE 3A: Vinca Alkaloid Biomarker Analysis:  in vitro vs clinical 
#       FIGURE 3B: Taxane Biomarker Analysis:  in vitro vs clinical 
#       FIGURE 3C: Antrhacycline Biomarker Analysis:  in vitro vs clinical
#
#################################################################################################################################################################################################################################################################### 


####################################################################################################################################################################################################################################################################
#0 LOAD DATA
####################################################################################################################################################################################################################################################################



##############################
#SENSITIVITY DATA:
###############################
#CLINICAL DATA + DRUGS TO BE ANALYZED
load("NCI1970meta_Database.RData")
rm(NCI1970meta_total)
rm(NCI1970meta_respond)


#INDIVIDUAL CELL LINES:
load("cell-line_dbs/invitro_AUC-lines-individ_dbs.RData")
rm(gdsc1_auc_matrix)
rm(gdsc2_auc_matrix)
rm(prism_auc_matrix)
rm(prism_ec50_matrix)


##############################
#MECHANISM DATA:
###############################
#MECHANISM:  LOAD Corerlation Genes:
load("expr_literature/EXPR_causal-biomarkers.RData")
rm(DrugBank_metab_NCI1970)
rm(DrugBank_targets_NCI1970)
rm(DrugBank_transport_NCI1970)
rm(EFD_repair_NCI1970)



##############################
#GENE EXPRESSION:
###############################

#INDIVIDUAL LINE CCLE:
load("cell-line_dbs/CCLE_RNAseq.RData")

#AVG LINEAGE TCGAE:
load("clinical_dbs/TCGA_RNAseq_lineageAVG.RData")



#LINEAR REGRESSION:
linear_fit_plot<-function(x,y,pt_colors="black"){
        #linear regression:
        lin_reg<-lm(y~x)
        #make fit equation + r-squared:
        pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
        fit_coeff<-round(coef(lin_reg),6)
        r2 <- round(summary(lin_reg)$r.squared, 2)
        rmse <- round(sqrt(mean(resid(lin_reg)^2)), 2)
        eq_r2<-paste("y = ", fit_coeff[2], "x + ", fit_coeff[1]," "," "," ","r^2 = ", r2," "," "," ","rmse = ",rmse)
        #plot data:
        plot(x,y, pch = 16, cex = 1.3, col = pt_colors)
        #add fit-line and fit-equation:
        abline(lin_reg)
        mtext(eq_r2, 3, line=-2,cex=1)
}



####################################################################################################################################################################################################################################################################
#FIGURE 3A:  Vinca Alkalod Biomarkers 
####################################################################################################################################################################################################################################################################
#NOTE:   vindesine (most data)  vincristine (most resist genes)

##################
#CLINCIAL CORR:
##################
#RECONCILE CANCER LABELS:
vinblastine_ORR_vector_raw<-NCI1970meta_rate["Vindesine",-which(is.na(NCI1970meta_rate["Vindesine",])==T)]
names(vinblastine_ORR_vector_raw)<-gsub("nsclc","luad",gsub("lymph","dlbc",gsub("leuk","laml",names(vinblastine_ORR_vector_raw))))


#EXTRACT BIOMARKERS:
vinb_resist_genes<-names(which(MASTER_MoA_matrix_efd["Vincristine",]==1))#18 resist genes

vinblastine_TCGA_moa_genes<-t(TCGA_RNAseq_lineageAVG[vinb_resist_genes,names(vinblastine_ORR_vector_raw)])


#CALCULATE CORRELATION:
dataframe_multiLinReg<-cbind("ORR"=vinblastine_ORR_vector_raw,vinblastine_TCGA_moa_genes[names(vinblastine_ORR_vector_raw),])

vinb_ORRgene_rank<-sort(cor(dataframe_multiLinReg)["ORR",],decreasing=T)
#        ORR       BRCA1       ABCC1        TUBB      RALBP1      TUBA4A      ABCC10      ABCB11       ABCB4       ABCB1      CYP3A4       ABCC2       ABCC3     SLC22A3      CYP3A5 
#1.00000000  0.53930040  0.39224267  0.34009033  0.17358267  0.08797618  0.03801962 -0.07544003 -0.28434574 -0.37432998 -0.43334934 -0.44793733 -0.49351752 -0.52203949 -0.61099564 
#CYP3A7       ABCG2 
#-0.64291334 -0.74443048



##################
#IN VITRO CORR:
##################
#EXTRACT DATA:
vinblastine_AUC_vector_raw<-16-ctrp_auc_matrix["vincristine",-which(is.na(ctrp_auc_matrix["vincristine",])==T)]

#EXTRACT BIOMARKERS:
vinb_resist_genes<-names(which(MASTER_MoA_matrix_efd["Vincristine",]==1))


#CALCULATE CORRELATION:
line_overlap<-names(vinblastine_AUC_vector_raw)[which(names(vinblastine_AUC_vector_raw) %in% colnames(ccle_RNAseq_filtered))]
vinblastine_CCLE_moa_genes<-t(ccle_RNAseq_filtered[vinb_resist_genes,line_overlap])



vinb_AUC_multiLinReg<-cbind("AUC"=vinblastine_AUC_vector_raw[line_overlap],vinblastine_CCLE_moa_genes[line_overlap,vinb_resist_genes])
vinb_AUCgene_rank<-sort(cor(vinb_AUC_multiLinReg)["AUC",],decreasing=T)
#        AUC       BRCA1        TUBB      RALBP1      ABCC10       ABCB4      CYP3A7      CYP3A4      ABCB11      TUBA4A       ABCC1     SLCO1B1       ABCB1      CYP3A5       ABCG2 
#1.00000000  0.26540999  0.11352532  0.10514109  0.06354762  0.05843371 -0.03108048 -0.04274429 -0.05376229 -0.08262648 -0.10880927 -0.11305797 -0.12767877 -0.22245782 -0.22577887 
#SLC22A3     SLCO1B3       ABCC2       ABCC3 
#-0.22684473 -0.22892966 -0.24562246 -0.50339445 

######################################################
#PLOTS: clinical & in vitro aggreement
######################################################
linear_fit_plot(vinb_AUCgene_rank,vinb_ORRgene_rank[names(vinb_AUCgene_rank)])

# Figure 3A

#PLOT CORR
plot(vinb_AUCgene_rank,vinb_ORRgene_rank[names(vinb_AUCgene_rank)],main="Vinca Alkaloids: in vitro vs clinical (r2 = 0.46)",
     pch=16,xlim=c(-0.65,0.65),ylim=c(-0.65,0.65),xlab="In vitro correlation",ylab="Clinical correlation")
abline(h=0,v=0,lty=2,lwd=0.5)
abline(a=-0.072,b=1.5,col="grey")

text(vinb_AUCgene_rank+0.08,vinb_ORRgene_rank[names(vinb_AUCgene_rank)],labels=names(vinb_AUCgene_rank))



######################################################
#PLOTS: multi-linear regression
######################################################
vinb_AUC_multiLinReg_df<-data.frame(vinb_AUC_multiLinReg)
#MULTIPLE LINEAR REGRESSION:
vinb_AUC_model<-lm(AUC~BRCA1+TUBB+RALBP1+ABCC10+ABCB4+CYP3A7+CYP3A4+ABCB11+TUBA4A+ABCC1+SLCO1B1+ABCB1+CYP3A5+ABCG2+SLC22A3+SLCO1B3+ABCC2+ABCC3,data=data.frame(vinb_AUC_multiLinReg))

summary(vinb_AUC_model)
#Coefficients:
        #  Estimate Std. Error t value Pr(>|t|)    
        #(Intercept)  9.393824   1.673263   5.614 2.73e-08 ***
        #BRCA1        0.837448   0.155505   5.385 9.53e-08 ***
        #TUBB        -0.450784   0.150922  -2.987  0.00291 ** 
        #TUBA4A       0.172050   0.053084   3.241  0.00124 ** 
        #ABCB1       -0.183544   0.070037  -2.621  0.00894 ** 
        #SLCO1B3     -0.230545   0.084095  -2.741  0.00625 ** 
        #ABCC3       -0.670582   0.060048 -11.167  < 2e-16 ***
        
model_predict<-vinb_AUC_multiLinReg_df$BRCA1*0.84-vinb_AUC_multiLinReg_df$TUBB*0.45-vinb_AUC_multiLinReg_df$ABCB1*0.18-0.67*vinb_AUC_multiLinReg_df$ABCC3



linear_fit_plot(model_predict,vinb_AUC_multiLinReg_df$AUC)








####################################################################################################################################################################################################################################################################
#FIGURE 3B:  Taxane Biomarkers Biomarkers 
####################################################################################################################################################################################################################################################################

##################
#CLINCIAL CORR:
##################
#RECONCILE CANCER LABELS:
taxol_ORR_vector_raw<-NCI1970meta_rate["Paclitaxel",-which(is.na(NCI1970meta_rate["Paclitaxel",])==T)]
names(taxol_ORR_vector_raw)<-gsub("nsclc","luad",gsub("lymph","dlbc",gsub("leuk","laml",names(taxol_ORR_vector_raw))))

#EXTRACT BIOMARKERS:
taxol_resist_genes<-c(names(which(MASTER_MoA_matrix_efd["Paclitaxel",]==1)),"TUBB")

taxol_TCGA_moa_genes<-t(TCGA_RNAseq_lineageAVG[taxol_resist_genes,names(taxol_ORR_vector_raw)])


#SETUP LINEAR REGRESSION:
dataframe_multiLinReg_taxol<-cbind("ORR"=taxol_ORR_vector_raw,taxol_TCGA_moa_genes[names(taxol_ORR_vector_raw),])

taxol_ORRgene_rank<-sort(cor(dataframe_multiLinReg_taxol)["ORR",],decreasing=T)
#        ORR       BRCA1        TUBB       ABCC1      ABCC10     SLCO1B3     CYP19A1      ABCB11        BCL2      CYP1B1        MAP4       TUBB1       NR1I2        MAPT      CYP2C8 
#1.00000000  0.72906943  0.68954791  0.33012755  0.08429380  0.06243158 -0.02759440 -0.09922976 -0.11696356 -0.28497052 -0.28829713 -0.28892485 -0.32474957 -0.42504056 -0.56084312 
#ABCC2       ABCB1      CYP3A4      CYP3A5        MAP2      CYP3A7 
#-0.60646117 -0.66669470 -0.69845768 -0.71572287 -0.76615441 -0.80235210 



##################
#IN VITRO CORR: 
##################
#EXTRACT DATA:
taxol_AUC_vector_raw<-16-ctrp_auc_matrix["paclitaxel",-which(is.na(ctrp_auc_matrix["paclitaxel",])==T)]

#EXTRACT BIOMARKERS:
taxol_resist_genes<-c(names(which(MASTER_MoA_matrix_efd["Paclitaxel",]==1)),"TUBB")


#CALCULATE CORRELATION:
line_overlap<-names(taxol_AUC_vector_raw)[which(names(taxol_AUC_vector_raw) %in% colnames(ccle_RNAseq_filtered))]
taxol_CCLE_moa_genes<-t(ccle_RNAseq_filtered[taxol_resist_genes,line_overlap])


taxol_AUC_multiLinReg<-cbind("AUC"=taxol_AUC_vector_raw[line_overlap],taxol_CCLE_moa_genes[line_overlap,taxol_resist_genes])
taxol_AUCgene_rank<-sort(cor(taxol_AUC_multiLinReg)["AUC",],decreasing=T)
#AUC         BCL2        BRCA1        TUBB1         TUBB       ABCC10        NR1I2       CYP2C8         MAPT        ABCC1       CYP3A7       ABCB11       CYP3A4      CYP19A1 
#1.000000000  0.291749633  0.235046045  0.125384316  0.065134398  0.028911675  0.006336948 -0.007620625 -0.011128508 -0.021995739 -0.022055626 -0.030678385 -0.046117704 -0.063146738 
#MAP2      SLCO1B3         MAP4       CYP1B1       CYP3A5        ABCB1        ABCC2 
#-0.076388545 -0.111999614 -0.149774088 -0.165596154 -0.169897156 -0.194089126 -0.235534541






######################################################
#PLOTS: clinical & in vitro aggreement
######################################################
linear_fit_plot(taxol_AUCgene_rank,taxol_ORRgene_rank[names(taxol_AUCgene_rank)])#r2 = 0.27

# Figure 3B

#PLOT CORR
plot(taxol_AUCgene_rank,taxol_ORRgene_rank[names(taxol_AUCgene_rank)],main="Taxanes: in vitro vs clinical (r2 = 0.27)",
     pch=16,xlim=c(-0.35,0.35),ylim=c(-0.8,0.8),xlab="In vitro correlation",ylab="Clinical correlation")
abline(h=0,v=0,lty=2,lwd=0.5)
abline(a=-0.2,b=1.76,col="grey")

text(taxol_AUCgene_rank+0.05,taxol_ORRgene_rank[names(taxol_AUCgene_rank)],labels=names(taxol_AUCgene_rank))




######################################################
#PLOTS: multi-linear regression
######################################################
taxol_AUC_multiLinReg_df<-data.frame(taxol_AUC_multiLinReg)


taxol_AUC_model<-lm(AUC~NR1I2+TUBB1+MAPT+BCL2+MAP2+MAP4+ABCB1+ABCC2+ABCB11+ABCC10+SLCO1B3+ABCC1+CYP3A4+CYP3A5+CYP1B1+CYP2C8+CYP3A7+CYP19A1+BRCA1+TUBB,data=data.frame(taxol_AUC_multiLinReg))

summary(taxol_AUC_model)
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  9.14645    1.84002   4.971 8.26e-07 ***
#NR1I2        0.46232    0.18432   2.508 0.012345 *  
#TUBB1        0.70438    0.34689   2.031 0.042648 *  
#BCL2         0.43421    0.09290   4.674 3.50e-06 ***
#MAP4        -0.47917    0.18257  -2.625 0.008853 ** 
#ABCB1       -0.40797    0.08326  -4.900 1.18e-06 ***
#ABCC2       -0.29874    0.06815  -4.384 1.33e-05 ***
#CYP3A5      -0.22223    0.10709  -2.075 0.038306 *  
#BRCA1        0.66285    0.19195   3.453 0.000585 ***

taxol_model_predict<-(9.14645+taxol_AUC_multiLinReg_df$TUBB1*0.70438+taxol_AUC_multiLinReg_df$BCL2*0.43421+taxol_AUC_multiLinReg_df$MAP4*-0.47917+taxol_AUC_multiLinReg_df$ABCB1*-0.40797+taxol_AUC_multiLinReg_df$ABCC2*-0.29874+taxol_AUC_multiLinReg_df$CYP3A5*-0.22223+taxol_AUC_multiLinReg_df$BRCA1*0.66285)



linear_fit_plot(taxol_model_predict,taxol_AUC_multiLinReg_df$AUC)





####################################################################################################################################################################################################################################################################
#FIGURE 3C:  Anthracycline Biomarkers 
####################################################################################################################################################################################################################################################################


##################
#CLINCIAL CORR:
##################
#RECONCILE CANCER LABELS:
dox_ORR_vector_raw<-NCI1970meta_rate["Doxorubicin",-which(is.na(NCI1970meta_rate["Doxorubicin",])==T)]
names(dox_ORR_vector_raw)<-gsub("nsclc","luad",gsub("lymph","dlbc",gsub("leuk","laml",names(dox_ORR_vector_raw))))

#EXTRACT BIOMARKERS:
dox_resist_genes<-names(which(MASTER_MoA_matrix_efd["Doxorubicin",]==1))

dox_TCGA_moa_genes<-t(TCGA_RNAseq_lineageAVG[dox_resist_genes,names(dox_ORR_vector_raw)])


#CALCULATE CORRELATION:
dataframe_multiLinReg_dox<-cbind("ORR"=dox_ORR_vector_raw,dox_TCGA_moa_genes[names(dox_ORR_vector_raw),])

dox_ORRgene_rank<-sort(cor(dataframe_multiLinReg_dox)["ORR",],decreasing=T)
#ORR       CHEK2    SLC22A16        MSH2       RAD51        TP53        MLH1       TOP2A       NOLC1       PRKDC         ATM       XRCC6       ABCC1      NDUFS3      CYP2D6 
#1.00000000  0.62488635  0.59048558  0.58996769  0.48892027  0.44992355  0.39555184  0.34441185  0.32490970  0.31477197  0.31423193  0.22082683  0.19751746  0.14844874  0.12343663 
#AKR1A1        LIG4      CYP1B1      ABCC10      RALBP1        TDP2       ABCB8       ABCC6      ABCB11      NDUFS7      NDUFS2      CYP3A4       RRM2B       ABCB4       ABCB1 
#0.11775352  0.11569707  0.08841421  0.08786283 -0.07045425 -0.08193288 -0.08945862 -0.09365632 -0.10703404 -0.14253261 -0.14441611 -0.14861304 -0.16753224 -0.22869237 -0.24728943 
#ABCC2        CBR3       ABCG2       ABCB5      AKR1C3         POR       ABCC3        NOS2        NOS3        NQO1        CBR1 
#-0.27010280 -0.28884777 -0.30067904 -0.35476856 -0.35588765 -0.41885995 -0.49351272 -0.54520006 -0.56931782 -0.61493928 -0.64452539







##################
#IN VITRO: ctrp 804 lines
##################
#EXTRACT DATA:
dox_AUC_vector_raw<-16-ctrp_auc_matrix["doxorubicin",-which(is.na(ctrp_auc_matrix["doxorubicin",])==T)]

#EXTRACT BIOMARKERS:
dox_resist_genes<-names(which(MASTER_MoA_matrix_efd["Doxorubicin",]==1))

#CALCULATE CORRELATION:
line_overlap<-names(dox_AUC_vector_raw)[which(names(dox_AUC_vector_raw) %in% colnames(ccle_RNAseq_filtered))]
dox_CCLE_moa_genes<-t(ccle_RNAseq_filtered[dox_resist_genes,line_overlap])



dox_AUC_multiLinReg<-cbind("AUC"=dox_AUC_vector_raw[line_overlap],dox_CCLE_moa_genes[line_overlap,dox_resist_genes])
dox_AUCgene_rank<-sort(cor(dox_AUC_multiLinReg)["AUC",],decreasing=T)
#         AUC     SLC22A16          ATM       NDUFS7        CHEK2        RAD51        PRKDC       AKR1A1        NOLC1       CYP2D6         TP53         LIG4        TOP2A       NDUFS3 
#1.000000000  0.333633518  0.266793430  0.265862468  0.236430090  0.227346565  0.211715397  0.196765578  0.183868559  0.183410906  0.182825896  0.168043442  0.167875537  0.160929672 
#NDUFS2         MSH2        ABCB4         NOS3         MLH1       ABCC10        XRCC6       RALBP1         TDP2         NOS2        RRM2B       CYP3A4         NOS1        ABCB5 
#0.113468643  0.104126162  0.102041102  0.095386545  0.084039976  0.073605660  0.072047803  0.069416212  0.065377653  0.061632385  0.046570551  0.007052794 -0.033313191 -0.057460111 
#ABCB11        ABCB8       CYP2B6        ABCB1        ABCC1        ABCC6        ABCG2          XDH       CYP1B1        ABCC2         CBR3          POR         CBR1       AKR1C3 
#-0.061778231 -0.071872462 -0.147159291 -0.169071858 -0.181513903 -0.203335523 -0.220169715 -0.221599154 -0.242323828 -0.251427016 -0.268329687 -0.306015531 -0.337617758 -0.347641360 
#NQO1        ABCC3 
#-0.478837816 -0.524162445 





######################################################
#PLOTS: clinical & in vitro aggreement
######################################################
linear_fit_plot(dox_AUCgene_rank,dox_ORRgene_rank[names(dox_AUCgene_rank)])#r2 = 0.48

#Figure 3C

#PLOT CORR
plot(dox_AUCgene_rank,dox_ORRgene_rank[names(dox_AUCgene_rank)],main="Anthracyclines: in vitro vs clinical (r2 = 0.48)",
     pch=16,xlim=c(-0.65,0.65),ylim=c(-0.65,0.65),xlab="In vitro correlation",ylab="Clinical correlation")
abline(h=0,v=0,lty=2,lwd=0.5)
abline(a=-0.02,b=1.11,col="grey")

text(dox_AUCgene_rank+0.08,dox_ORRgene_rank[names(dox_AUCgene_rank)],labels=names(dox_AUCgene_rank))




######################################################
#PLOTS: multi-linear regression
######################################################
dox_AUC_multiLinReg_df<-data.frame(dox_AUC_multiLinReg)



dox_AUC_model<-lm(AUC~SLC22A16+ATM+NDUFS7+CHEK2+RAD51+PRKDC+AKR1A1+NOLC1+CYP2D6+TP53+LIG4+TOP2A+
                          NQO1+ABCC3+ABCC2+CBR3+POR+CBR1+AKR1C3+ABCB1+ABCC1+ABCC6+ABCG2+XDH+CYP1B1,data=data.frame(dox_AUC_multiLinReg))

summary(dox_AUC_model)
#Coefficients:
#Estimate Std. Error t value Pr(>|t|)    
#(Intercept) -8.589964   1.288800  -6.665 5.01e-11 ***
#SLC22A16     0.443820   0.087666   5.063 5.17e-07 ***

#NDUFS7       0.249498   0.112780   2.212 0.027239 *  
#CYP2D6       0.582237   0.122452   4.755 2.37e-06 ***
#TP53         0.107837   0.051317   2.101 0.035929 *  
#TOP2A        0.114680   0.099004   1.158 0.247084    

#ABCB1       -0.161288   0.041694  -3.868 0.000119 ***

#NQO1        -0.119559   0.040073  -2.984 0.002939 ** 
#ABCC3       -0.219590   0.042111  -5.215 2.36e-07 ***
#CBR1        -0.104872   0.043757  -2.397 0.016780 *  
#ABCC1       -0.322242   0.081555  -3.951 8.48e-05 ***
#ABCC6       -0.168003   0.065705  -2.557 0.010749 *  


dox_model_predict<-(dox_AUC_multiLinReg_df$SLC22A16*0.443820+dox_AUC_multiLinReg_df$NDUFS7*0.249498+dox_AUC_multiLinReg_df$TP53*0.107837+dox_AUC_multiLinReg_df$TOP2A*0.114680
                    +dox_AUC_multiLinReg_df$ABCB1*-0.161288+dox_AUC_multiLinReg_df$NQO1*-0.119559+dox_AUC_multiLinReg_df$ABCC3*-0.219590+dox_AUC_multiLinReg_df$CBR1*-0.104872+dox_AUC_multiLinReg_df$ABCC1*-0.322242+dox_AUC_multiLinReg_df$ABCC6*-0.168003)

linear_fit_plot(dox_model_predict,dox_AUC_multiLinReg_df$AUC)








####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##           ########    ##    #######                                               
##     ##          ##    ##                          ##   ##          ##           ##   ##     ##         
##     ##          ##   ##                      #################     ##           ##   ##     ##           
##     ########    ##   ##   ####                    ##   ##           ########    ##   ##     ##                            
##     ##          ##   ##      ##              #################             ##   ##   ##     ##           
##     ##          ##    ##     ##                  ##   ##                   ##   ##   ##     ##
##     ##          ##     #######                   ##   ##            ########    ##    #######                       
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########      ######     ######   ##    ##  #######   #######      ######                                                             
##     ##     ##    ##    ##   ##    ##  ##    ##  ##        ##    ##    ##                                                 
##     ##     ##    ##    ##   ##    ##  ##    ##  ##        ##    ##    ##                                              
##     ########     ##    ##   ##    ##  ##    ##  #######   #######      ######                                                     
##     ##    ##     ##    ##   ##    ##   ##  ##   ##        ##   ##           ##                                        
##     ##     ##    ##    ##   ##    ##    ####    ##        ##    ##          ##                                          
##     ##      ##    ######     ######      ##     #######   ##     ##    ######                                                      
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
##   ##     ##   ##   ###     ##   ########  ##########  ##     #######   #######                                                                                                      
##   ##    ##    ##   ####    ##   ##            ##      ##    ##        ##                                                                                  
##   ##   ##     ##   ## ##   ##   ##            ##      ##   ##         ##                                                                                     
##   ######      ##   ##  ##  ##   #######       ##      ##   ##          ######                                                                                           
##   ##   ##     ##   ##   ## ##   ##            ##      ##   ##               ##                                                                               
##   ##    ##    ##   ##    ####   ##            ##      ##    ##              ##                                                                                
##   ##     ##   ##   ##     ###   #######       ##      ##     ######    ######                                                                                            
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##    ##     ########                              
##    ##    ##          
##    ##   ##            ####   ####     
##    ##   ##            #     ##  ##
##    ##   ##            ####  ##  ##
##    ##    ##              ## ##  ##   
##    ##      ########   ####   ####                                    
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 



#KINETICS DRUG-UPTAKE DATA AND PLOTTING FUNCTION"
ida_transport<-read.csv("/Users/eugenedouglass/Dropbox/00 Independent work/--- 02 ORIGINAL RESEARCH/2021-06-11 - MDR biochem/04 IDArubicin data extract/raw-idarubicin-KINETIC-data-v05.csv",as.is=T,header=T,row.names = 1)[1:8,]

smooth_influx<-function(x,y,color="black",width=1){
        fitmodel <- nls(y~y_change*(1-exp(-k_trans*x)), start=list(y_change=50,k_trans=0.1),control=nls.control(maxiter = 100, warnOnly = T))
        xl <- seq(min(x),max(x), (max(x) - min(x))/100)
        lines(xl, predict(fitmodel,list(x=xl)), col=color, lwd=width)
        
}
smooth_influx_dotted<-function(x,y,color="black",width=1){
        fitmodel <- nls(y~y_change*(1-exp(-k_trans*x)), start=list(y_change=50,k_trans=0.1),control=nls.control(maxiter = 100, warnOnly = T))
        xl <- seq(min(x),max(x), (max(x) - min(x))/100)
        lines(xl, predict(fitmodel,list(x=xl)), col=color, lwd=width,lty=2)
        
}


#IC50 DATA AND PLOTTING FUNCTION"
ida_IC50<-read.csv("anthracycline_data/raw-idarubicin-IC50-data-v05.csv",as.is=T,header=T,row.names = 1)[-10,]

smooth_hill<-function(x,y,color="black",width=1){
        fitmodel <- nls(y~y_change*IC50/(IC50 + x)+y_min, start=list(y_change=1,IC50=0.1,y_min=0.1))
        xl <- 10^seq(log10(min(x)),log10(max(x)), (log10(max(x)) - log10(min(x)))/100)
        lines(xl, predict(fitmodel,list(x=xl)), col=color, lwd=width)
        
}


#######################################################################################
#Figure S10C:  IDARUBICIN UPDATE IN 3 CELL-LINES:
#######################################################################################



x=as.numeric(rownames(ida_transport))


#S
y=ida_transport$S_IDA
plot(x,y,pch=16,main="Idarubicin Influx:  3 cell lines",col="darkgreen")
smooth_influx(x,y,color="darkgreen",width=2)

#R7:
y=ida_transport$R7_IDA
points(x,y,pch=16,col="navy")
smooth_influx(x,y,color="navy",width=2)

#Dox40:
y=ida_transport$D40_IDA
points(x,y,pch=16,col="darkred")
smooth_influx(x,y,col="darkred",width=2)



#############################
#IDARUBICIN + VERAPAMIL
#############################

#S
y=ida_transport$S_IDAver
points(x,y,pch=1,col="darkgreen",cex=0.6)
smooth_influx_dotted(x,y,color="darkgreen")

#R7:
y=ida_transport$R7_IDAver
points(x,y,pch=1,col="navy",cex=0.6)
smooth_influx_dotted(x,y,color="navy")

#Dox40:
y=ida_transport$D40_IDAver
points(x,y,pch=1,col="darkred",cex=0.6)
smooth_influx_dotted(x,y,col="darkred")





#######################################################################################
#Figure S10D:  IDARUBICIN IC50's +/- verapamil IN 3 CELL-LINES:
#######################################################################################


x=as.numeric(rownames(ida_IC50))

#######################################################################################
#8226-Dox40:
#######################################################################################

#DOX40:  Dox
y=ida_IC50$Dox40_DNR/90
plot(x,y,log="x",ylim=c(0,1.2),pch=16,main="8226-Dox40:  Daunorubicin +/- Verapamil")
smooth_hill(x,y)
#DOX40:  Dox + Verapamil
y=ida_IC50$Dox40_DNRver/65
points(x,y,pch=16)
smooth_hill(x,y)





#######################################################################################
#8226-R7:
#######################################################################################

#DOX40:  Dox
y=ida_IC50$R_DNR/80
plot(x,y,log="x",ylim=c(0,1.2),pch=16,main="8226-R7:  Daunorubicin +/- Verapamil")
smooth_hill(x,y)
#DOX40:  Dox + Verapamil
y=ida_IC50$R_DNRver/90
points(x,y,pch=16)
smooth_hill(x,y)




#######################################################################################
#8226-S:
#######################################################################################

#DOX40:  Dox
y=ida_IC50$S_DNR/60
plot(x,y,log="x",ylim=c(0,1.2),pch=16,main="8226-S:  Daunorubicin +/- Verapamil")
smooth_hill(x,y)
#DOX40:  Dox + Verapamil
y=ida_IC50$S_DNRver/60
points(x,y,pch=16)
smooth_hill(x,y)






####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##           ########                                                            
##     ##          ##    ##                          ##   ##                 ##                                   
##     ##          ##   ##                      #################            ##                                       
##     ########    ##   ##   ####                    ##   ##           ########                                              
##     ##          ##   ##      ##              #################            ##                                              
##     ##          ##    ##     ##                  ##   ##                  ##                                      
##     ##          ##     #######                   ##   ##           ########                                                     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     #######    ########                                                                                                          
##           ##   ##     ##                       
##           ##   ##      ##                                      
##     #######    ##      ##     
##    ##          ##      ##                               
##    ##          #      ##                                     
##     #######    ########                                                         
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########     ########    ##     ##     ########             ##        ##  #########    #######   ##     ##                                                                                                               
##    ##     ##    ##     ##   ##     ##    ##                    ###      ###  ##          ##         ##     ##                                     
##    ##      ##   ##     ##   ##     ##   ##                     ####    ####  ##         ##          ##     ##                                                         
##    ##      ##   ########    ##     ##   ##    ####    ######   ## ##  ## ##  ########   ##          #########                  
##    ##      ##   ##    ##    ##     ##   ##       ##            ##  ####  ##  ##         ##          ##     ##                                                
##    ##     ##    ##     ##   ##     ##    ##      ##            ##   ##   ##  ##          ##         ##     ##                                                       
##    ########     ##      ##   #######      ########             ##        ##  #########    #######   ##     ##                                                                 
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################




library("gplots")
library("RColorBrewer")
library("Hmisc")

#SUPPORTING FUNCTIONS:
make_lineage_avg<-function(matrix){
        lineage_ordered_ccle<-sapply(strsplit(colnames(matrix),"_"), function(x) x[2])
        
        lineage_table_ccle<-table(lineage_ordered_ccle)
        
        lineage_filtered<-names(which(lineage_table_ccle>10))
        
        #MAKE & POPULATE MATRIX:
        matrix_LINEAGE<-matrix(ncol=0,nrow=nrow(matrix))
        rownames(matrix_LINEAGE)<-rownames(matrix)
        
        for (lineage in lineage_filtered){
                tmp_idx<-which(lineage_ordered_ccle==lineage)
                tmp_vector<-rowMeans(matrix[,tmp_idx],na.rm=T)
                
                matrix_LINEAGE<-cbind(matrix_LINEAGE,tmp_vector[rownames(matrix_LINEAGE)])
        }
        colnames(matrix_LINEAGE)<-lineage_filtered
        
        return(matrix_LINEAGE)
}

make_interLINEAGE_analysis<-function(matrix){
        lineage_ordered_ccle<-sapply(strsplit(colnames(matrix),"_"), function(x) x[2])
        
        lineage_table_ccle<-table(lineage_ordered_ccle)
        
        lineage_filtered<-names(which(lineage_table_ccle>10))
        ########################
        #MASTER-LINEAGE-LIST:
        ########################
        tmp_list<-list()
        
        ########################
        #LINEAGE AVG:
        ########################
        matrix_LINEAGE<-matrix(ncol=0,nrow=nrow(matrix))
        rownames(matrix_LINEAGE)<-rownames(matrix)
        
        for (lineage in lineage_filtered){
                tmp_idx<-which(lineage_ordered_ccle==lineage)
                tmp_vector<-rowMeans(matrix[,tmp_idx],na.rm=T)
                
                matrix_LINEAGE<-cbind(matrix_LINEAGE,tmp_vector[rownames(matrix_LINEAGE)])
        }
        colnames(matrix_LINEAGE)<-lineage_filtered
        
        tmp_list[["avg"]]<-matrix_LINEAGE
        
        ########################
        #LINEAGE STD DEV:
        ########################
        #STD DEV:
        stdev_LINEAGE<-matrix(ncol=length(lineage_filtered),nrow=nrow(matrix))
        rownames(stdev_LINEAGE)<-rownames(matrix)
        colnames(stdev_LINEAGE)<-lineage_filtered
        
        #STD ERROR:
        error_LINEAGE<-matrix(ncol=length(lineage_filtered),nrow=nrow(matrix))
        rownames(error_LINEAGE)<-rownames(matrix)
        colnames(error_LINEAGE)<-lineage_filtered
        
        
        #lineage="brca"
        for (lineage in lineage_filtered){
                lineage_idx<-which(lineage_ordered_ccle==lineage)
                
                #row_name=rownames(stdev_LINEAGE)[1]
                for (row_name in rownames(stdev_LINEAGE)){
                        tmp_dev_vector<-matrix[row_name,lineage_idx]
                        tmp_stdev<-sd(tmp_dev_vector,na.rm=T)
                        tmp_err<-tmp_stdev/sqrt(length(tmp_dev_vector))
                        
                        stdev_LINEAGE[row_name,lineage]<-tmp_stdev
                        error_LINEAGE[row_name,lineage]<-tmp_err
                }
                
        }
        
        tmp_list[["dev"]]<-stdev_LINEAGE
        tmp_list[["error"]]<-error_LINEAGE
        
        ########################
        #LINEAGE STD ERROR::
        ########################
        
        return(tmp_list)
}

color.bar <- function(lut, min, max=-min, ticks=seq(min, max, len=nticks), title='') {
        nticks=(round(max,0)-round(min,0)+1)/2
        scale = (length(lut)-1)/(max-min)
        
        #dev.new(width=1.75, height=5)
        plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
        axis(2, ticks, las=1)
        for (i in 1:(length(lut)-1)) {
                y = (i-1)/scale + min
                rect(0,y,10,y+1/scale, col=lut[i], border=NA)
        }
}

#LOAD DATA:

#RNAseq and RNAseq-SUMS:
load("moa_transform/MODEL-indiv_ccle_tcga.RData")

load("moa_transform/RNAseq-MDR_ccle_tcga.RData")#MADE prev section


#Drug-Sensitivity:
load("cell-line_dbs/AUC_INDIV_zscore.RData")
load("cell-line_dbs/AUC_AVG-zscore.RData")


heatmap_Cluster<-function(matrix,title,row_size=0.8,label="zscore(AUC)"){
        color_gradient<-colorRampPalette(c("blue", "white","red"))(n = 20)
        
        #SORT GENES:
        row_means<-rowMeans(matrix,na.rm=T)#take sum of rows
        row_means_sort<-sort(row_means,decreasing=T)#sort
        
        row_names_sorted<-names(row_means_sort)#extract names
        
        #SORT CANCERS:
        col_means<-colMeans(matrix,na.rm=T)#take sum of columns
        col_means_sort<-sort(col_means,decreasing=T)#sort
        
        col_names_sorted<-names(col_means_sort)
        
        ################
        #SORT PLOT:
        ################
        
        
        heatmap.2(matrix[row_names_sorted,col_names_sorted],col=color_gradient,cexCol =1,cexRow = row_size,keysize = 1,
                  symm=F,density.info = "none",trace="none",margin=c(4,6),main=title,
                  dendrogram = "both",
                  key.xlab = label)
}


#PLOT metab & efflux in black and white:  TCGA as don't have individual cancer drug sensitivity
plot_2D_moa_BW<-function(data=CCLE_efflux_metab,drug="Doxorubicin",type=c("avg","indiv","hybrid"),error=c("se","sd"),hybrid="brca"){
        #LOG TRANSFORM DATA:
        CCLE_efflux_metab_LOG<-log2(data+1)
        CCLE_efflux_metab_LOG_lineage<-make_lineage_avg(CCLE_efflux_metab_LOG)
        
        if(type=="indiv"){
                ########################################################################
                #1  BW INDIVIDUAL
                ########################################################################
                plot(CCLE_efflux_metab_LOG["Doxorubicin_metab",],CCLE_efflux_metab_LOG["Doxorubicin_efflux",],pch=16,cex=0.5)
        }else{
                ########################################################################
                #2  BW AVG:
                ########################################################################
                #CALCULATE AVERAGES & ERRORS
                CCLE_model_ANALYSES<-make_interLINEAGE_analysis(CCLE_efflux_metab_LOG)
                
                #extract drug-data for plotting
                x=CCLE_model_ANALYSES$avg[paste0(drug,"_metab"),]
                y=CCLE_model_ANALYSES$avg[paste0(drug,"_efflux"),]
                if (error=="sd"){
                        x_dev=CCLE_model_ANALYSES$dev[paste0(drug,"_metab"),]
                        y_dev=CCLE_model_ANALYSES$dev[paste0(drug,"_efflux"),]
                }else{
                        x_dev=CCLE_model_ANALYSES$error[paste0(drug,"_metab"),]
                        y_dev=CCLE_model_ANALYSES$error[paste0(drug,"_efflux"),]
                }
                
                
                #PLOT:  data pts, error bard & cancer-labels
                plot(CCLE_efflux_metab_LOG[paste0(drug,"_metab"),],CCLE_efflux_metab_LOG[paste0(drug,"_efflux"),],col="white")
                points(x,y,pch=16,cex=1)
                arrows(x, y-y_dev, x, y+y_dev, length=0.05, angle=90, code=3)
                arrows(x-x_dev, y, x+x_dev, y, length=0.05, angle=90, code=3)
                text(x-0.05,y,labels=names(CCLE_model_ANALYSES$avg[paste0(drug,"_metab"),]))
                
                
                if (type=="hybrid"){
                        type_idx<-grep(paste0("_",hybrid),colnames(CCLE_efflux_metab_LOG))
                        points(CCLE_efflux_metab_LOG[paste0(drug,"_metab"),type_idx],
                               CCLE_efflux_metab_LOG[paste0(drug,"_efflux"),type_idx],col=rgb(0,0,0,0.3),pch=16,cex=0.5)
                        
                }
        }
        
        
        
        
        
        
        
}

#PLOT drug-sensitivity colored 2D plot:
plot_2D_moa_COLOR<-function(data=CCLE_efflux_metab,auc=ctrp_zscore,aucAVG=ctrp_avg_matrix,drug="Doxorubicin",type=c("avg","indiv","hybrid"),error=c("se","sd"),hybrid="brca"){
        #LOG TRANSFORM DATA:
        CCLE_efflux_metab_LOG<-log2(data+1)
        CCLE_efflux_metab_LOG_lineage<-make_lineage_avg(CCLE_efflux_metab_LOG)
        
        
        
        
        if(type=="indiv"){
                ########################################################################
                #1  COLOR-auc INDIVIDUAL
                ########################################################################
                #ID individual lines w/ RNAseq & drug sensitivity
                ccle_lines<-sapply(strsplit(colnames(CCLE_efflux_metab),"_"),function(x) x[1])
                ccle_line_lineages<-sapply(strsplit(colnames(CCLE_efflux_metab),"_"),function(x) x[2])
                
                CCLE_efflux_metab_LOG_names<-CCLE_efflux_metab_LOG
                colnames(CCLE_efflux_metab_LOG_names)<-ccle_lines
                
                #OVERLAPPING LINES:
                lines_overlap<-intersect(colnames(auc),ccle_lines)
                
                #EXRACT DRUGS-sensitivity
                AUC_zscore_na<-auc[tolower(drug),lines_overlap]
                AUC_zscore<-AUC_zscore_na[-which(is.na(AUC_zscore_na)==T)]
                
                #MAKE DRUG-sensitivity COLOR BAR
                rbPal <- colorRampPalette(c('blue','blue','white','red','red'))
                color_scheme_gene <- rbPal(100)[as.numeric(cut(c(AUC_zscore),breaks = 100))]
                
                #PLOT MODEL DATA w/ Drug-sensitivity COLOR BAR
                plot(CCLE_efflux_metab_LOG_names[paste0(drug,"_metab"),names(AUC_zscore)],CCLE_efflux_metab_LOG_names[paste0(drug,"_efflux"),names(AUC_zscore)],pch=16)
                plot(CCLE_efflux_metab_LOG_names[paste0(drug,"_metab"),names(AUC_zscore)],CCLE_efflux_metab_LOG_names[paste0(drug,"_efflux"),names(AUC_zscore)],pch=16,cex=1,
                     col=color_scheme_gene)
                
                
                
                subplot(color.bar(rbPal(100),min(round(AUC_zscore,2))), x=6.4, y=4.5,size=c(0.2,5))
                
        }else{
                ########################################################################
                #2 COLOR-AUC AVG:
                ########################################################################
                CCLE_model_ANALYSES<-make_interLINEAGE_analysis(CCLE_efflux_metab_LOG)
                
                #OVERLAPPING LINES:
                lineage_overlap<-intersect(colnames(aucAVG),colnames(CCLE_efflux_metab_LOG_lineage))
                
                #EXRACT DRUGS-sensitivity
                AUC_zscore_avg<-aucAVG[tolower(drug),lineage_overlap]
                
                
                #MAKE DRUG-sensitivity COLOR BAR
                rbPal <- colorRampPalette(c('darkblue','blue','white','red','darkred'))
                color_scheme_avg <- rbPal(100)[as.numeric(cut(c(AUC_zscore_avg),breaks = 100))]
                
                #PLOT MODEL DATA w/ Drug-sensitivity COLOR BAR
                x=CCLE_model_ANALYSES$avg[paste0(drug,"_metab"),lineage_overlap]
                y=CCLE_model_ANALYSES$avg[paste0(drug,"_efflux"),lineage_overlap]
                if (error=="sd"){
                        x_dev=CCLE_model_ANALYSES$dev[paste0(drug,"_metab"),lineage_overlap]
                        y_dev=CCLE_model_ANALYSES$dev[paste0(drug,"_efflux"),lineage_overlap]
                }else{
                        x_dev=CCLE_model_ANALYSES$error[paste0(drug,"_metab"),lineage_overlap]
                        y_dev=CCLE_model_ANALYSES$error[paste0(drug,"_efflux"),lineage_overlap]
                }
                
                
                #MAKE PLOT OF ALL INDIV DATA:
                plot(CCLE_efflux_metab_LOG[paste0(drug,"_metab"),],CCLE_efflux_metab_LOG[paste0(drug,"_efflux"),],col="white")
                points(x,y,pch=16,cex=1,col=color_scheme_avg)
                arrows(x, y-y_dev, x, y+y_dev, length=0.05, angle=90, code=3,col=color_scheme_avg)
                arrows(x-x_dev, y, x+x_dev, y, length=0.05, angle=90, code=3,col=color_scheme_avg)
                text(0.05+x,y,labels=lineage_overlap,col=color_scheme_avg)
                
                
                subplot(color.bar(rbPal(100),min(round(AUC_zscore_avg,2))), x=6.4, y=4.5,size=c(0.2,5))
                
                
                
                ##################
                #ADD ONE LINEAGE:
                #################
                
                
                if (type=="hybrid"){
                        type_idx<-grep(paste0("_",hybrid),colnames(CCLE_efflux_metab_LOG))
                        points(CCLE_efflux_metab_LOG[paste0(drug,"_metab"),type_idx],
                               CCLE_efflux_metab_LOG[paste0(drug,"_efflux"),type_idx],col=rgb(0,0,0,0.3),pch=16,cex=0.5)
                        
                }
        }
        
        
        
        
        
        
        
}




################################################################################################################################################
#Figure S11:  Heatmaps of metabolism & efflux transformations of 
################################################################################################################################################




heatmap_Cluster(CCLE_efflux_metab_zscore_lineage,"CCLE efflux/metab Model Heatmap",row_size=0.6)

heatmap_Cluster(TCGA_efflux_metab_zscore_lineage,"TCGA efflux/metab Model Heatmap",row_size=0.6)






################################################################################################################################################
#Figure 3A-B:  DOXORUBUCIN CCLE data color-coded by average drug-AUC:
################################################################################################################################################

#DOXORUBICIN:
plot_2D_moa_COLOR(data=CCLE_efflux_metab,
                  auc=ctrp_zscore,
                  aucAVG=ctrp_avg_matrix,
                  drug="Doxorubicin",
                  type="hybrid",
                  error="se",
                  hybrid="brca")

##################################################################################################################################################################################
#Figure 3C-D:  DOXORUBICIN TCGA data 
##################################################################################################################################################################################

#CCLE:  avg
plot_2D_moa_BW(CCLE_efflux_metab,drug="Doxorubicin",type="hybrid",error="se",hybrid="brca")


#TCGA:  avg
plot_2D_moa_BW(TCGA_efflux_metab,drug="Doxorubicin",type="hybrid",error="se",hybrid="brca")




####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##         ##     ##                                                          
##     ##          ##    ##                          ##   ##         ##     ##                                      
##     ##          ##   ##                      #################    ##     ##                                         
##     ########    ##   ##   ####                    ##   ##         ###########                                                 
##     ##          ##   ##      ##              #################           ##                                               
##     ##          ##    ##     ##                  ##   ##                 ##                                     
##     ##          ##     #######                   ##   ##                 ##                                             
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########   ####    ##     ########    #######             ###########  ##    ##   #######    #######                                             
##     ##      ##  ## ##   ##    ##          ##     ##                ##        ##  ##    ##    ##   ##                         
##     ##      ##  ##  ##  ##   ##           ##     ##                ##         ####     ##    ##   ##                                
##     ##      ##  ##   ## ##   ##           ##     ##   #######      ##          ##      #######    #######                                    
##     ##      ##  ##    ####   ##           ##     ##                ##          ##      ##         ##                        
##     ##      ##  ##     ###    ##          ##     ##                ##          ##      ##         ##                      
##      ########   ##      ##     ########    #######                 ##          ##      ##         #######                                          
####################################################################################################################################################################################################################################################################
#################################################################################################################################################################################################################################################################### 
##      ########    #######     ########      #######    ##              ##      ########### ########    ########   ####    ##                                                                            
##     ##          ##     ##    ##     ##     ##    ##   ##             ####         ##         ##      ##      ##  ## ##   ##                
##    ##           ##     ##    ##     ##     ##    ##   ##            ##  ##        ##         ##      ##      ##  ##  ##  ##                         
##    ##           ##     ##    ########      #######    ##           ##    ##       ##         ##      ##      ##  ##   ## ##                 
##    ##           ##     ##    ##    ##      ##   ##    ##          ##########      ##         ##      ##      ##  ##    ####                    
##     ##          ##     ##    ##     ##     ##    ##   ##         ##        ##     ##         ##      ##      ##  ##     ###                  
##      ########    #######     ##      ##    ##     ##  ########  ##          ##    ##      ########    ########   ##      ##                             
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################



#######################################
#SIGNATURES:
#######################################

oncotypeDX<-c("BIRC5","AURKA","MKI67","TFRC","MYBL2","CCNB1","GAPDH","RPLP0","ACTB","GUSB","BAG1","MMP11","BCL2","GRB7","GSTM1","ERBB2","SCUBE2","ESR1","PGR")


mammaprint<-c("CDCA7","DIAPH3","NMU","TSPYL5","MELK","NDC80","MMP9","CENPA","TMEM65","HRASLS","ORC6","PRC1","ESM1","RFC4","CCNE2","ECT2","NUSAP1","MCM6","MTDH",
              "DTL","GMPS","MSANTD3","UCHL5","QSOX2","EXT1","COL4A2","GPR180","DCK","CMC2","PITRM1","OXCT1","EGLN1","LPCAT1","GNAZ","FBXO31","RAB6A","FLT1",
              "SLC2A3","CDC42BPA","DHX58","WISP1","TMEM74B","FGF18","IGFBP5","AP2B1","ECI2","RUNDC1","BBC3","TGFB3","RTN4RL1","EBF4","ALDH4A1","SMIM5","GSTM3",
              "STK32B","ZNF385B","MS4A7","SCUBE2")



pam50<-c("FOXC1","ANLN","SFRP1","CCNE1","PTTG1","CDC20","EXO1","MELK","NDC80","PHGDH","UBE2C","KIF2C","CEP55","BIRC5","ORC6","RRM2","CENPF","CDH3","NUF2","ACTR3B",
         "MYC","KRT5","MKI67","MYBL2","EGFR","UBE2T","CCNB1","CDC6","TYMS","MDM2","BAG1","MMP11","BCL2","FGFR4","BLVRA","GRB7","GPR160","CXXC5","ERBB2","SLC39A6",
         "TMEM45B","MLPH","MAPT","FOXA1","ESR1","NAT1","PGR")

brca_sigs_union<-unique(c(oncotypeDX,mammaprint,pam50))
brca_sigs_intersect<-sort(table(c(oncotypeDX,mammaprint,pam50)),decreasing=T)


#######################################
#PLOTTING FUNCTIONS:
#######################################

library("gplots")
library("RColorBrewer")

#4color gradient
hcluster_unsort<-function(matrix,type,col_size=1,row_size=0.4,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=NA,Colv=NA, cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white","red"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab = units,symm=F,dendrogram = "none")
}

#CLUSTERED:
hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
        #distance matrix:
        dist_columns <- dist(t(matrix))
        dist_rows <- dist(matrix)
        
        #hclustering:
        hclust_columns<-hclust(dist_columns, method="average")
        hclust_rows<-hclust(dist_rows, method="average")
        
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), 
                  cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "lightcyan","red"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab =units)
}



hcluster_corr<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
        #distance matrix:
        cor_columns <- cor(matrix,method="spearman")
        cor_rows <- cor(t(matrix),method="spearman")
        
        
        #distance matrix:
        dist_columns <- as.dist(1-cor_columns)
        dist_rows <- as.dist(1-cor_rows)
        
        #hclustering:
        hclust_columns<-hclust(dist_columns, method="average")
        hclust_rows<-hclust(dist_rows, method="average")
        
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}


#######################################
#DATA
#######################################



load("cell-line_dbs/invitro_AUC-lines-individ_dbs.RData")
rownames(gdsc1_auc_matrix)<-tolower(rownames(gdsc1_auc_matrix))
rownames(gdsc2_auc_matrix)<-tolower(rownames(gdsc2_auc_matrix))

load("cell-line_dbs/CCLE_RNAseq.RData")
load("expr_literature/MASTER_DrugBankPittQM_v5_sider.RData")
rownames(MASTER_DrugBankPittQM_v5_sider)<-tolower(rownames(MASTER_DrugBankPittQM_v5_sider))


#######################################
#PLOTS:
#######################################



#1  CANCER-SCREEN DRUG OVERLAP:  PRISM (136), CTRP (58), GDSC1 (50), GDSC2 (45)
cancer_drugs_drugBank<-names(which((MASTER_DrugBankPittQM_v5_sider[,"disease"]=="ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS") & 
                                           (MASTER_DrugBankPittQM_v5_sider[,"type"]=="small molecule")))

prism_drugBank_overlap_idx<-which(cancer_drugs_drugBank %in% rownames(prism_auc_matrix))
cancer_drugs_drugBank_prism<-cancer_drugs_drugBank[prism_drugBank_overlap_idx]

prism_auc_cancer<-prism_auc_matrix[cancer_drugs_drugBank_prism,]#PRISM:  136 drugs
ctrp_auc_cancer<-ctrp_auc_matrix[which(rownames(ctrp_auc_matrix) %in% cancer_drugs_drugBank_prism),]#CTRP:  58 drugs
gdsc1_auc_cancer<-gdsc1_auc_matrix[which(rownames(gdsc1_auc_matrix) %in% cancer_drugs_drugBank_prism),]#GDSC1:  50 drugs
gdsc2_auc_cancer<-gdsc2_auc_matrix[which(rownames(gdsc2_auc_matrix) %in% cancer_drugs_drugBank_prism),]#GDSC2:  45 drugs


#2 BREAST CANCER OVERLAP:  GDSC(45), CTRP (40), PRISM (22)
breast_lines<-rownames(ccle_lineage_filtered)[which(ccle_lineage_filtered$lineage=="breast")]

overlap_lines_prism<-colnames(prism_auc_cancer)[which(colnames(prism_auc_cancer) %in% breast_lines)]
overlap_lines_ctrp<-colnames(ctrp_auc_cancer)[which(colnames(ctrp_auc_cancer) %in% breast_lines)]
overlap_lines_gdsc1<-colnames(gdsc1_auc_cancer)[which(colnames(gdsc1_auc_cancer) %in% breast_lines)]
overlap_lines_gdsc2<-colnames(gdsc2_auc_cancer)[which(colnames(gdsc2_auc_cancer) %in% breast_lines)]


#3 CORRELATION MATRICES:

brca_prism_RNA_corr<-cor(-t(prism_auc_cancer[,overlap_lines_prism]),
                         t(ccle_RNAseq_filtered[brca_sigs_union,overlap_lines_prism]),use= "pairwise.complete.obs")

brca_ctrp_RNA_corr<-cor(-t(ctrp_auc_cancer[,overlap_lines_ctrp]),
                        t(ccle_RNAseq_filtered[brca_sigs_union,overlap_lines_ctrp]),use= "pairwise.complete.obs")


brca_gdsc1_RNA_corr<-cor(-t(gdsc1_auc_cancer[,overlap_lines_gdsc1]),
                         t(ccle_RNAseq_filtered[brca_sigs_union,overlap_lines_gdsc1]),use= "pairwise.complete.obs")

brca_gdsc2_RNA_corr<-cor(-t(gdsc2_auc_cancer[,overlap_lines_gdsc2]),
                         t(ccle_RNAseq_filtered[brca_sigs_union,overlap_lines_gdsc2]),use= "pairwise.complete.obs")
brca_gdsc_RNA_corr<-rbind(brca_gdsc1_RNA_corr[,brca_sigs_union],brca_gdsc2_RNA_corr[,brca_sigs_union])


hcluster_corr(brca_prism_RNA_corr,row_size = 0.2,col_size = 0.2,title = "PRISM-ONCOTYPE CORR:")

hcluster_corr(brca_ctrp_RNA_corr,row_size = 0.5,col_size = 0.2,title = "CTRP-ONCOTYPE CORR:")

hcluster_corr(brca_gdsc_RNA_corr,row_size = 0.2,col_size = 0.2,title = "GDSC-ONCOTYPE CORR:")



####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##         ##     ##                                                          
##     ##          ##    ##                          ##   ##         ##     ##                                      
##     ##          ##   ##                      #################    ##     ##                                         
##     ########    ##   ##   ####                    ##   ##         ###########                                                 
##     ##          ##   ##      ##              #################           ##                                               
##     ##          ##    ##     ##                  ##   ##                 ##                                     
##     ##          ##     #######                   ##   ##                 ##                                             
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ########     ########    ##     ##     ########             ##        ##  #########    #######   ##     ##                                                                                                               
##    ##     ##    ##     ##   ##     ##    ##                    ###      ###  ##          ##         ##     ##                                     
##    ##      ##   ##     ##   ##     ##   ##                     ####    ####  ##         ##          ##     ##                                                         
##    ##      ##   ########    ##     ##   ##    ####    ######   ## ##  ## ##  ########   ##          #########                  
##    ##      ##   ##    ##    ##     ##   ##       ##            ##  ####  ##  ##         ##          ##     ##                                                
##    ##     ##    ##     ##   ##     ##    ##      ##            ##   ##   ##  ##          ##         ##     ##                                                       
##    ########     ##      ##   #######      ########             ##        ##  #########    #######   ##     ##                                                                 
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##    ##      ##     ######      ##          ##        ###      ###    #######     #######     ##     ##     #######                                                                                                                                  
##    ##      ##    ##    ##     ##          ##        ####    ####   ##     ##    ##    ##    ##    ##     ##                                                                                      
##    ##      ##   ##      ##    ##          ##        ## ##  ## ##  ##       ##   ##    ##    ##   ##      ##                                                                                                           
##    ##########   ##########    ##          ##        ##  ####  ##  ###########   #######     ######        #######                                                                                                 
##    ##      ##   ##      ##    ##          ##        ##   ##   ##  ##       ##   ##    ##    ##   ##             ##                                                                              
##    ##      ##   ##      ##    ##          ##        ##        ##  ##       ##   ##     ##   ##    ##            ##                                                                            
##    ##      ##   ##      ##    #########   ########  ##        ##  ##       ##   ##      ##  ##     ##    ########                                                                                                   
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################






library("gplots")
library("RColorBrewer")

#4color gradient
hcluster_unsort<-function(matrix,type,col_size=1,row_size=0.4,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), 
                  Rowv=NA,Colv=NA, cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white","red"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab = units,symm=F,dendrogram = "none")
}

#CLUSTERED:
hcluster<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
        #distance matrix:
        dist_columns <- dist(t(matrix))
        dist_rows <- dist(matrix)
        
        #hclustering:
        hclust_columns<-hclust(dist_columns, method="average")
        hclust_rows<-hclust(dist_rows, method="average")
        
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), 
                  cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "lightcyan","red"))(n = 20),
                  density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab =units)
}

hcluster_corr<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering"){
        #distance matrix:
        cor_columns <- cor(matrix,method="spearman")
        cor_rows <- cor(t(matrix),method="spearman")
        
        
        #distance matrix:
        dist_columns <- as.dist(1-cor_columns)
        dist_rows <- as.dist(1-cor_rows)
        
        #hclustering:
        hclust_columns<-hclust(dist_columns, method="average")
        hclust_rows<-hclust(dist_rows, method="average")
        
        
        #PLOT HEATMAP:
        heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white", "red"))(n = 20),density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = "NES",key.xlab = "NES")
}

#INVITRO SENSITIVITY CORRELATIONS
lineage_biomarkers_SIGNATURE<-function(drug="doxorubicin",cancer="breast",gene_num=10,signature,
                                       row_size=0.5,col_size=1){
        #################################
        #1.1 EXTRACT & NORMALIZE biomarkers w/r to CCLE-reference
        #################################
        
        #BRCA LINES:
        brca_lines<-colnames(ccle_brca_sigs_zscore)
        drug_lines<-unique(c(colnames(ctrp_zscore),colnames(gdsc1_zscore),colnames(prism_zscore),colnames(gdsc2_zscore)))
        brca_drug_lines<-brca_lines[which(brca_lines %in% drug_lines)]
        
        #BIOMARKERS
        drug_resist_genes_raw<-signature
        drug_resist_genes<-drug_resist_genes_raw[which(drug_resist_genes_raw %in% rownames(ccle_brca_sigs_zscore))]
        
        
        #CCLE z-score:
        brca_RNAseq_CCLEzscore<-t(ccle_brca_sigs_zscore[drug_resist_genes,])
        
        #hcluster(brca_RNAseq_CCLEzscore,col_size = 0.7,row_size = 0.7,title=paste0(drug," biomarkers:  CCLE-norm"))
        
        
        #################################
        #1.2 EXTRACT & NORMALIZE drug sensitivity
        #################################
        #ID Datasets that have drug:
        all_datasets<-c("CTRP","PRISM","GDSC1","GDSC2")
        
        drug_idx<-c()
        if(drug %in% rownames(ctrp_zscore)){drug_idx<-c(drug_idx,1)}
        if(drug %in% rownames(prism_zscore)){drug_idx<-c(drug_idx,2)}
        if(drug %in% rownames(gdsc1_zscore)){drug_idx<-c(drug_idx,3)}
        if(drug %in% rownames(gdsc2_zscore)){drug_idx<-c(drug_idx,4)}
        
        drug_datasets<-all_datasets[drug_idx]
        
        #POPULATION AUC MATRIX:
        auc_matrix<-matrix(ncol=length(drug_datasets),nrow=length(brca_drug_lines))
        colnames(auc_matrix)<-drug_datasets
        rownames(auc_matrix)<-brca_drug_lines
        
        #line="DU4475"
        for (line in brca_drug_lines){
                if ("CTRP" %in% drug_datasets){ 
                        if (line %in% colnames(ctrp_zscore)) {auc_matrix[line,"CTRP"]<-ctrp_zscore[drug,line]} }
                if ("PRISM" %in% drug_datasets){ 
                        if (line %in% colnames(prism_zscore)) {auc_matrix[line,"PRISM"]<-prism_zscore[drug,line]} }
                if ("GDSC1" %in% drug_datasets){ 
                        if (line %in% colnames(gdsc1_zscore)) {auc_matrix[line,"GDSC1"]<-gdsc1_zscore[drug,line]} }
                if ("GDSC2" %in% drug_datasets){ 
                        if (line %in% colnames(gdsc2_zscore)) {auc_matrix[line,"GDSC2"]<-gdsc2_zscore[drug,line]} }
                
        }
        
        auc_matrix_mean<-cbind(auc_matrix,"AVG"=rowMeans(auc_matrix,na.rm=T))
        
        
        
        rowSort_auc<-names(sort(auc_matrix_mean[,"AVG"],decreasing=T))
        #par(las=1)
        #hcluster_unsort(auc_matrix_mean[rowSort_auc,],col_size = 0.7,row_size = 0.7,title=paste0(cancer," cancer ",drug," AUC's"))
        
        
        
        #################################
        #1.3 BIOMARKER CORRELATION
        #################################
        line_overlap<-intersect(rownames(auc_matrix_mean),rownames(brca_RNAseq_CCLEzscore))
        drug_AUC_multiLinReg<-cbind(auc_matrix_mean[line_overlap,],brca_RNAseq_CCLEzscore[line_overlap,])
        
        
        gene_num_final<-min(gene_num,length(drug_resist_genes))
        
        #SORT:
        colSort_corr_raw<-names(sort(cor(drug_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],decreasing=T))
        
        abs_corr_top<-names(sort(abs(cor(drug_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",]),decreasing=T))[1:(gene_num_final+length(drug_datasets))]
        
        colSort_corr<-colSort_corr_raw[which(colSort_corr_raw %in% abs_corr_top)]
        
        
        #PLOT:
        drug_AUCgene_rank<-cor(drug_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",][colSort_corr[length(colSort_corr):1]]
        
        
        #par(las=2)
        #barplot(drug_AUCgene_rank,horiz = T,cex.names=col_size,main=paste0(drug,": biomarker-AUC corr in ",cancer))
        
        
        
        #ADD CORR TO COL-NAMES:
        name2corr<-round(cor(drug_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],2)
        
        sorted_matrix<-drug_AUC_multiLinReg[rowSort_auc,colSort_corr]
        colnames(sorted_matrix)<-paste0(colnames(sorted_matrix)," (",name2corr[colnames(sorted_matrix)],")")
        
        
        
        par(las=1)
        hcluster_unsort(sorted_matrix,col_size = col_size,row_size = row_size,title=paste0(drug," biomarkers(",cancer,")"))
        
        return(drug_AUC_multiLinReg[rowSort_auc,colSort_corr])
}


#CLINICAL SENSITIVITY CORRELATIONS:
lineage_biomarkers_CLINIC_sig<-function(signature,dataset,metadata,row_size=0.2,col_size=1,gene_num="all"){
        
        
        #################################
        #1.1 EXTRACT & NORMALIZE biomarkers w/r to CCLE-reference
        #################################
        
        #BIOMARKERS
        drug_genes_overlap<-signature[which(signature %in% rownames(dataset))]
        
        #EXTRACT RNAseq
        brca_RNAseq<-t(dataset[drug_genes_overlap,])
        
        
        #INTERNAL z-score
        brca_RNAseq_INTRNLzscore<-scale(brca_RNAseq,center=T,scale=T)
        #hcluster(brca_RNAseq_INTRNLzscore,col_size = 0.7,row_size = 0.7)
        
        #CCLE z-score:   NEED TO CHANGE w/ TCGA
        #brca_RNAseq_CCLEzscore<-scale(brca_RNAseq,
        #                              center=apply(ccle_RNAseq_filtered[colnames(brca_RNAseq),],1,mean),
        #                              scale=apply(ccle_RNAseq_filtered[colnames(brca_RNAseq),],1,sd))
        
        #hcluster(brca_RNAseq_CCLEzscore,col_size = 0.7,row_size = 0.7,title=paste0(drug," biomarkers:  CCLE-norm"))
        
        
        #################################
        #1.2 EXTRACT & NORMALIZE drug sensitivity
        #################################
        
        all_patients<-as.character(metadata[,"sample"])
        responders_idx<-which(metadata[,"other"]=="postchemo: no invasive cancer")
        
        response_vector<-rep(0,length(all_patients))
        names(response_vector)<-all_patients
        response_vector[responders_idx]<-1
        
        response_vector_sort<-sort(response_vector,decreasing=T)
        
        #################################
        #1.3 BIOMARKER CORRELATION
        #################################
        line_overlap<-intersect(names(response_vector_sort),rownames(brca_RNAseq_INTRNLzscore))
        dox_AUC_multiLinReg<-cbind("AVG"=response_vector_sort,brca_RNAseq_INTRNLzscore[names(response_vector_sort),])
        
        
        #RANK OR SELECT TOP& BOTTOM BIOMAKRES
        if (gene_num=="all"){
                #SORT:
                colSort_corr<-names(sort(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],decreasing=T))
                
                #PLOT:
                dox_AUCgene_rank<-sort(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],decreasing=F)
                
                #par(las=2)
                #barplot(dox_AUCgene_rank,horiz = T,cex.names=0.3,main="Biomarkers Correlation")
                
                
        }else{
                #SORT:
                colSort_corr_raw<-names(sort(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],decreasing=T))
                
                abs_corr_top<-names(sort(abs(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",]),decreasing=T))[1:(gene_num+3)]
                
                colSort_corr<-colSort_corr_raw[which(colSort_corr_raw %in% abs_corr_top)]
                
                
                #PLOT:
                dox_AUCgene_rank<-cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",][colSort_corr[length(colSort_corr):1]]
                
                
                #par(las=2)
                #barplot(dox_AUCgene_rank,horiz = T,cex.names=col_size,main="Biomarkers Correlation")
                
                
                
                
        }
        
        
        
        
        
        
        #################################
        #1.4 Cluster N
        #################################
        resist_data<-dox_AUC_multiLinReg[which(dox_AUC_multiLinReg[,"AVG"]==1),colSort_corr]
        sense_data<-dox_AUC_multiLinReg[which(dox_AUC_multiLinReg[,"AVG"]==0),colSort_corr]
        
        resist_cor_rows <- cor(t(resist_data),method="spearman")
        sense_cor_rows <- cor(t(sense_data),method="spearman")
        
        
        #distance matrix:
        dist_resist <- as.dist(1-resist_cor_rows)
        dist_sense <- as.dist(1-sense_cor_rows)
        
        #hclustering:
        hclust_resist<-hclust(dist_resist, method="average")
        hclust_sense<-hclust(dist_sense, method="average")
        
        
        #order:
        hclust_resist<-rownames(resist_data)[hclust_resist$order]
        hclust_sense<-rownames(sense_data)[hclust_sense$order]
        
        
        
        
        
        #ADD CORR TO COL-NAMES:
        name2corr<-round(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],2)
        
        sorted_matrix<-dox_AUC_multiLinReg[c(hclust_resist,hclust_sense),colSort_corr]
        colnames(sorted_matrix)<-paste0(colnames(sorted_matrix)," (",name2corr[colnames(sorted_matrix)],")")
        
        
        
        par(las=1)
        hcluster_unsort(sorted_matrix,col_size = 0.7,row_size = row_size,title=paste0("Signature biomarkers"))
        
        return(dox_AUC_multiLinReg[c(hclust_resist,hclust_sense),colSort_corr])
}
lineage_biomarkers_CLINIC<-function(drug_set="doxorubicin",dataset,metadata,row_size=0.2,col_size=1,gene_num="all"){
        
        
        tmp_biomarkers<-chemo_BIOMARKERS_matrix_DDR
        rownames(tmp_biomarkers)<-tolower(rownames(tmp_biomarkers))
        
        #################################
        #1.1 EXTRACT & NORMALIZE biomarkers w/r to CCLE-reference
        #################################
        
        #BIOMARKERS
        if (length(drug_set)==1){
                drug=drug_set
                
                drug_genes<-unique(names(which(tmp_biomarkers[drug,]==1)))
                drug_genes_overlap<-drug_genes[which(drug_genes %in% rownames(dataset))]
        }else{
                drug_genes_overlap<-c()
                for (drug in drug_set){
                        tmp_genes<-unique(names(which(tmp_biomarkers[drug,]==1)))
                        tmp_genes_overlap<-tmp_genes[which(tmp_genes %in% rownames(dataset))]
                        
                        drug_genes_overlap<-unique(c(drug_genes_overlap,tmp_genes_overlap))
                }
                
        }
        
        
        #EXTRACT RNAseq
        brca_RNAseq<-t(dataset[drug_genes_overlap,])
        
        
        #INTERNAL z-score
        brca_RNAseq_INTRNLzscore<-scale(brca_RNAseq,center=T,scale=T)
        #hcluster(brca_RNAseq_INTRNLzscore,col_size = 0.7,row_size = 0.7)
        
        #CCLE z-score:   NEED TO CHANGE w/ TCGA
        #brca_RNAseq_CCLEzscore<-scale(brca_RNAseq,
        #                              center=apply(ccle_RNAseq_filtered[colnames(brca_RNAseq),],1,mean),
        #                              scale=apply(ccle_RNAseq_filtered[colnames(brca_RNAseq),],1,sd))
        
        #hcluster(brca_RNAseq_CCLEzscore,col_size = 0.7,row_size = 0.7,title=paste0(drug," biomarkers:  CCLE-norm"))
        
        
        #################################
        #1.2 EXTRACT & NORMALIZE drug sensitivity
        #################################
        
        all_patients<-as.character(metadata[,"sample"])
        responders_idx<-which(metadata[,"other"]=="postchemo: no invasive cancer")
        
        response_vector<-rep(0,length(all_patients))
        names(response_vector)<-all_patients
        response_vector[responders_idx]<-1
        
        response_vector_sort<-sort(response_vector,decreasing=T)
        
        #################################
        #1.3 BIOMARKER CORRELATION
        #################################
        line_overlap<-intersect(names(response_vector_sort),rownames(brca_RNAseq_INTRNLzscore))
        dox_AUC_multiLinReg<-cbind("AVG"=response_vector_sort,brca_RNAseq_INTRNLzscore[names(response_vector_sort),])
        
        
        #RANK OR SELECT TOP& BOTTOM BIOMAKRES
        if (gene_num=="all"){
                #SORT:
                colSort_corr<-names(sort(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],decreasing=T))
                
                #PLOT:
                dox_AUCgene_rank<-sort(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],decreasing=F)
                
                #par(las=2)
                #barplot(dox_AUCgene_rank,horiz = T,cex.names=0.3,main="Biomarkers Correlation")
                
                
        }else{
                #SORT:
                colSort_corr_raw<-names(sort(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",],decreasing=T))
                
                abs_corr_top<-names(sort(abs(cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",]),decreasing=T))[1:(gene_num+3)]
                
                colSort_corr<-colSort_corr_raw[which(colSort_corr_raw %in% abs_corr_top)]
                
                
                #PLOT:
                dox_AUCgene_rank<-cor(dox_AUC_multiLinReg,use = "pairwise.complete.obs")["AVG",][colSort_corr[length(colSort_corr):1]]
                
                
                #par(las=2)
                #barplot(dox_AUCgene_rank,horiz = T,cex.names=col_size,main="Biomarkers Correlation")
                
                
                
                
        }
        
        #################################
        #1.4 Cluster N
        #################################
        resist_data<-dox_AUC_multiLinReg[which(dox_AUC_multiLinReg[,"AVG"]==1),colSort_corr]
        sense_data<-dox_AUC_multiLinReg[which(dox_AUC_multiLinReg[,"AVG"]==0),colSort_corr]
        
        resist_cor_rows <- cor(t(resist_data),method="spearman")
        sense_cor_rows <- cor(t(sense_data),method="spearman")
        
        
        #distance matrix:
        dist_resist <- as.dist(1-resist_cor_rows)
        dist_sense <- as.dist(1-sense_cor_rows)
        
        #hclustering:
        hclust_resist<-hclust(dist_resist, method="average")
        hclust_sense<-hclust(dist_sense, method="average")
        
        
        #order:
        hclust_resist<-rownames(resist_data)[hclust_resist$order]
        hclust_sense<-rownames(sense_data)[hclust_sense$order]
        
        
        
        par(las=1)
        hcluster_unsort(dox_AUC_multiLinReg[c(hclust_resist,hclust_sense),colSort_corr],col_size = 0.7,row_size = row_size,title=paste0(drug," biomarkers"))
        
        return(dox_AUC_multiLinReg[c(hclust_resist,hclust_sense),colSort_corr])
}

#######################
#Cell-line data: RNA & drug-sensitivity
#######################
load("cell-line_dbs/CCLE_RNAseq.RData")

load("cell-line_dbs/AUC_INDIV_zscore.RData")


#######################
#MDACC Patient data: 
#######################
load("clinical_dbs/MDACC_IGR_brca_clinical.RData")

load("clinical_dbs/USO_02103_brca_clinical.RData")




#######################
#SIGNATURE TRANSFORMED:   tcga, ccle, mdacc data (made in #0 and #1 below)
#######################
load("clinical_dbs/DRUG_HALLMARK-transformed-RNAseq.RData")

load("expr_literature/BIOMARKER-signs_matrices.RData")
rm(chemo_BIOMARKERS_matrix)
rm(chemo_ddr_matrix)
rm(chemo_transport_matrix)
rm(chemo_metab_matrix)
rm(chemo_target_matrix)
rm(chemo_BIOMARKERS_matrix_noTargets)



drug_signature<-unique(c(rownames(igr_brca_sigs_zscore)[grep("Doxorubicin",rownames(igr_brca_sigs_zscore))],
                         rownames(igr_brca_sigs_zscore)[grep("Cyclophosphamide",rownames(igr_brca_sigs_zscore))],
                         rownames(igr_brca_sigs_zscore)[grep("Fluorouracil",rownames(igr_brca_sigs_zscore))],
                         "MITOTIC_SPINDLE","ESTROGEN_RESPONSE_EARLY","ESTROGEN_RESPONSE_LATE","ESTROGEN_RESPONSE_EARLY","G2M_CHECKPOINT","E2F_TARGETS","DNA_REPAIR","XENOBIOTIC_METABOLISM","XENOBIOTIC_METABOLISM","APOPTOSIS","REACTIVE_OXIGEN_SPECIES_PATHWAY"))


#################################################################################################################################################################
#FIGURE 4C:  Proliferation correlations IN VITRO
#################################################################################################################################################################
CCLE_assocations<-sort(cor(t(ccle_brca_sigs_zscore))["E2F_TARGETS",drug_signature])

barplot(CCLE_assocations,horiz=T)


#################################################################################################################################################################
#FIGURE 4D:  DOXORUBICIN correlations across 49 breast cancers
#################################################################################################################################################################
otdx_glm<-lineage_biomarkers_SIGNATURE(drug="doxorubicin",cancer="breast",row_size = 0.3,gene_num = 20,
                                       signature=drug_signature,col_size = 0.8)





#################################################################################################################################################################
#FIGURE 4E:  Proliferation correlations Clincial
#################################################################################################################################################################

TCGA_assocations<-sort(cor(t(tcga_brca_sigs_zscore))["E2F_TARGETS",drug_signature])

barplot(TCGA_assocations,horiz=T)



#################################################################################################################################################################
#FIGURE 4F & S12A:  Neoadjuvant CAF corrleations:   
#################################################################################################################################################################

OTDX_glm<-lineage_biomarkers_CLINIC_sig(signature=drug_signature,
                                        dataset=igr_brca_sigs_zscore,
                                        metadata=MDACC_IGR_metadata,gene_num = 20)
#################################################################################################################################################################
#FIGURE S12B:  Neoadjuvant CAF corrleations:   
#################################################################################################################################################################

OTDX_glm<-lineage_biomarkers_CLINIC_sig(signature=drug_signature,
                                        dataset=uso_brca_sigs_zscore,
                                        metadata=USO_02103_metadata,gene_num = 20)














####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##     ########    ##     #######                    ##   ##            ##########                                                      
##     ##          ##    ##                          ##   ##            ##                                 
##     ##          ##   ##                      #################       ##                                      
##     ########    ##   ##   ####                    ##   ##             ######                                             
##     ##          ##   ##      ##              #################             ##                                              
##     ##          ##    ##     ##                  ##   ##                   ##                                    
##     ##          ##     #######                   ##   ##             #######                                           
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##      #######      ##      #######    #######    ##  #######   #######                                                                                    
##      ##    ##    ####     ##    ##   ##    ##   ##  ##        ##    ##                                                                                      
##      ##    ##   ##  ##    ##    ##   ##    ##   ##  ##        ##    ##                                                                                                           
##      #######   ##    ##   #######    #######    ##  #######   #######                                                                                      
##      ##    ##  ########   ##   ##    ##   ##    ##  ##        ##   ##                                                                                                       
##      ##    ##  ##    ##   ##    ##   ##    ##   ##  ##        ##    ##                                                                                                         
##      #######   ##    ##   ##     ##  ##     ##  ##  #######   ##     ##                                                                                                                        
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##        ##          ##    #######    ########     #########   ##         #########                                                              
##        ###        ###   ##     ##   ##     ##    ##          ##        ##                                                  
##        ####      ####   ##     ##   ##      ##   ##          ##        ##                                                  
##        ## ##    ## ##   ##     ##   ##      ##   #######     ##         #######                                                       
##        ##  ##  ##  ##   ##     ##   ##      ##   ##          ##               ##                                            
##        ##   ####   ##   ##     ##   ##     ##    ##          ##               ##                                            
##        ##    ##    ##    #######    ########     #########   ########   #######                                                           
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  



library("gplots")
library("RColorBrewer")


load("expr_literature/BIOMARKER-signs_matrices.RData")
rm(chemo_BIOMARKERS_matrix)
rm(chemo_BIOMARKERS_matrix_DDR)
rm(chemo_BIOMARKERS_matrix_noTargets)
rm(chemo_ddr_matrix)
rm(chemo_target_matrix)
chemo_metab_matrix<-(-chemo_metab_matrix)
chemo_transport_matrix<-(-chemo_transport_matrix)

################
#FUNCTIONS:
################

#MAKE MoA matrix:   efflux & metab
make_MoA_model<-function(matrix,transport,metab){
  #1 convert log2 to tpm:
  tmp_matrix<-(2^(matrix)-1)
  
  
  #2 POPULATE METAB:
  overlap_metab_genes<-intersect(colnames(metab),rownames(matrix))
  tmp_metab_model<-metab[,overlap_metab_genes] %*% matrix[overlap_metab_genes,]
  rownames(tmp_metab_model)<-paste(rownames(tmp_metab_model),"metab",sep="_")
  
  #3 POPULATE EFFLUX/INFLUX:
  overlap_efflux_genes<-intersect(colnames(transport),rownames(matrix))
  tmp_efflux_model<-transport[,overlap_efflux_genes] %*% matrix[overlap_efflux_genes,]
  
  rownames(tmp_efflux_model)<-paste(rownames(tmp_efflux_model),"efflux",sep="_")
  
  #4 COMBINE:
  MoA_matrix<-rbind(tmp_efflux_model,tmp_metab_model)
  
}
#MAKE MoA matrix: TOTAL
make_MoA_model_TOTAL<-function(matrix,transport,metab){
  #1 convert log2 to tpm:
  tmp_matrix<-(2^(matrix)-1)
  
  
  #2 MAKE MASTER MATRIX:
  all_genes<-unique(c(colnames(transport),colnames(metab)))
  all_drugs<-unique(c(rownames(transport),rownames(metab)))
  
  master_matrix<-matrix(0,ncol=length(all_genes),nrow=length(all_drugs))
  rownames(master_matrix)<-all_drugs
  colnames(master_matrix)<-all_genes
  
  
  
  
  for (gene in all_genes){
    for (drug in all_drugs){
      #CHECK METABOLISM:
      if (drug %in% rownames(metab)){
        if (gene %in% colnames(metab)){
          master_matrix[drug, gene]<-metab[drug, gene]
        }
      }
      #CHECK TRANSPORT
      if (drug %in% rownames(transport)){
        if (gene %in% colnames(transport)){
          master_matrix[drug, gene]<-transport[drug, gene]
        }
      }
    }
  }
  
  
  #2 POPULATE METAB:
  overlap_genes<-intersect(colnames(master_matrix),rownames(matrix))
  tmp_total_model<-master_matrix[,overlap_genes] %*% matrix[overlap_genes,]
  
  
  
}


#VISUALIZE
heatmap_Cluster<-function(matrix,title,row_size=0.8,col_size=1){
  color_gradient<-colorRampPalette(c("blue", "white","red"))(n = 20)
  
  #SORT GENES:
  row_means<-rowMeans(matrix,na.rm=T)#take sum of rows
  row_means_sort<-sort(row_means,decreasing=T)#sort
  
  row_names_sorted<-names(row_means_sort)#extract names
  
  #SORT CANCERS:
  col_means<-colMeans(matrix,na.rm=T)#take sum of columns
  col_means_sort<-sort(col_means,decreasing=T)#sort
  
  col_names_sorted<-names(col_means_sort)
  
  ################
  #SORT PLOT:
  ################
  
  
  heatmap.2(matrix[row_names_sorted,col_names_sorted],col=color_gradient,cexRow = row_size,cexCol = col_size,keysize = 1,
            symm=F,density.info = "none",trace="none",margin=c(9,10),main=title,
            dendrogram = "both",
            key.xlab = "zscore(SUM)")
}
hcluster_corr<-function(matrix,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",key="z-score"){
  #distance matrix:
  cor_columns <- cor(matrix,method="spearman")
  cor_rows <- cor(t(matrix),method="spearman")
  
  
  #distance matrix:
  dist_columns <- as.dist(1-cor_columns)
  dist_rows <- as.dist(1-cor_rows)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), 
            cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white","red"))(n = 20),
            density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = key,key.xlab = key)
}
hcluster<-function(matrix,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",key="z-score"){
  #distance matrix:
  dist_columns <- dist(t(matrix))
  dist_rows <- dist(matrix)
  
  #hclustering:
  hclust_columns<-hclust(dist_columns, method="average")
  hclust_rows<-hclust(dist_rows, method="average")
  
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=as.dendrogram(hclust_rows),Colv=as.dendrogram(hclust_columns), 
            cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("blue", "white","red"))(n = 20),
            density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = key,key.xlab = key)
}


################################################################################################################################################################
#FIGURE 5B:  estimate metabolism and transport at Barrier Tissues (using single-cell RNAseq data from Human Cell Atlas)
################################################################################################################################################################

################################################
#LOAD DATA
################################################
load("scRNAseq_db/HCA_matrix_tpm.RData")
dim(HCA_matrix_tpm)
hist(log10(HCA_matrix_tpm+1),breaks=100)

################################################
#SELECT & AVERAGE TPM Expression of Barrier Tissues
################################################

barrier_cells<-c(colnames(HCA_matrix_tpm)[grep("endothelial",colnames(HCA_matrix_tpm))],
                 colnames(HCA_matrix_tpm)[grep("enterocytes",colnames(HCA_matrix_tpm))],
                 colnames(HCA_matrix_tpm)[grep("tubular",colnames(HCA_matrix_tpm))],
                 colnames(HCA_matrix_tpm)[grep("hepatocytes",colnames(HCA_matrix_tpm))])

#MAKE cell-type avg:
lineage_ordered_ccle<-sapply(strsplit(colnames(HCA_matrix_tpm[,barrier_cells]),"_c-"), function(x) x[1])




lineage_filtered<-unique(lineage_ordered_ccle)

#MAKE & POPULATE MATRIX:
HCA_matrix_tpm_SUM<-matrix(ncol=0,nrow=nrow(HCA_matrix_tpm))
rownames(HCA_matrix_tpm_SUM)<-rownames(HCA_matrix_tpm)

#lineage=lineage_filtered[1]
for (lineage in lineage_filtered){
  tmp_idx<-which(lineage_ordered_ccle==lineage)
  if (length(tmp_idx)>1){
    tmp_vector<-rowMeans(HCA_matrix_tpm[,tmp_idx],na.rm=T)
    HCA_matrix_tpm_SUM<-cbind(HCA_matrix_tpm_SUM,tmp_vector[rownames(HCA_matrix_tpm_SUM)])
  }else{
    tmp_vector<-HCA_matrix_tpm[,tmp_idx]
    HCA_matrix_tpm_SUM<-cbind(HCA_matrix_tpm_SUM,tmp_vector[rownames(HCA_matrix_tpm_SUM)])
  }
  
}
colnames(HCA_matrix_tpm_SUM)<-lineage_filtered

#Rename to simplify labels for figure (too long otherwise)
colnames(HCA_matrix_tpm_SUM)<-gsub("proximal enterocytes_","oral-",
                                   gsub("distal enterocytes_","oral-",
                                        gsub("tubular cells_kidney","renal",
                                             gsub("_liver","",
                                                  gsub("endothelial cells_",
                                                       "blood-",
                                                       gsub("peritubular cells_","",gsub("eye","brain",gsub("skeletal ","",colnames(HCA_matrix_tpm_SUM)))))))))


################################################
#CALC efflux & metab by summing TPM
################################################

#RAW SUM:  actual rate estimate
Barrier_efflux_metab<-make_MoA_model(HCA_matrix_tpm_SUM,
                                     chemo_transport_matrix,
                                     chemo_metab_matrix)
#TAKE Z-score so units similar
Barrier_efflux_metab_zscore<-t(scale(t(Barrier_efflux_metab[which(rowSums(Barrier_efflux_metab)>0),]),center=T,scale=T))

#FIGURE 5B:
hcluster(Barrier_efflux_metab_zscore,title="Transport of Chemo at Barrier Tissues",row_size=0.7,col_size = 1)









################################################################################################################################################################
#FIGURE 5D:
################################################################################################################################################################

#######################
#LOAD DATA & Simplify Labels
#######################

load("scRNAseq_db/monaco_matrix_tpm.RData")

dim(monaco_matrix_tpm)#30 immune types
hist(log10(monaco_matrix_tpm+1),breaks=100)#confirm tpm units (log normal distributed)



colnames(monaco_matrix_tpm)<-gsub("CD4 T-cell","CD4",gsub("CD8 T-cell","CD8",
                                                          gsub("Exhausted","Exh",
                                                               gsub("Effector","Eff",
                                                                    gsub("memory","mem",
                                                                         gsub("Memory","Mem",
                                                                              gsub("CD4 T-cell T","T",
                                                                                   gsub("mediate","",
                                                                                        gsub("classical","class",
                                                                                             gsub("effector","eff",
                                                                                                  gsub("Terminal","Term",gsub("witched","wit",colnames(monaco_matrix_tpm)))))))))))))

########################
#CALCULATE METAB/EFFLUX
########################

#RAW SUM:  actual rate estimate
Immune_efflux_metab<-make_MoA_model(monaco_matrix_tpm,
                                    chemo_transport_matrix,
                                    chemo_metab_matrix)

Immune_efflux_metab_zscore<-t(scale(t(Immune_efflux_metab[which(rowSums(Immune_efflux_metab)>0),]),center=T,scale=T))
hcluster(Immune_efflux_metab_zscore,title="Immune efflux/metab Model",row_size=0.7,col_size = 1)
hcluster_corr(Immune_efflux_metab_zscore,title="Immune cell Chemotherapy Transport",row_size=0.7,col_size = 1)


























