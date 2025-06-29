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
#       FIGURE #2:      Comparison of clinical (NCI1970meta), in vitro chemo-sensitivity data with FDA indications
#       FIGURE #3:      Evaluation of AUC/IC50 correlation with ORR, Cmax & TRmin
#       FIGURE #4:      Biomarker definitions and correlation with clinical (an invitro) chemo-sensitivity
#       FIGURE #S7-9    Multi-linear regression of expr-biomarkers demonstrated independent correlations with drug-resposne
#       FIGURE #S10     Literature experimental data on anthracycline intracellular concentration and IC50 changes (to validate chemo-resistance model)
#       FIGURE #5:      Mouse BBB exact prediction of KO-effects
#
#################################################################################################################################################################################################################################################################### 





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
##     ##          ##   ##      ##              #################            ##                                              
##     ##          ##    ##     ##                  ##   ##                  ##                                      
##     ##          ##     #######                   ##   ##           ########                                                     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##       #########     ########                                                            ########                                         
##       ##           ##                                                                  ##                            
##       ##          ##          ######   #####               ##      ##   #####         ##                                                 
##       #########   ##          ##      ##   ##               ##    ##   ##             ##          ###  ###    ##    ##    ##                                     
##       ##          ##          ######  ##   ##                ##  ##     #####         ##          ## ## ##  ##  ##    ####                                    
##       ##           ##             ##  ##   ##                 ####          ##         ##         ##    ##  ######    ####                          
##       #########     ########  ######   #####                   ##       #####           ########  ##    ##  ##  ##  ##    ##                                                  
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################

load("lab_v_clinic/AUC_ORR_compare.RData")


linear_fit_plot_pvalue<-function(x,y,text,pt_colors="black",title="",xlab="",ylab=""){
  #linear regression:
  lin_reg<-lm(y~x)
  #make fit equation + r-squared:
  fit_coeff<-round(coef(lin_reg),2)
  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  r2 <- round(sqrt(summary(lin_reg)$r.squared), 2)
  eq_r2<-paste("r = ", r2," "," "," ","p = ",signif(pvalue_corr,2))
  #plot data:
  plot(x,y, pch = 16, cex = 0.7, col = pt_colors, main=title,xlab=xlab,ylab=ylab)
  #add fit-line and fit-equation:
  abline(lin_reg)
  mtext(eq_r2, 3, line=-2,cex=1)
  text(x,y,text,cex=0.3)
}

linear_fit_plot_pvalue_format<-function(x,y,text,pt_colors="black",title="",xlab="",ylab="",
                                        ylim=c(0,0.7),xlim=c(0.15,1.35)){
  #linear regression:
  lin_reg<-lm(y~x)
  #make fit equation + r-squared:
  fit_coeff<-round(coef(lin_reg),2)
  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  r2 <- round(sqrt(summary(lin_reg)$r.squared), 2)
  eq_r2<-paste("r = ", r2," "," "," ","p = ",signif(pvalue_corr,2))
  #plot data:
  plot(x,y, pch = 16, cex = 0.7, col = pt_colors, main=title,xlab=xlab,ylab=ylab,ylim=ylim,xlim=xlim)
  #add fit-line and fit-equation:
  abline(lin_reg)
  mtext(eq_r2, 3, line=-2,cex=1)
  text(x,y,text,cex=0.3)
}


############################################################
# AUC vs CLINICAL RESPONSE RATES:
############################################################

linear_fit_plot_pvalue_format((1-AUC_ORR_compare$PRISM),AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)")


prism_avg<-tapply(AUC_ORR_compare$PRISM, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
prism_avg<-prism_avg[!is.na(prism_avg)]
orr_avg<-tapply(AUC_ORR_compare$ORR, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
orr_avg<-orr_avg[!is.na(orr_avg)]


linear_fit_plot_pvalue_format((1-prism_avg[names(prism_avg)]),orr_avg[names(prism_avg)],names(prism_avg),
                              title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)",ylim=c(0.08,0.35))


############################################################
# Cmax
############################################################


tmp_Cmax<-read.csv("lab_v_clinic/Cmax_database-v04.csv",as.is=T,header=T)

#MAKE NUMERIC:
tmp_Cmax$AUC.ng.hr.mL.<-as.numeric(tmp_Cmax$AUC.ng.hr.mL.)
tmp_Cmax$Cmax.ng.mL.<-as.numeric(tmp_Cmax$Cmax.ng.mL.)
tmp_Cmax$Cmax.μmol.L.<-as.numeric(tmp_Cmax$Cmax.μmol.L.)
tmp_Cmax$Proteinbinding<-(1-as.numeric(tmp_Cmax$Proteinbinding))


#NEW METRICS:
tmp_Cmax$Cmax_unbound<-(tmp_Cmax$Proteinbinding * tmp_Cmax$Cmax.μmol.L.)
tmp_Cmax$AUC_uM<-tmp_Cmax$AUC.ng.hr.mL.* ( tmp_Cmax$Cmax.μmol.L./ tmp_Cmax$Cmax.ng.mL.)
tmp_Cmax$AUC_uM_fu<-(tmp_Cmax$Proteinbinding * tmp_Cmax$AUC_uM)


#############################################
#PERFORM COMPARISONS
#############################################


#COMPARISON VECTOR:
Cmax_vector<-log10(tmp_Cmax$Cmax_unbound)
names(Cmax_vector)<-tolower(tmp_Cmax$Genericname)


###############
#PRISM:
###############
prism_avg<-AVG_responses$PRISM
overlap_drugs<-intersect(names(Cmax_vector),names(prism_avg))

linear_fit_plot_pvalue(prism_avg[overlap_drugs],Cmax_vector[overlap_drugs],overlap_drugs, title="PRISM vs Cmax",ylab="Clinical Conc (Cmax)",xlab="In vitro Response (1-AUC)")

linear_fit_plot_pvalue_format((1-prism_avg[overlap_drugs]),Cmax_vector[overlap_drugs],overlap_drugs,
                              title="PRISM vs Cmax",ylab="Clinical Cmax Log10(uM)",xlab="In vitro Response (AUC)",ylim=c(-3,3))





############################################################
# Cmax
############################################################


load("lab_v_clinic/MASTER_pk_database_MW.RData")

rownames(PKDB_mw)<-tolower(rownames(PKDB_mw))

#############################################
#PERFORM COMPARISONS
#############################################


#COMPARISON VECTOR:
Cmax_vector<-log10(PKDB_mw[,"thera_min"]/(1000*PKDB_mw[,"MW"]))
names(Cmax_vector)<-rownames(PKDB_mw)


names(AVG_responses$PRISM)[-which(names(AVG_responses$PRISM) %in% names(Cmax_vector))]

###############
#PRISM: 18 present
###############
prism_avg<-AVG_responses$PRISM
overlap_drugs<-intersect(names(Cmax_vector),names(prism_avg))

linear_fit_plot_pvalue(prism_avg[overlap_drugs],Cmax_vector[overlap_drugs],overlap_drugs, title="PRISM vs Cmax",ylab="Clinical Conc (Cmax)",xlab="In vitro Response (1-AUC)")

linear_fit_plot_pvalue_format((1-prism_avg[overlap_drugs]),Cmax_vector[overlap_drugs],overlap_drugs,
                              title="PRISM vs TRmin",ylab="Clinical TRmin Log10(M)",xlab="In vitro Response (AUC)",ylim=c(-10,-4))






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
ida_transport<-read.csv("anthracycline_data/raw-idarubicin-KINETIC-data-v05.csv",as.is=T,header=T,row.names = 1)[1:8,]

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
##     ########    ##     #######                    ##   ##            ##########                                                      
##     ##          ##    ##                          ##   ##            ##                                 
##     ##          ##   ##                      #################       ##                                      
##     ########    ##   ##   ####                    ##   ##             ######                                             
##     ##          ##   ##      ##              #################             ##                                              
##     ##          ##    ##     ##                  ##   ##                   ##                                    
##     ##          ##     #######                   ##   ##             #######                                           
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
#EVALUATE JOINT-INDEPENDENT MICHAELIS-FLUXES:
# KO-ratio-G2 + KO-ratio-B1 = KO-ratio-TKO
#
#  EQ:  [OUT]/[IN] = kefflux / kdiff   = (kg1 + kb1) / kdiff = (kg2/kdiff) + (kb1/kdiff)
#
#GOAL:   show that individual BR's are sufficient to explain TKO-BR (= validation of chemo-resist model trying to publish in NCI1970 paper)
#
#RESULT:  CORRELATION NEARLY PERFECT (log scale and ALL linear scales!!!!!!!!!!!!!!!!!!)


load("mouse_db/BBB_matrix.RData")



#NEEDED LIBRARIES:
library("gplots")
library("RColorBrewer")

#LINEAR REGRESSION:
linear_fit_plot<-function(x,y,pt_colors="black",xlab=NA,ylab=NA,title=NA,xlim=NA,ylim=NA){
  #linear regression:
  lin_reg<-lm(y~x)
  #make fit equation + r-squared:
  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  fit_coeff<-round(coef(lin_reg),6)
  r2 <- round(summary(lin_reg)$r.squared, 2)
  rmse <- round(sqrt(mean(resid(lin_reg)^2)), 2)
  eq_r2<-paste("r = ", r2," "," "," ","p-value = ",signif(pvalue_corr,2))
  #plot data:
  plot(x,y, pch = 16, cex = 0.8, col = pt_colors,xlab=xlab,ylab=ylab,main=title,xlim=xlim,ylim=ylim)
  #add fit-line and fit-equation:
  abline(lin_reg,col=rgb(1,0,0,0.5))
  mtext(eq_r2, 3, line=-2,cex=1.5,col=rgb(1,0,0,0.5))
}


####################################################################################################
#EVALUATE 1/g2 + 1/b1 = 1/wt    (previous section evaluated:  tko/g2 + tko/b1 = tko/wt)
####################################################################################################


#ABCB1 function data:
invert_b1<-BBB_matrix[,"abcb1a -/-"]^-1
invert_g2<-BBB_matrix[,"abcg2 -/-"]^-1
invert_wt<-BBB_matrix[,"WT"]^-1

invert_sum<-invert_b1+invert_g2


linear_fit_plot(log10(invert_sum),log10(invert_wt),xlim=c(-2,3),ylim=c(-2,3),
                xlab="Abcb1 (1/g2) + Abcg2 (1/b1)",ylab="TKO (1/wt)",title="Parallel Independent Model (1/g2 + 1/b1) ")
text(log10(invert_sum),log10(invert_wt),names(invert_sum),cex=0.5)
abline(a=0,b=1,lty=2)



#ALL B1 model:
linear_fit_plot(log10(invert_g2),log10(invert_wt),xlim=c(-2,3),ylim=c(-2,3),
                xlab="Abcb1 (1/g2)",ylab="TKO (1/wt)",title="All B1 Model (1/g2 = 1/wt)")
text(log10(invert_g2),log10(invert_wt),names(invert_sum),cex=0.5)
abline(a=0,b=1,lty=2)


#ALL G2 model:
linear_fit_plot(log10(invert_b1),log10(invert_wt),xlim=c(-2,3),ylim=c(-2,3),
                xlab="Abcg2 (1/b1)",ylab="TKO (1/wt)",title="All G2 Model (1/b1 = 1/wt)")
text(log10(invert_b1),log10(invert_wt),names(invert_sum),cex=0.5)
abline(a=0,b=1,lty=2)



linear_fit_plot(log10(invert_sum),log10(invert_wt),xlim=c(-0.5,3),ylim=c(-0.5,3),
                xlab="Abcb1 (1/g2) + Abcg2 (1/b1)",ylab="TKO (1/wt)",title="Parallel Independent Model (no TKO data) ")
text(log10(invert_sum),log10(invert_wt),names(invert_sum),cex=0.5)


















