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
#DATE:   11-18-2025
#GOAL:   Update NCI1970 paper for publication (based on Jon's feedback)
#STRUCTURE:
#       FIGURE #1:      Data Gap Overview
#       FIGURE #2:      Single-Agent History data
#       FIGURE #3:      Rosetta Stone:   Clinic vs FDA vs In vitro
#       FIGURE #4:      Evaluation of AUC/IC50 correlation with ORR, Cmax & TRmin
#       FIGURE #5:      Global Biomarker Evaluation
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


##############################################################################################################################################
#LOAD DATA AND STATS:
##############################################################################################################################################


library(tidyverse)

# Read the CSV
df <- read.csv("clinical_trials.csv", stringsAsFactors = FALSE)

# Ensure Year is numeric
df$Year <- as.numeric(df$Year)
pre1970_idx<-which(df$Seconday.Source=="Single Agents in Cancer Chemotherapy")

#######################
#NCI STATS:
#######################
length(unique(paste(df$Journal[pre1970_idx],df$Year[pre1970_idx],df$Volume[pre1970_idx],df$Page[pre1970_idx])))
sum(df$total.patients[pre1970_idx])

#######################
#ALL STATS:
#######################

#TOTAL NUMBER PAPERS:
length(unique(paste(df$Journal,df$Year,df$Volume,df$Page)))

sum(df$total.patients)


##############################################################################################################################################
#FIG #2 pieces:   compilation of historical data
##############################################################################################################################################





# ---- 1. Histogram: Number of publications per year ----
ggplot(df, aes(x = Year)) +
  geom_histogram(binwidth = 1, fill = "lightgrey", color = "white") +
  labs(title = "Number of Publications per Year",
       x = "Year",
       y = "Number of Publications") +
  theme_minimal(base_size = 14)

# ---- 2. Boxplot: Range of publication years for top 20 drugs (sorted by mean) ----

# Identify top 20 drugs by number of publications
top_drugs <- df %>%
  count(drug, sort = TRUE) %>%
  slice_head(n = 30) %>%
  pull(drug)

# Filter to those drugs
df_top <- df %>% filter(drug %in% top_drugs)

# Compute mean publication year for sorting
drug_means <- df_top %>%
  group_by(drug) %>%
  summarize(mean_year = mean(Year, na.rm = TRUE))

# Join to main dataset
df_top <- df_top %>%
  left_join(drug_means, by = "drug")

# Create boxplot sorted by mean year
ggplot(df_top, aes(x = reorder(drug, -mean_year), y = Year)) +
  geom_boxplot(fill = "lightgrey", outlier.shape = 16, outlier.size = 2) +
  coord_flip() +
  labs(title = "Range of Publication Years for Top Chemos",
       x = "Drug",
       y = "Publication Year") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 8),   # <-- reduce drug name font size
    axis.text.x = element_text(size = 10)   # optional: adjust x-axis size too
  )




##############################################################################################################################################
#FIG 2: complete figure
##############################################################################################################################################

library(tidyverse)
library(patchwork)

# ---- Read and clean ----
df <- read.csv("clinical_trials.csv", stringsAsFactors = FALSE)

# Coerce Year to numeric, suppressing non-numeric warnings
df <- df %>%
  mutate(
    Year = as.numeric(str_extract(Year, "\\d{4}")),  # extract 4-digit year if mixed text
  ) %>%
  filter(!is.na(Year) & Year > 1900 & Year < 2100)   # keep only plausible years

# ---- Histogram ----
p_hist <- ggplot(df, aes(x = Year)) +
  geom_histogram(
    aes(y = ..density..),       # scale histogram to density for overlay
    binwidth = 2,
    fill = "grey",
    color = "grey",
    boundary = 0.5,
    alpha = 1                 # make bars slightly transparent
  ) +
  geom_density(
    color = "black",
    size = 0.5,
    adjust = 1.2,               # smoothness; smaller = bumpier, larger = smoother
    na.rm = TRUE
  ) +
  labs(y = "Density (Publications per Year)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_blank(),
    plot.margin = margin(0, 10, 0, 10)
  )

# ---- Boxplot (top 20 drugs) ----
top_drugs <- df %>%
  count(drug, sort = TRUE) %>%
  slice_head(n = 30) %>%
  pull(drug)

df_top <- df %>% filter(drug %in% top_drugs)

# Compute mean years
drug_means <- df_top %>%
  group_by(drug) %>%
  summarize(mean_year = mean(Year, na.rm = TRUE), .groups = "drop")

# Merge and reorder factor levels by mean year (recent on right)
df_top <- df_top %>%
  left_join(drug_means, by = "drug") %>%
  mutate(drug = fct_reorder(drug, mean_year))

# ---- Boxplot (vertical, shared x-axis) ----
p_box <- ggplot(df_top, aes(x = Year, y = drug)) +
  geom_boxplot(fill = "grey", outlier.shape = 16, outlier.size = 1.5) +
  labs(x = "Publication Year", y = "Drug") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.margin = margin(0, 10, 0, 10)
  )

# ---- Combine ----
combined_plot <- p_hist / p_box +
  plot_layout(heights = c(0.25, 0.75), guides = "collect") &
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

combined_plot <- combined_plot +
  plot_annotation(title = "History of Single-Agent Chemotherapy Publications")

# Display
combined_plot





##############################################################################################################################################
#SUPP FIGURES:   numbers of papers per drug (does it correlate with time?)
##############################################################################################################################################

library(ggplot2)

num_papers <- df %>%
  count(Drug, sort = TRUE) 

#num papers vector:
num_papers_vector<-num_papers$n
names(num_papers_vector)<-num_papers$Drug


barplot(num_papers_vector[1:20],horiz=T,las=2,cex.names = 0.5)

#drug year vector:
year_drug_vector<-drug_means$mean_year
names(year_drug_vector)<-drug_means$Drug

plot(year_drug_vector,num_papers_vector[names(year_drug_vector)], pch=16, ylim =c(0,130))
text(year_drug_vector,num_papers_vector[names(year_drug_vector)],names(year_drug_vector),cex=0.5)


linear_fit_plot<-function(x,y,pt_colors="black"){
  #linear regression:
  lin_reg<-lm(y~x)
  #make fit equation + r-squared:
  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  fit_coeff<-round(coef(lin_reg),2)
  r2 <- round(summary(lin_reg)$r.squared, 2)
  rmse <- round(sqrt(mean(resid(lin_reg)^2)), 2)
  eq_r2<-paste("y = ", fit_coeff[2], "x + ", fit_coeff[1]," "," "," ","r^2 = ", r2," "," "," ","rmse = ",rmse)
  #plot data:
  plot(x,y, pch = 16, cex = 1.3, col = pt_colors)
  #add fit-line and fit-equation:
  abline(lin_reg)
  mtext(eq_r2, 3, line=-2,cex=1)
}
linear_fit_plot_pvalue<-function(x,y,pt_colors="black"){
  #linear regression:
  lin_reg<-lm(y~x)
  #make fit equation + r-squared:
  fit_coeff<-round(coef(lin_reg),2)
  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  r2 <- round(sqrt(summary(lin_reg)$r.squared), 2)
  eq_r2<-paste("r = ", r2," "," "," ","p-value = ",round(pvalue_corr,3))
  #plot data:
  plot(x,y, pch = 16, cex = 1.3, col = pt_colors,
       xlab="Chemistry-Grouping P-value Rank",ylab="Average Correlation Drug-Targets",cex.lab=1.5,cex.axis=1.5,
       main="Heterogenity in VIPER-effectors = F(Heterogentity Targets)",cex.main=2)
  #add fit-line and fit-equation:
  abline(lin_reg)
  mtext(eq_r2, 3, line=-2,cex=2)
}




linear_fit_plot(year_drug_vector[c(-5,-13)],num_papers_vector[names(year_drug_vector)][c(-5,-13)])
linear_fit_plot_pvalue(year_drug_vector[c(-5,-13)],num_papers_vector[names(year_drug_vector)][c(-5,-13)])


#GGPLOT PLOTS:
year_pubs_df<-data.frame("names"=names(year_drug_vector),
                         "year"=year_drug_vector,
                         "pubs"=num_papers_vector[names(year_drug_vector)])

ggplot(year_pubs_df[-which(year_pubs_df$names %in% c("cyclophosphamide","fluorouracil")),],aes_string(x="year",y="pubs"))+
  geom_point(size=1)+
  geom_smooth(method="lm")












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
#SETUP:  packages, functions, data and row/column order
####################################################################################################################################################################################################################################################################  

#------------------------------------------------------------------------------
#LOAD LIBRARIES AND FUNCTIONS:
#------------------------------------------------------------------------------


#NEEDED LIBRARIES:
library("gplots")
library("RColorBrewer")
library(patchwork)
library(tidyverse)
library(ggplot2)
library(dplyr)

#PLOTTING FUNCTION:
hcluster_unsortBW<-function(matrix,type,col_size=1,row_size=1,x_label=NA,ylabel=NA,title="Hierarchical Clustering",units="zscore"){
  
  #PLOT HEATMAP:
  heatmap.2(matrix,margin=c(7,10), Rowv=NA,Colv=NA, cexCol =col_size,cexRow = row_size,col=colorRampPalette(c("white", "black"))(n = 20),
            density.info = "none",trace="none",xlab=x_label,main=title,keysize=1,key.title = units,key.xlab = units,symm=F,dendrogram="none")
  
  return(rowMeans(matrix))
}
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
boot_r2 <- function(formula, data, R = 1000) {
  r2_vals <- numeric(R)
  n <- nrow(data)
  
  for (i in 1:R) {
    idx <- sample(seq_len(n), replace = TRUE)
    d <- data[idx, ]
    r2_vals[i] <- summary(lm(formula, data = d))$r.squared
  }
  
  tibble(
    mean = mean(r2_vals),
    se   = sd(r2_vals),
    lwr  = quantile(r2_vals, 0.025),
    upr  = quantile(r2_vals, 0.975)
  )
}


#LOAD DATA Construction
load("early_db/NCI1970-monotherapy-DB_simple.RData")
load("early_db/NCI1970-META_monotherapyDB.RData")
load("NCI1970meta_Database.RData")
#LOAD Clincial Context:
load("clinical_dbs/SEER_cancer_incidence.RData")
load("clinical_dbs/FDA_indictations.RData")

#LOAD IN VITRO DATA:
load("cell-line_dbs/invitro_AUC-lineage-AVG_dbs.RData")



####################################################################################################################################################################################################################################################################  
#FIG S3:   Data Construction
####################################################################################################################################################################################################################################################################  

#------------------------------------------------------------------------------
#BUILD MATRICES WITH GAPS TO SHOW DATA ADDITION
#------------------------------------------------------------------------------


#sort drugs & cancers for BW-heatmap
sorted_cancers<-names(sort(colSums(NCI1970_META_numbers),decreasing=T))
sorted_drugs<-names(sort(rowSums(NCI1970_META_numbers),decreasing=T))
drugs_na_num<-sort(rowSums(is.na(data.frame(NCI1970_META_ORrate))))
good_drugs<-names(drugs_na_num)[1:34]

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


#------------------------------------------------------------------------------
#STATISTICS ON MATRICES FOR PAPER DISCUSSION
#------------------------------------------------------------------------------


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

#Patients per cancer
cancer_num<-sort(colSums(NCI1970_update_number_FINAL),decreasing=T)
#brca nsclc  coad lymph    ov  leuk  skcm  hnsc  sarc  cesc  paad  tgct  kirc  stad  blca  esca  ucec   gbm  lihc 
#7070  6136  5555  5230  4982  4668  3004  2892  1422  1093  1065  1061  1056  1030   856   635   425   418   404 


#------------------------------------------------------------------------------
#FINAL PLOTS
#------------------------------------------------------------------------------


#Figure S3A
par(las=1)
hcluster_unsortBW(t(log10(NCI1970_orig_number_FINAL+1)),units="log10(patient #)",title="NCI1970 Textbook alone",
                  col_size = 0.7,row_size = 1.2)



#Figure S3B
par(las=1)
hcluster_unsortBW(t(log10(NCI1970_update_number_FINAL+1)),units="log10(patient #)",title="NCI1970 + Meta Reviews",
                  col_size = 0.7,row_size = 1.2)


#Figure S3C
par(las=2,mar=c(5, 4, 4, 2)) # make label text perpendicular to axis
barplot(sort(colSums(NCI1970_update_number_FINAL),decreasing=F), main="NCI1970:  # Patients / Cancer", 
        beside=F,cex.axis = 1.5,cex=1.4,cex.main=2,cex.names=1,ylab=c(0,100),col="black",horiz=T)

#------------------------------------------------------------------------------
#SEER EXPLAINS PATIENT FEQUENCY
#------------------------------------------------------------------------------

seer_vector<-SEER_cancer_incidence[names(cancer_num)]

#BEST FIT LINE:  y = 96.44x + 861.98       r^2 = 0.58    rmse = 1436.51
linear_fit_plot(seer_vector,cancer_num)

# Figure S3D
plot(seer_vector,cancer_num,pch=16,xlim=c(0,80),ylim=c(0,8000),xlab="Cancer Incidence (Millions)",ylab="# patients")
text(seer_vector+2.5,cancer_num,labels=names(cancer_num))
abline(a=0,b=121.18)





####################################################################################################################################################################################################################################################################  
#FIG 1A-B:   NCI1970 Data Understanding
####################################################################################################################################################################################################################################################################  

#------------------------------------------------------------------------------
#Organization for NCI1970 matrix (filter & sort rows and columns)
#------------------------------------------------------------------------------


#Filter Matrix by Top 18 cancers w/ least NA's
cancers_na_num<-sort(colSums(is.na(data.frame(NCI1970meta_rate))))
good_cancers<-names(cancers_na_num)[1:18]
filtered_ORR<-NCI1970meta_rate[,good_cancers]

#Sort Rows and Columns by avg-ORR
orr_sorted_cancers<-names(sort(colMeans(filtered_ORR,na.rm=T),decreasing=T))
orr_sorted_drugs<-names(sort(rowMeans(filtered_ORR,na.rm=T),decreasing=T))

#Filter by Drugs in FDA-indications matrix
drug_overlap<-orr_sorted_drugs[which(orr_sorted_drugs %in% rownames(FDA_indictations))]

#------------------------------------------------------------------------------
#Fig 3A:   ORR heatmap
#------------------------------------------------------------------------------

#Figure 2A
hcluster_unsort(t(NCI1970meta_rate[drug_overlap,orr_sorted_cancers]),title="NCI1970meta Response Rates",units="% Response",col_size = 0.7,row_size = 0.7)


#------------------------------------------------------------------------------
#Fig 3B:  Boxplot + simple-ANOVA (transpose = drugs)
#------------------------------------------------------------------------------
#CANCER ANOVA:  F = 12.33, p-value = <2e-16
#DRUG ANOVA:  F = 1.942, p-value = 0.00331


# Assume your matrix is NCI1970meta_rate

# 1. Convert matrix → long form
df_long <- as.data.frame(NCI1970meta_rate[drug_overlap,orr_sorted_cancers]) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  drop_na()  # remove NAs

# 2. Compute column means (ignoring NAs)
col_means <- df_long %>%
  group_by(variable) %>%
  summarize(mean_val = mean(value, na.rm = TRUE))

# 3. Reorder factor levels by mean
df_long <- df_long %>%
  mutate(variable = factor(variable, levels = col_means$variable[order(col_means$mean_val)]))

# 4. Boxplot sorted by mean
ggplot(df_long, aes(x = variable, y = value)) +
  geom_boxplot(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "",
    y = "Value",
    title = "Boxplots of NCI1970meta_rate Columns\n(sorted by column mean)"
  ) +
  theme_bw()

anova_result <- aov(value ~ variable, data = df_long)
summary(anova_result)

#------------------------------------------------------------------------------
#Fig 3B-TEST:  two way anova
#------------------------------------------------------------------------------
#F-statistics:  cancer (13.43) >> drug (3.29)
#Percent Variance:  Cancer (37.3%), >> drug (15.6%)
#p-value:  cancer (<2.2e-16) << drug (1.43e-7)

# suppose your matrix is called M
# rows = e.g., genes, columns = samples

# turn matrix into long format
df <- as.data.frame(as.table(NCI1970meta_rate[drug_overlap,orr_sorted_cancers]))
colnames(df) <- c("drug", "cancer", "value")

# drop NAs if present
df <- na.omit(df)

# treat row and col as factors
df$rdrug <- factor(df$drug)
df$cancer <- factor(df$cancer)

# two-way ANOVA without interaction (we mostly care about SS decomposition)
fit <- aov(value ~ drug + cancer, data = df)
anova_tab <- anova(fit)
anova_tab

#Percent of the variance
ss <- anova_tab$"Sum Sq"
names(ss) <- rownames(anova_tab)
prop_var <- ss / sum(ss)
prop_var


######################################################################################################################################################################
#FIGURE 1C-D:  Heatmap of FDA-indications (A) correlation of avg-ORR with total # FDA-indications (B) 
######################################################################################################################################################################



# Figure 2C

#------------------------------------------------------------------------------
#Fig 3C:  FDA-indications heatmap
#------------------------------------------------------------------------------

hcluster_unsortFDA(FDA_indictations[drug_overlap,orr_sorted_cancers],title="FDA Indications",units="binary",col_size = 1.3,row_size = 0.7)

#------------------------------------------------------------------------------
#Fig 3C-SUPP:  2-way anova (analogous to NCI1970)
#------------------------------------------------------------------------------
#F-statistics:  cancer (5.55) >> drug (2.19)
#Percent Variance:  Cancer (14.5%), >> drug (9.7%)
#p-value:  cancer (1.67e-11) << drug (4.32e-4)

# turn matrix into long format
df <- as.data.frame(as.table(FDA_indictations[drug_overlap,orr_sorted_cancers]))
colnames(df) <- c("drug", "cancer", "value")

# drop NAs if present
df <- na.omit(df)

# treat row and col as factors
df$drug <- factor(df$drug)
df$cancer <- factor(df$cancer)

# two-way ANOVA without interaction (we mostly care about SS decomposition)
fit <- aov(value ~ drug + cancer, data = df)
anova_tab <- anova(fit)
anova_tab

#Percent of the variance
ss <- anova_tab$"Sum Sq"
names(ss) <- rownames(anova_tab)
prop_var <- ss / sum(ss)
prop_var



#------------------------------------------------------------------------------
# FIG 3D:  R² bar plot with SE bars for univariate and bivariate models
#------------------------------------------------------------------------------

#MAKE DATAFRAME
df<-data.frame("mean_ORR"=colMeans(filtered_ORR[,orr_sorted_cancers],na.rm=T),
                           "sum_FDAindic"=colSums(FDA_indictations[,orr_sorted_cancers],na.rm=T),
                           "sum_NCIpatients"=colSums(NCI1970meta_total[,orr_sorted_cancers]),
                           "SEER"=seer_vector[orr_sorted_cancers])
                           
#MAKE MODELS
m1 <- lm(sum_FDAindic ~ mean_ORR, data = df)
m2 <- lm(sum_FDAindic ~ SEER, data = df)
m3 <- lm(sum_FDAindic ~ mean_ORR + SEER, data = df)


#COMPUTE R2
r2_df <- bind_rows(
  boot_r2(sum_FDAindic ~ mean_ORR, df) %>% mutate(model = "avg-ORR"),
  boot_r2(sum_FDAindic ~ SEER, df) %>% mutate(model = "Cancer Freq"),
  boot_r2(sum_FDAindic ~ mean_ORR + SEER, df) %>% mutate(model = "Both")
)


#PLOT
ggplot(r2_df, aes(y = model, x = mean)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(aes(xmin = mean - se, xmax = mean + se),
                width = 0.15, size = 0.8) +
  geom_text(aes(label = sprintf("%.2f", mean)), 
            hjust = -0.3, size = 4) +
  xlim(0, max(r2_df$mean + r2_df$se) * 1.25) +
  labs(title = "Panel A: R² Comparison (Bootstrap SE)",
       x = "R²", y = "") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12))

#JON NOTE:  CHECK WEIGHTED LINEAR REGRESSION ACROSS RAW ORR = f(FDA, SEER)
#.        -> I did this at the end of Figure 4's code cause I had the AUC_ORR_compare df properly setup
#.  RESULT:
#Call:
#  lm(formula = ORR ~ SEER + FDA, data = AUC_ORR_compare, weights = patient_num)
#
#Weighted Residuals:
#  Min      1Q  Median      3Q     Max 
#-5.9745 -0.9266 -0.4027  0.2291  6.9757 
#
#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  0.2200589  0.0168469  13.062  < 2e-16 ***
#  SEER        -0.0009432  0.0003475  -2.714  0.00699 ** 
#  FDA          0.1616000  0.0161286  10.019  < 2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 1.638 on 332 degrees of freedom
#(205 observations deleted due to missingness)
#Multiple R-squared:  0.2546,	Adjusted R-squared:  0.2501 
#F-statistic:  56.7 on 2 and 332 DF,  p-value: < 2.2e-16



######################################################################################################################################################################
#FIGURE 3E-F & S4 & S5:  Heatmap of CTRPv2 
######################################################################################################################################################################

#------------------------------------------------------------------------------
# FIG S4 & FIG 3E:  CTRP:  heatmap + anova
#------------------------------------------------------------------------------
#Percent Var: drug (85.2%) >> Cancer (8.6%)
#F-statistic:  drug (212.2) >> cancer (21.60)
#P-value:  both <2.2e-16

#reconcile formatting drug/cancer names
rownames(ctrp_avg_matrix)<-tolower(rownames(ctrp_avg_matrix))
ctrp_cancer_overlap<-orr_sorted_cancers[which(orr_sorted_cancers %in% colnames(ctrp_avg_matrix))]
ctrp_drug_overlap<-tolower(orr_sorted_drugs[which(tolower(orr_sorted_drugs) %in% rownames(ctrp_avg_matrix))])
NCI1970_META_ORrate_lowercase<-NCI1970meta_rate
rownames(NCI1970_META_ORrate_lowercase)<-tolower(rownames(NCI1970meta_rate))

#Filter matrices
norm_CTRP<-(16-ctrp_avg_matrix[ctrp_drug_overlap,ctrp_cancer_overlap])/16
norm_NCI_ctrp<-NCI1970_META_ORrate_lowercase[ctrp_drug_overlap,ctrp_cancer_overlap]

#PLOT:
hcluster_unsort(norm_CTRP,title="CTRP:  lineage-average",units="(16-AUC)/16   / drug-rt-mean-sq",col_size = 1.3,row_size = 0.9)

#ANOVA:
df <- as.data.frame(as.table(norm_CTRP))
colnames(df) <- c("drug", "cancer", "value")
df <- na.omit(df)
df$drug <- factor(df$drug)
df$cancer <- factor(df$cancer)
fit <- aov(value ~ drug + cancer, data = df)
anova_tab <- anova(fit)
anova_tab
ss <- anova_tab$"Sum Sq"
names(ss) <- rownames(anova_tab)
prop_var <- ss / sum(ss)
prop_var

#------------------------------------------------------------------------------
# FIG S4 and 3F:  Correlation cancer rates
#------------------------------------------------------------------------------


#EXTRACT AVGS:
mean_drug_sens_CTRP<-colMeans(norm_CTRP,na.rm=T)
ctrp_se <- apply(norm_CTRP, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

mean_drug_sens_NCI<-colMeans(norm_NCI_ctrp,na.rm=T)
nci_se <- apply(norm_NCI_ctrp, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))


#BEST FIT LINE:  y = 1.46x -0.2       r^2 = 0.68    rmse = 0.06
linear_fit_plot(mean_drug_sens_CTRP,mean_drug_sens_NCI)


#ACTUAL PLOT:
plot(mean_drug_sens_CTRP,mean_drug_sens_NCI,pch=16,xlim=c(0.2,0.45),ylim=c(0.05,0.5),xlab="AVG in vitro sensitivity (1-AUC)",ylab="Clinical Response Rate")
text(mean_drug_sens_CTRP+0.008,mean_drug_sens_NCI,labels=names(mean_drug_sens_NCI))
abline(a=-0.2,b=1.46)
arrows(mean_drug_sens_CTRP, mean_drug_sens_NCI-nci_se, 
       mean_drug_sens_CTRP, mean_drug_sens_NCI+nci_se, length=0.05, angle=90, code=3,lwd=0.25)
arrows(mean_drug_sens_CTRP-ctrp_se, mean_drug_sens_NCI, 
       mean_drug_sens_CTRP+ctrp_se, mean_drug_sens_NCI, length=0.05, angle=90, code=3,lwd=0.25)



#------------------------------------------------------------------------------
# FIG S4:  Correlation drug rates
#------------------------------------------------------------------------------


#EXTRACT AVGS:
mean_drug_sens_CTRP<-rowMeans(norm_CTRP,na.rm=T)
ctrp_se <- apply(norm_CTRP, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

mean_drug_sens_NCI<-rowMeans(norm_NCI_ctrp,na.rm=T)
nci_se <- apply(norm_NCI_ctrp, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))


#BEST FIT LINE:    r^2 = 0.00    rmse = 0.07
linear_fit_plot(mean_drug_sens_CTRP,mean_drug_sens_NCI)


#ACTUAL PLOT:
plot(mean_drug_sens_CTRP,mean_drug_sens_NCI,pch=16,xlab="AVG in vitro sensitivity (1-AUC)",ylab="Clinical Response Rate")
text(mean_drug_sens_CTRP+0.008,mean_drug_sens_NCI,labels=names(mean_drug_sens_NCI))
abline(a=0.2,b=0.02)
arrows(mean_drug_sens_CTRP, mean_drug_sens_NCI-nci_se, 
       mean_drug_sens_CTRP, mean_drug_sens_NCI+nci_se, length=0.05, angle=90, code=3,lwd=0.25)
arrows(mean_drug_sens_CTRP-ctrp_se, mean_drug_sens_NCI, 
       mean_drug_sens_CTRP+ctrp_se, mean_drug_sens_NCI, length=0.05, angle=90, code=3,lwd=0.25)


######################################################################################################################################################################
#FIGURE S4:  Heatmap & clinical correlations of GDSC
######################################################################################################################################################################

#------------------------------------------------------------------------------
# FIG S4 & FIG 3E:  CTRP:  heatmap + anova
#------------------------------------------------------------------------------
#Percent Var: drug (80.8%) >> Cancer (6.6%)
#F-statistic:  drug (101.9) >> cancer 7.8)
#P-value:  both <2.2e-16

#reconcile formatting drug/cancer names
rownames(gdsc1_avg_matrix)<-tolower(rownames(gdsc1_avg_matrix))

gdsc1_cancer_overlap<-orr_sorted_cancers[which(orr_sorted_cancers %in% colnames(gdsc1_avg_matrix))]
gdsc1_drug_overlap<-tolower(orr_sorted_drugs[which(tolower(orr_sorted_drugs) %in% rownames(gdsc1_avg_matrix))])

#Filter matrices
norm_CTRP<-(1-gdsc1_avg_matrix[gdsc1_drug_overlap,gdsc1_cancer_overlap])
norm_NCI_ctrp<-NCI1970_META_ORrate_lowercase[gdsc1_drug_overlap,gdsc1_cancer_overlap]

#PLOT:
hcluster_unsort(norm_CTRP,title="GDSC:  lineage-average",units="1-AUC",col_size = 1.3,row_size = 0.9)

#ANOVA:
df <- as.data.frame(as.table(norm_CTRP))
colnames(df) <- c("drug", "cancer", "value")
df <- na.omit(df)
df$drug <- factor(df$drug)
df$cancer <- factor(df$cancer)
fit <- aov(value ~ drug + cancer, data = df)
anova_tab <- anova(fit)
anova_tab
ss <- anova_tab$"Sum Sq"
names(ss) <- rownames(anova_tab)
prop_var <- ss / sum(ss)
prop_var

#------------------------------------------------------------------------------
# FIG S4 and 3F:  Correlation cancer rates
#------------------------------------------------------------------------------


#EXTRACT AVGS:
mean_drug_sens_CTRP<-colMeans(norm_CTRP,na.rm=T)
ctrp_se <- apply(norm_CTRP, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

mean_drug_sens_NCI<-colMeans(norm_NCI_ctrp,na.rm=T)
nci_se <- apply(norm_NCI_ctrp, 2, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))


#BEST FIT LINE:  y = 1.46x -0.2       r^2 = 0.68    rmse = 0.06
linear_fit_plot(mean_drug_sens_CTRP,mean_drug_sens_NCI)


#ACTUAL PLOT:
plot(mean_drug_sens_CTRP,mean_drug_sens_NCI,pch=16,xlim=c(0.2,0.45),ylim=c(0.05,0.5),xlab="AVG in vitro sensitivity (1-AUC)",ylab="Clinical Response Rate")
text(mean_drug_sens_CTRP+0.01,mean_drug_sens_NCI,labels=names(mean_drug_sens_NCI),cex=0.9)
abline(a=-0.2,b=1.46)
arrows(mean_drug_sens_CTRP, mean_drug_sens_NCI-nci_se, 
       mean_drug_sens_CTRP, mean_drug_sens_NCI+nci_se, length=0.05, angle=90, code=3,lwd=0.25)
arrows(mean_drug_sens_CTRP-ctrp_se, mean_drug_sens_NCI, 
       mean_drug_sens_CTRP+ctrp_se, mean_drug_sens_NCI, length=0.05, angle=90, code=3,lwd=0.25)



#------------------------------------------------------------------------------
# FIG S4:  Correlation drug rates
#------------------------------------------------------------------------------


#EXTRACT AVGS:
mean_drug_sens_CTRP<-rowMeans(norm_CTRP,na.rm=T)
ctrp_se <- apply(norm_CTRP, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))

mean_drug_sens_NCI<-rowMeans(norm_NCI_ctrp,na.rm=T)
nci_se <- apply(norm_NCI_ctrp, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))


#BEST FIT LINE:    r^2 = 0.00    rmse = 0.07
linear_fit_plot(mean_drug_sens_CTRP,mean_drug_sens_NCI)


#ACTUAL PLOT:
plot(mean_drug_sens_CTRP,mean_drug_sens_NCI,pch=16,xlab="AVG in vitro sensitivity (1-AUC)",ylab="Clinical Response Rate")
text(mean_drug_sens_CTRP+0.03,mean_drug_sens_NCI,labels=names(mean_drug_sens_NCI),cex=0.9)
abline(a=0.22,b=(-0.08))
arrows(mean_drug_sens_CTRP, mean_drug_sens_NCI-nci_se, 
       mean_drug_sens_CTRP, mean_drug_sens_NCI+nci_se, length=0.05, angle=90, code=3,lwd=0.25)
arrows(mean_drug_sens_CTRP-ctrp_se, mean_drug_sens_NCI, 
       mean_drug_sens_CTRP+ctrp_se, mean_drug_sens_NCI, length=0.05, angle=90, code=3,lwd=0.25)









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
##       #########     ########                                                            ########                                         
##       ##           ##                                                                  ##                            
##       ##          ##          ######   #####               ##      ##   #####         ##                                                 
##       #########   ##          ##      ##   ##               ##    ##   ##             ##          ###  ###    ##    ##    ##                                     
##       ##          ##          ######  ##   ##                ##  ##     #####         ##          ## ## ##  ##  ##    ####                                    
##       ##           ##             ##  ##   ##                 ####          ##         ##         ##    ##  ######    ####                          
##       #########     ########  ######   #####                   ##       #####           ########  ##    ##  ##  ##  ##    ##                                                  
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################


library(dplyr)
library(tibble)

load("lab_v_clinic/AUC_ORR_compare_COMPLETE.RData")

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
  print(summary(lin_reg))
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
  print(summary(lin_reg))
}

library(MBESS)

linear_fit_plot_ci<-function(x,y,text,pt_colors="black",title="",xlab="",ylab="",labels=F){
  #linear regression:
  lin_reg<-lm(y~x)

  pvalue_corr <- summary(lin_reg)$coefficients["x","Pr(>|t|)"] 
  
  #########################
  #Confidence Interval
  #########################
  summary_fit <- summary(lin_reg)
  R2 <- summary_fit$r.squared
  N  <- length(y)
  p  <- 1   # one predictor x
  
  ci <- ci.R2(R2 = R2, N = N, K = p, conf.level = 0.95)
  ci
  
  eq_r2<-paste0("r2 = ", round(R2,2),"(",round(ci$Lower.Conf.Limit.R2,2),"-",round(ci$Upper.Conf.Limit.R2,2),")"," p = ",signif(pvalue_corr,2))
  
  
  if (labels==T){
    text(x,y,text,cex=0.3)
  }
  
  #plot data:
  plot(x,y, pch = 16, cex = 0.7, col = pt_colors, main=eq_r2,xlab=xlab,ylab=ylab,cex.main=1)
  #add fit-line and fit-equation:
  abline(lin_reg)
  
}


####################################################################################################################################################################################################################################################################
#0 MAKE ALL DATA:  loaded above
####################################################################################################################################################################################################################################################################

{
#------------------------------------------------------------------------------
# ORIGINAL COMPARISON DATASET:
#------------------------------------------------------------------------------


#LOAD DATA
load("lab_v_clinic/AUC_ORR_compare.RData")

#CORRECT AUC "correction": (I saved 1-prism, 16-ctrp/16 and 1-gdsc in the file above (for ORR comparisons in previous figure))
AUC_ORR_compare$PRISM<-(1-AUC_ORR_compare$PRISM)
AUC_ORR_compare$CTRP<-(1-AUC_ORR_compare$CTRP)
AUC_ORR_compare$GDSC<-(1-AUC_ORR_compare$GDSC)
AVG_responses$PRISM<-(1-AVG_responses$PRISM)
AVG_responses$GDSC<-(1-AVG_responses$GDSC)
AVG_responses$CTRP<-(1-AVG_responses$CTRP)
rm(AVG_responses)



#------------------------------------------------------------------------------
# Add number of patients
#------------------------------------------------------------------------------

load("NCI1970meta_Database.RData")
rownames(NCI1970meta_total)<-tolower(rownames(NCI1970meta_total))

AUC_ORR_compare$patient_num<-NA
for (drug in unique(AUC_ORR_compare$drug)){
  for (cancer in unique(AUC_ORR_compare$cancer)){
    tmp_idx<-intersect(which(AUC_ORR_compare$drug==drug),which(AUC_ORR_compare$cancer==cancer))
    if (drug %in% rownames(NCI1970meta_total) && cancer %in% colnames(NCI1970meta_total)){
      AUC_ORR_compare$patient_num[tmp_idx]<-NCI1970meta_total[drug,cancer]
    }
  }
}

load("clinical_dbs/SEER_cancer_incidence.RData")
AUC_ORR_compare$SEER<-NA
for (cancer in names(SEER_cancer_incidence)){
  tmp_idx<-which(AUC_ORR_compare$cancer==cancer)
  if (length(tmp_idx)>0){
    AUC_ORR_compare$SEER[tmp_idx]<-SEER_cancer_incidence[cancer]
  }
}


#------------------------------------------------------------------------------
# ADD EC50 data
#------------------------------------------------------------------------------
#######################
#LOAD NEEDED DATA
#######################

load("cell-line_dbs/invitro_AUC-lines-individ_dbs.RData")
rm(ctrp_auc_matrix)
rm(gdsc1_auc_matrix)
rm(gdsc2_auc_matrix)
rm(prism_auc_matrix)
hist(prism_ec50_matrix,breaks=100)
prism_ec50_matrix[which(prism_ec50_matrix>100)]<-100
hist(log10(prism_ec50_matrix+1.160622e-07),breaks=100)
prism_logEC50_matrix<-log10(prism_ec50_matrix+1.160622e-07)


load("cell-line_dbs/CCLE_RNAseq.RData")
rm(ccle_RNAseq_filtered)

#######################
#CALC EC50 average per lineage:
#######################

names(sort(table(ccle_lineage_filtered$lineage),decreasing=T))
lineage_to_tcga <- c(
  lung                      = "nsclc",   # or "lusc"
  blood                     = "leuk",
  lymphocyte                = "lymph",
  colorectal                = "coad",   # or "read"
  central_nervous_system    = "gbm",    # or "lgg"
  skin                      = "skcm",
  breast                    = "brca",
  ovary                     = "ov",
  gastric                   = "stad",
  pancreas                  = "paad",
  upper_aerodigestive       = "hnsc",
  uterus                    = "ucec",
  esophagus                 = "esca",
  plasma_cell               = "mm",     # multiple myeloma (not TCGA)
  soft_tissue               = "sarc",
  urinary_tract             = "blca",
  liver                     = "lihc",
  kidney                    = "kirc",   # could be kirc/kirp/kich
  bone                      = "os",     # osteosarcoma (not TCGA)
  peripheral_nervous_system = "pcpg",
  cervix                    = "cesc",
  thyroid                   = "thca",
  bile_duct                 = "chol",
  fibroblast                = "meso",   # closest TCGA-like code
  prostate                  = "prad",
  adrenal_cortex            = "acc"
)

line2lineage<-lineage_to_tcga[ccle_lineage_filtered$lineage]
names(line2lineage)<-rownames(ccle_lineage_filtered)

lineages_ordered<-line2lineage[colnames(prism_logEC50_matrix)]
table(lineages_ordered)
#######################
#POPULATE
#######################
drug="doxorubicin"
lineage="sarc"

AUC_ORR_compare$PRISM_ec50<-NA
for (drug in unique(AUC_ORR_compare$drug)){
  for (lineage in unique(AUC_ORR_compare$cancer)){
  if (drug %in% rownames(prism_logEC50_matrix) && lineage %in% lineages_ordered){
    #EXTRACT MEAN:
    col_idx<-which(lineages_ordered==lineage)
    mean_ec50<-mean(prism_logEC50_matrix[drug,col_idx],na.rm=T)
    
    #ADD TO MASTER COMPARE DF:
    add_idx<-intersect(which(AUC_ORR_compare$drug==drug),which(AUC_ORR_compare$cancer==lineage))
    AUC_ORR_compare$PRISM_ec50[add_idx]<-mean_ec50
  }
  }
}

linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare))
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare))


#------------------------------------------------------------------------------
# ADD Cmax data
#------------------------------------------------------------------------------
tmp_Cmax<-read.csv("lab_v_clinic/Cmax_database-v04.csv",as.is=T,header=T)

rownames(tmp_Cmax)<-make.names(tolower(tmp_Cmax$Genericname),unique=T)

#MAKE NUMERIC:
tmp_Cmax$AUC.ng.hr.mL.<-as.numeric(tmp_Cmax$AUC.ng.hr.mL.)
tmp_Cmax$Cmax.ng.mL.<-as.numeric(tmp_Cmax$Cmax.ng.mL.)
tmp_Cmax$Cmax.μmol.L.<-as.numeric(tmp_Cmax$Cmax.μmol.L.)
tmp_Cmax$Proteinbinding<-(1-as.numeric(tmp_Cmax$Proteinbinding))


#NEW METRICS:
tmp_Cmax$Cmax_unbound<-(tmp_Cmax$Proteinbinding * tmp_Cmax$Cmax.μmol.L.)
tmp_Cmax$AUC_uM<-tmp_Cmax$AUC.ng.hr.mL.* ( tmp_Cmax$Cmax.μmol.L./ tmp_Cmax$Cmax.ng.mL.)
tmp_Cmax$AUC_uM_fu<-(tmp_Cmax$Proteinbinding * tmp_Cmax$AUC_uM)

#POPULATE MASTER DF:
AUC_ORR_compare$Cmax<-NA
AUC_ORR_compare$Cmax_u<-NA
AUC_ORR_compare$AUC<-NA
AUC_ORR_compare$AUC_u<-NA


for (drug in unique(AUC_ORR_compare$drug)){
  if (drug %in% rownames(tmp_Cmax)){
    drug_idx<-which(AUC_ORR_compare$drug==drug)
    
    AUC_ORR_compare$Cmax[drug_idx]<-log10(tmp_Cmax[drug,"Cmax.μmol.L."])
    AUC_ORR_compare$Cmax_u[drug_idx]<-log10(tmp_Cmax[drug,"Cmax_unbound"])
    AUC_ORR_compare$AUC[drug_idx]<-log10(tmp_Cmax[drug,"AUC_uM"])
    AUC_ORR_compare$AUC_u[drug_idx]<-log10(tmp_Cmax[drug,"AUC_uM_fu"])
  }
}

linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$Cmax,rownames(AUC_ORR_compare),
                       title="PRISM vs Cmax",ylab="Log(Cmax uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="PRISM vs Cmax_u",ylab="Log(Cmax_u uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$AUC,rownames(AUC_ORR_compare),
                       title="PRISM vs AUC",ylab="Log(AUC)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       title="PRISM vs AUC_u",ylab="Log(AUC)",xlab="In vitro Response (AUC)")


linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$Cmax,rownames(AUC_ORR_compare),
                       title="PRISM vs Cmax",ylab="Log(Cmax uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="PRISM vs Cmax_u",ylab="Log(Cmax_u uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$AUC,rownames(AUC_ORR_compare),
                       title="PRISM vs AUC",ylab="Log(AUC)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       title="PRISM vs AUC_u",ylab="Log(AUC)",xlab="In vitro Response (AUC)")

linear_fit_plot_pvalue(AUC_ORR_compare$CTRP,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       title="PRISM vs AUC_u",ylab="Log(AUC)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$GDSC,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       title="PRISM vs AUC_u",ylab="Log(AUC)",xlab="In vitro Response (AUC)")

#------------------------------------------------------------------------------
# Add Therapeutic range data
#------------------------------------------------------------------------------


load("lab_v_clinic/MASTER_pk_database_MW.RData")

rownames(PKDB_mw)<-tolower(rownames(PKDB_mw))


#############################################
#PERFORM COMPARISONS
#############################################

PKDB_df<-as.data.frame(PKDB_mw)

#COMPARISON VECTOR:
PKDB_df$TRmin<-log10(PKDB_mw[,"thera_min"]/(1000*PKDB_mw[,"MW"]))
PKDB_df$TRmax<-log10(PKDB_mw[,"thera_max"]/(1000*PKDB_mw[,"MW"]))



###############
#ADD TR-range TO AUC_ORR_compare: Don't add others as not enough data points
###############

AUC_ORR_compare$TRmin<-NA
AUC_ORR_compare$TRmax<-NA


for (drug in unique(AUC_ORR_compare$drug)){
  if (drug %in% rownames(PKDB_df)){
    drug_idx<-which(AUC_ORR_compare$drug==drug)
    
    AUC_ORR_compare$TRmin[drug_idx]<-PKDB_df[drug,"TRmin"]
    AUC_ORR_compare$TRmax[drug_idx]<-PKDB_df[drug,"TRmax"]

  }
}

linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                       title="PRISM vs TRmin",ylab="Log(TRmin uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$TRmax,rownames(AUC_ORR_compare),
                       title="PRISM vs TRmax",ylab="Log(TRmax uM)",xlab="In vitro Response (AUC)")


linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                       title="PRISM vs TRmin",ylab="Log(TRmin uM)",xlab="In vitro Response (AUC)")


save(AUC_ORR_compare,file="AUC_ORR_compare_COMPLETE.RData")

}

####################################################################################################################################################################################################################################################################
#PLOTS:
####################################################################################################################################################################################################################################################################

#------------------------------------------------------------------------------
# RAW LAB VS RAW CLINIC:  no correlation
#------------------------------------------------------------------------------


linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$GDSC,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="GDSC vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$CTRP,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="CTRP vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="PRISM_ec50 vs ORR",ylab="Clinical Response (%)",xlab="In vitro EC50 Log10(uM)")


prism_avg<-tapply(AUC_ORR_compare$PRISM, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
prism_avg<-prism_avg[!is.na(prism_avg)]
orr_avg<-tapply(AUC_ORR_compare$ORR, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
orr_avg<-orr_avg[!is.na(orr_avg)]


#FINAL FIGURE:
linear_fit_plot_ci(AUC_ORR_compare$PRISM,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="PRISM_ec50 vs ORR",ylab="Clinical Response (%)",xlab="In vitro EC50 Log10(uM)")


linear_fit_plot_ci((prism_avg[names(prism_avg)]),orr_avg[names(prism_avg)],names(prism_avg),
                              title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)")

#------------------------------------------------------------------------------
# CLINICAL EXPOSURE:  in vitro AUC correlates best
#------------------------------------------------------------------------------

#PRISM: AUC vs EC50
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="PRISM (EC50) vs Cmax",ylab="Log(Cmax uM)",xlab="In vitro EC50 (log(uM))")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="PRISM (AUC) vs Cmax",ylab="Log(Cmax uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       title="PRISM (AUC) vs AUCu",ylab="Clinical Log(AUC)",xlab="In vitro Response (AUC)")

#PRISM vs GDSC vs CTRP
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="PRISM (AUC) vs Cmax",ylab="Log(Cmax uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$GDSC,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="GDSC (AUC) vs Cmax",ylab="Log(Cmax uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$CTRP,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="CTRP (AUC) vs Cmax",ylab="Log(Cmax uM)",xlab="In vitro Response (AUC)")



prism_avg<-tapply(AUC_ORR_compare$PRISM, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
prism_avg<-prism_avg[!is.na(prism_avg)]
orr_avg<-tapply(AUC_ORR_compare$Cmax_u, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
orr_avg<-orr_avg[!is.na(orr_avg)]


#FINAL FIGURE:
linear_fit_plot_ci(AUC_ORR_compare$PRISM,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       title="PRISM vs Cmax_u",ylab="Log(Cmax_u uM)",xlab="In vitro Response (AUC)")

linear_fit_plot_pvalue((prism_avg[names(prism_avg)]),orr_avg[names(prism_avg)],names(prism_avg),
                       title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)")

linear_fit_plot_ci((prism_avg[names(prism_avg)]),orr_avg[names(prism_avg)],names(prism_avg),
                       title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC)")


#------------------------------------------------------------------------------
# THERAPEUTIC RANGE:
#------------------------------------------------------------------------------
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                       title="PRISM vs TRmin",ylab="Log(TRmin uM)",xlab="In vitro EC50 (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                       title="PRISM vs TRmin",ylab="Log(TRmin uM)",xlab="In vitro Response (AUC)")
linear_fit_plot_pvalue(AUC_ORR_compare$PRISM,AUC_ORR_compare$TRmax,rownames(AUC_ORR_compare),
                       title="PRISM vs TRmax",ylab="Log(TRmax uM)",xlab="In vitro Response (AUC)")


prism_avg<-tapply(AUC_ORR_compare$PRISM, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
prism_avg<-prism_avg[!is.na(prism_avg)]
orr_avg<-tapply(AUC_ORR_compare$TRmin, AUC_ORR_compare$drug,function(x) mean(x, na.rm=T))
orr_avg<-orr_avg[!is.na(orr_avg)]


#FINAL FIGURE:
linear_fit_plot_ci(AUC_ORR_compare$PRISM,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro EC50 Log10(uM)")


linear_fit_plot_ci((prism_avg[names(prism_avg)]),orr_avg[names(prism_avg)],names(prism_avg),
                   title="PRISM vs ORR",ylab="TRmin (log(uM))",xlab="In vitro Response (AUC)")


#------------------------------------------------------------------------------
#FIG S6:  RAW SENSITIVITY ANALYSES
#------------------------------------------------------------------------------
#TEXT:  We performed predefined sensitivity analyses varying datasets (CTRP, GDSC, PRISM), in vitro potency metrics (EC50, AUC), and exposure definitions (Cmax,u, AUC,u, TRmin,u). Results were directionally consistent, with the strength of association increasing systematically as exposure metrics approached therapeutic concentrations.
#######################
#ORR:
#######################
pdf(file="ORR_ec50.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare),
                   ylab="Clinic ORR",xlab="PRISM (EC50)")
dev.off()
pdf(file="ORR_prism.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare),
                   ylab="Clinic ORR",xlab="PRISM (AUC)")
dev.off()
pdf(file="ORR_gdsc.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$GDSC,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare),
                   ylab="Clinic ORR",xlab="GDSC (AUC)")
dev.off()
pdf(file="ORR_ctrp.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$CTRP,AUC_ORR_compare$ORR,rownames(AUC_ORR_compare),
                   ylab="Clinic ORR",xlab="CTRP (AUC)")
dev.off()


#######################
#Cmax_u:
#######################
pdf(file="Cmax_ec50.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                   ylab="Log(Cmax uM)",xlab="PRISM (EC50)")
dev.off()
pdf(file="Cmax_prism.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       ylab="Log(Cmax uM)",xlab="PRISM (AUC)")
dev.off()
pdf(file="Cmax_gdsc.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$GDSC,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       ylab="Log(Cmax uM)",xlab="GDSC (AUC)")
dev.off()
pdf(file="Cmax_ctrp.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$CTRP,AUC_ORR_compare$Cmax_u,rownames(AUC_ORR_compare),
                       ylab="Log(Cmax uM)",xlab="CTRP (AUC)")
dev.off()



#######################
#AUC_u
#######################
pdf(file="AUC_ec50.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                   ylab="Log(Cmax uM)",xlab="PRISM (EC50)")
dev.off()
pdf(file="AUC_prism.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       ylab="Log(Cmax uM)",xlab="PRISM (AUC)")
dev.off()
pdf(file="AUC_gdsc.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$GDSC,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       ylab="Log(Cmax uM)",xlab="GDSC (AUC)")
dev.off()
pdf(file="AUC_ctrp.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$CTRP,AUC_ORR_compare$AUC_u,rownames(AUC_ORR_compare),
                       ylab="Log(Cmax uM)",xlab="CTRP (AUC)")
dev.off()


#######################
#TRmin
#######################
pdf(file="TR_ec50.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM_ec50,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                   ylab="Log(TR uM)",xlab="PRISM (EC50)")
dev.off()
pdf(file="TR_prism.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$PRISM,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                       ylab="Log(TR uM)",xlab="PRISM (AUC)")
dev.off()
pdf(file="TR_gdsc.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$GDSC,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                       ylab="Log(TR uM)",xlab="GDSC (AUC)")
dev.off()
pdf(file="TR_ctrp.pdf",width=3.5,height=3.5)
linear_fit_plot_ci(AUC_ORR_compare$CTRP,AUC_ORR_compare$TRmin,rownames(AUC_ORR_compare),
                       ylab="Log(TR uM)",xlab="CTRP (AUC)")
dev.off()



#------------------------------------------------------------------------------
#FIG S7: normalize AUC by Cmax
#------------------------------------------------------------------------------


####################
#Normalize AUC by Cmax
####################
#NEED TO RE-NORMALIZE LOG-SCALE range to 0-1 (so apples-to-apples w/ AUC)
AUC_ORR_compare$C_norm1<-(AUC_ORR_compare$Cmax_u - min(AUC_ORR_compare$Cmax_u, na.rm = TRUE)) / (max(AUC_ORR_compare$Cmax_u, na.rm = TRUE) - min(AUC_ORR_compare$Cmax_u, na.rm = TRUE))
AUC_ORR_compare$AUC_Cmax<-(AUC_ORR_compare$PRISM-AUC_ORR_compare$C_norm1)


linear_fit_plot_pvalue((AUC_ORR_compare$AUC_Cmax),AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC/Cmax)")

#PLOT w/o gemcitabine:
gem_idx<-which(AUC_ORR_compare$drug=="gemcitabine")
linear_fit_plot_pvalue((AUC_ORR_compare$AUC_Cmax[-gem_idx]),AUC_ORR_compare$ORR[-gem_idx],rownames(AUC_ORR_compare[-gem_idx,]), title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC/Cmax)")

orr_avg <- AUC_ORR_compare %>%
  group_by(drug) %>%
  summarise(mean_val = mean(ORR, na.rm = TRUE)) %>%
  deframe()

Cnorm_avg <- AUC_ORR_compare %>%
  group_by(drug) %>%
  summarise(mean_val = mean(AUC_Cmax, na.rm = TRUE)) %>%
  deframe()

gem_idx<-which(names(orr_avg)=="gemcitabine")
linear_fit_plot_pvalue(Cnorm_avg[-gem_idx],orr_avg[-gem_idx],names(orr_avg)[-gem_idx],
                       title="AUC-Cmax vs Cmax",ylab="Clinical Response",xlab="AUC/Cmax")

#CHECK WEIGHTED VS UNWEIGHTED:
AUCnorm_un<- lm(ORR~AUC_Cmax, data = AUC_ORR_compare) #Fstat 98.28 
summary(AUCnorm_un)

AUCnorm_wt<- lm(ORR~AUC_Cmax, data = AUC_ORR_compare, weights = patient_num) #Fstat 98.28 
summary(AUCnorm_wt)

####################
#Normalize AUC by TRmin
####################
AUC_ORR_compare$TR_norm1<-(AUC_ORR_compare$TRmin - min(AUC_ORR_compare$TRmin, na.rm = TRUE)) / (max(AUC_ORR_compare$TRmin, na.rm = TRUE) - min(AUC_ORR_compare$TRmin, na.rm = TRUE))
AUC_ORR_compare$AUC_TR<-(AUC_ORR_compare$PRISM-AUC_ORR_compare$TR_norm1)


linear_fit_plot_pvalue((AUC_ORR_compare$AUC_TR),AUC_ORR_compare$ORR,rownames(AUC_ORR_compare), title="PRISM vs ORR",ylab="Clinical Response (%)",xlab="In vitro Response (AUC/Cmax)")

orr_avg <- AUC_ORR_compare %>%
  group_by(drug) %>%
  summarise(mean_val = mean(ORR, na.rm = TRUE)) %>%
  deframe()

Cnorm_avg <- AUC_ORR_compare %>%
  group_by(drug) %>%
  summarise(mean_val = mean(AUC_TR, na.rm = TRUE)) %>%
  deframe()

linear_fit_plot_pvalue(Cnorm_avg,orr_avg,names(orr_avg),
                       title="AUC-Cmax vs Cmax",ylab="Clinical Response",xlab="AUC/Cmax")

#CHECK WEIGHTED VS UNWEIGHTED:
AUCnorm_un<- lm(ORR~AUC_TR, data = AUC_ORR_compare) #Fstat 98.28 
summary(AUCnorm_un)

AUCnorm_wt<- lm(ORR~AUC_TR, data = AUC_ORR_compare, weights = patient_num) #Fstat 98.28 
summary(AUCnorm_wt)



####################################################################################################################################################################################################################################################################
#2 CHECK CONSISTENCY WITH PATIENT-number weighting
####################################################################################################################################################################################################################################################################


#------------------------------------------------------------------------------
# FDA Association
#------------------------------------------------------------------------------

#MAKE MODELS
fda_un<- lm(ORR~FDA, data = AUC_ORR_compare) #Fstat 98.28 
fda_wt<-lm(ORR ~ FDA, data = AUC_ORR_compare, weights = patient_num)#Fstat 104.1
summary(fda_un)
summary(fda_wt)

fda_un<- lm(ORR~SEER+FDA, data = AUC_ORR_compare) #Fstat 98.28 
summary(fda_un)

fda_wt<-lm(ORR ~ SEER+FDA, data = AUC_ORR_compare, weights = patient_num)#Fstat 104.1
summary(fda_wt)


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



library(gplots)
library(RColorBrewer)

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


#------------------------------------------------------------------------------
#  LOAD DATA
#------------------------------------------------------------------------------

#Mechanism of Action:   
load("expr_literature/BIOMARKER-signs_matrices.RData")
rm(chemo_BIOMARKERS_matrix)
rm(chemo_BIOMARKERS_matrix_noTargets)
rm(chemo_ddr_matrix)
rm(chemo_metab_matrix)
rm(chemo_target_matrix)
rm(chemo_transport_matrix)
rownames(chemo_BIOMARKERS_matrix_DDR)<-tolower(rownames(chemo_BIOMARKERS_matrix_DDR))

#LOAD LINEAGE AVG DRUGS:
load("NCI1970meta_Database.RData")
rownames(NCI1970meta_rate)<-tolower(rownames(NCI1970meta_rate))
rm(NCI1970meta_respond)
rm(NCI1970meta_total)

#LOAD INDIVIDUAL CELL-LINES:
load("cell-line_dbs/AUC_INDIV_zscore.RData")

load("cell-line_dbs/CCLE_RNAseq.RData")
ccle_RNAseq_zscore<-t(scale(t(ccle_RNAseq_filtered),scale=T,center=T))
rm(ccle_lineage_filtered)


#------------------------------------------------------------------------------
#  CALC CTRP CORR MATRIX:
#------------------------------------------------------------------------------
#synchronize all drugs
overlap_genes<-intersect(colnames(chemo_BIOMARKERS_matrix_DDR),rownames(ccle_RNAseq_zscore))
overlap_lines<-intersect(colnames(ctrp_zscore),colnames(ccle_RNAseq_zscore))
overlap_drugs<-intersect(rownames(chemo_BIOMARKERS_matrix_DDR),rownames(ctrp_zscore))

biomarkers_filtered<-chemo_BIOMARKERS_matrix_DDR[,overlap_genes]

#Calc global correlation matrix
tcga_corr_all<-cor(t(ctrp_zscore[overlap_drugs,overlap_lines]),
                   t(ccle_RNAseq_zscore[overlap_genes,overlap_lines]),use="pairwise.complete.obs")



corr_df<-as.data.frame(matrix(ncol=4,nrow=0))
colnames(corr_df)<-c("drugs","gene","sign","corr")


drug="paclitaxel"
for (drug in overlap_drugs){
  tmp_biomarkers<-biomarkers_filtered[drug,which(abs(biomarkers_filtered[drug,])==1)]
  
  if (length(tmp_biomarkers)>1){
    tmp_corr<-tcga_corr_all[drug,names(tmp_biomarkers)]
    
    tmp_df<-cbind(rep(drug,length(tmp_corr)),
                  names(tmp_biomarkers),
                  tmp_biomarkers,
                  tmp_corr)
    colnames(tmp_df)<-c("drugs","gene","sign","corr")
    
    corr_df<-rbind(corr_df,tmp_df)
  }
  
}


corr_df$final<-as.numeric(corr_df$sign)*as.numeric(corr_df$corr)

library(dplyr)
library(ggplot2)

ctrp_df <- corr_df %>%
  mutate(drugs = reorder(drugs, final, FUN = median, na.rm = TRUE))

rownames(ctrp_df)<-paste0(ctrp_df$drugs,"_",ctrp_df$gene)



ggplot(ctrp_df, aes(x = drugs, y = final)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linewidth = 1, color = "black") +
  coord_flip() +
  theme_bw() +
  xlab("Drug") +
  ylab("Final value")


#ALL MATRIX CORR:
red_vals<-correlation_drugs_fuction(biomarkers_filtered,
                                    ccle_RNAseq_filtered[overlap_genes,overlap_lines],
                                    ctrp_zscore[overlap_drugs,overlap_lines])

red_df <- data.frame(
  drugs = names(red_vals),
  red_point = as.numeric(red_vals)
)

ctrp_df$drugs <- factor(
  ctrp_df$drugs,
  levels = red_df$drugs[order(red_df$red_point)]
)

# also reorder red_df to match
red_df$drugs <- factor(
  red_df$drugs,
  levels = levels(ctrp_df$drugs)
)


ggplot(ctrp_df, aes(x = drugs, y = final)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_hline(yintercept = 0, linewidth = 1, color = "black") +
  geom_point(
    data = red_df,
    aes(x = drugs, y = red_point),
    color = "red",
    size = 3
  ) +
  coord_flip() +
  theme_bw() +
  xlab("Drug") +
  ylab("Final value")





df <- ctrp_df %>%
  mutate(sign = factor(sign))

df$corr<-as.numeric(df$corr)

ggplot(df, aes(x = sign, y = corr)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  xlab("Sign") +
  ylab("Correlation")



#------------------------------------------------------------------------------
#  CALC CLINCIAL CORR MATRIX:
#------------------------------------------------------------------------------
#LOAD LINEAGE AVG RNA:
load("cell-line_dbs/RNAseq_AVG-zscore.RData")
colnames(tcga_avg_matrix)<-gsub("luad","nsclc",gsub("dlbc","lymph",gsub("laml","leuk",colnames(tcga_avg_matrix))))
load("cell-line_dbs/RNAseq_AVG-tpm.RData")
colnames(tcga_avg_tpm)<-gsub("luad","nsclc",gsub("dlbc","lymph",gsub("laml","leuk",colnames(tcga_avg_tpm))))


#synchronize all drugs
overlap_genes<-intersect(colnames(chemo_BIOMARKERS_matrix_DDR),rownames(tcga_avg_tpm))
overlap_lines<-intersect(colnames(NCI1970meta_rate),colnames(tcga_avg_tpm))

biomarkers_filtered<-chemo_BIOMARKERS_matrix_DDR[,overlap_genes]



#Calc global correlation matrix
tcga_corr_all<-cor(t(NCI1970meta_rate[overlap_drugs,overlap_lines]),
                     t(tcga_avg_tpm[overlap_genes,overlap_lines]),use="pairwise.complete.obs")



corr_df<-as.data.frame(matrix(ncol=4,nrow=0))
colnames(corr_df)<-c("drugs","gene","sign","corr")


drug="paclitaxel"
for (drug in overlap_drugs){
  tmp_biomarkers<-biomarkers_filtered[drug,which(abs(biomarkers_filtered[drug,])==1)]
  
  if (length(tmp_biomarkers)>1){
    tmp_corr<-tcga_corr_all[drug,names(tmp_biomarkers)]
    
    tmp_df<-cbind(rep(drug,length(tmp_corr)),
                  names(tmp_biomarkers),
                  tmp_biomarkers,
                  tmp_corr)
    colnames(tmp_df)<-c("drugs","gene","sign","corr")
    
    corr_df<-rbind(corr_df,tmp_df)
  }
  
}


corr_df$final<-as.numeric(corr_df$sign)*as.numeric(corr_df$corr)

library(dplyr)
library(ggplot2)

nci_df <- corr_df %>%
  mutate(drugs = reorder(drugs, final, FUN = median, na.rm = TRUE))

rownames(nci_df)<-paste0(nci_df$drugs,"_",nci_df$gene)

ggplot(nci_df, aes(x = drugs, y = final)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linewidth = 1, color = "black") +
  coord_flip() +
  theme_bw() +
  xlab("Drug") +
  ylab("Final value")


#ALL MATRIX CORR:
red_vals<-correlation_drugs_fuction(biomarkers_filtered,
                                    tcga_avg_tpm[overlap_genes,overlap_lines],
                                    NCI1970meta_rate[overlap_drugs,overlap_lines])

red_df <- data.frame(
  drugs = names(red_vals),
  red_point = as.numeric(red_vals)
)

nci_df$drugs <- factor(
  nci_df$drugs,
  levels = red_df$drugs[order(red_df$red_point)]
)

# also reorder red_df to match
red_df$drugs <- factor(
  red_df$drugs,
  levels = levels(nci_df$drugs)
)


ggplot(nci_df, aes(x = drugs, y = final)) +
  geom_boxplot(fill = "grey", color = "black") +
  geom_hline(yintercept = 0, linewidth = 1, color = "black") +
  geom_point(
    data = red_df,
    aes(x = drugs, y = red_point),
    color = "red",
    size = 3
  ) +
  coord_flip() +
  theme_bw() +
  xlab("Drug") +
  ylab("Final value")




df <- nci_df %>%
  mutate(sign = factor(sign))

df$corr<-as.numeric(df$corr)

ggplot(df, aes(x = sign, y = corr)) +
  geom_boxplot() +
  coord_flip() +
  theme_bw() +
  xlab("Sign") +
  ylab("Correlation")

#------------------------------------------------------------------------------
#  OVERAL BIOMARKER CORRELATION:
#------------------------------------------------------------------------------

overlap_names<-intersect(rownames(nci_df),rownames(ctrp_df))


color<-gsub("1","black",gsub("-1","red",nci_df[overlap_names,"sign"]))

plot(nci_df[overlap_names,"corr"],ctrp_df[overlap_names,"corr"],col=color,pch=16,cex=0.6,
     ylim=c(-0.6,0.6),xlim=c(-0.8,0.8))
abline(h=0,lty=2,lwd=0.5)
abline(v=0,lty=2,lwd=0.5)

#R = .41, p = 2.8e-14
linear_fit_plot_pvalue(as.numeric(nci_df[overlap_names,"corr"]),
                       as.numeric(ctrp_df[overlap_names,"corr"]),overlap_names,pt_colors = color)
abline(a=0,b=1)

#R = .41, p = 2.8e-14
linear_fit_plot_ci(as.numeric(nci_df[overlap_names,"corr"]),
                       as.numeric(ctrp_df[overlap_names,"corr"]),overlap_names,pt_colors = color)
abline(a=0,b=1)


#------------------------------------------------------------------------------
#  ODDS RATIOS:  
#------------------------------------------------------------------------------

lit<-gsub("1","sen",gsub("-1","res",nci_df[overlap_names,"sign"]))

clinic<-nci_df[overlap_names,"corr"]
clinic[which(clinic>0)]<-"sen"
clinic[which(clinic<0)]<-"res"

lab<-ctrp_df[overlap_names,"corr"]
lab[which(lab>0)]<-"sen"
lab[which(lab<0)]<-"res"


#LABORATORY STATISTICS:  OR = 3.7, p = 2.36e-8
lab_cm<-table(lit,lab)
fisher.test(lab_cm)

#CLINIC STATISTICS:  OR = 3.9, p = 6.52e-9
clin_cm<-table(lit,clinic)
fisher.test(clin_cm)


#LAB-CLINIC CORR:  OR = 3.8, p = 1.44e-8
lab_clinic<-table(lab,clinic)
fisher.test(lab_clinic)


#------------------------------------------------------------------------------
#  ODDS RATIOS:   stringent
#------------------------------------------------------------------------------

lit_align<-nci_df[overlap_names,"sign"]
lab_align<-ctrp_df[overlap_names,"corr"]
clinic_align<-nci_df[overlap_names,"corr"]

stringent_idx<-intersect(which(abs(as.numeric(lab_align))>0.1),which(abs(as.numeric(clinic_align))>0.1))

lit_prune<-lit_align[stringent_idx]
lab_prune<-lab_align[stringent_idx]
clinic_prune<-clinic_align[stringent_idx]

lit<-gsub("1","sen",gsub("-1","res",lit_prune))

clinic<-clinic_prune
clinic[which(clinic>0)]<-"sen"
clinic[which(clinic<0)]<-"res"

lab<-lab_prune
lab[which(lab>0)]<-"sen"
lab[which(lab<0)]<-"res"


#LABORATORY STATISTICS:  OR = 13.3
lab_cm<-table(lit,lab)
fisher.test(lab_cm)

#CLINIC STATISTICS:  OR = 8.4
clin_cm<-table(lit,clinic)
fisher.test(clin_cm)


#LAB-CLINIC CORR:  OR = 14.2
lab_clinic<-table(lab,clinic)
fisher.test(lab_clinic)

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


#ALL G2 modelÈ
linear_fit_plot(log10(invert_b1),log10(invert_wt),xlim=c(-2,3),ylim=c(-2,3),
                xlab="Abcg2 (1/b1)",ylab="TKO (1/wt)",title="All G2 Model (1/b1 = 1/wt)")
text(log10(invert_b1),log10(invert_wt),names(invert_sum),cex=0.5)
abline(a=0,b=1,lty=2)



linear_fit_plot(log10(invert_sum),log10(invert_wt),xlim=c(-0.5,3),ylim=c(-0.5,3),
                xlab="Abcb1 (1/g2) + Abcg2 (1/b1)",ylab="TKO (1/wt)",title="Parallel Independent Model (no TKO data) ")
text(log10(invert_sum),log10(invert_wt),names(invert_sum),cex=0.5)







####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################
##        ##          ##  #######   ############      ##                                                       
##        ###        ###  ##             ##          ####                             
##        ####      ####  ##             ##         ##  ##                          
##        ## ##    ## ##  #######        ##        ##    ##                                   
##        ##  ##  ##  ##  ##             ##       ##########                           
##        ##   ####   ##  ##             ##       ##      ##                   
##        ##    ##    ##  #######        ##       ##      ##                                     
####################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################  
#RANDOM EFFECTS META ANALYSES: 

##################
#RESULT TEXT #1:  
##################
#“Across all 540 drug–cancer meta-analyses, the pooled ORR had a median of 19.3% (IQR 10.0%–31.7%). Between-study heterogeneity was low (median τ² = 0; median I² = 0%), reflecting the limited number of historical trials available per drug–cancer pair.”
#We performed 540 independent random-effects meta-analyses of clinical response rates across 30 chemotherapies and 18 cancers. The median number of studies per drug–cancer pair was X (range Y–Z). The pooled ORRs exhibited broad variability (median M%, IQR A–B%), reflecting the heterogeneity of single-agent activity across drugs and diseases. Between-study variability was moderate overall (median I² = H%, IQR H1–H2%). Complete pooled estimates, heterogeneity metrics, and study counts for all drug–cancer pairs are provided in Supplementary Table S1.

##################
#METHOD TEXT:  
##################
#For each drug–cancer pair, we identified all clinical trials reporting single-agent response rates and extracted the number of responders and total evaluable patients. We then performed independent random-effects meta-analyses for each of 540 drug–cancer pairs using a binomial–normal model (logit-transformed proportion, REML estimator). This approach accounts for both within-study sampling variance and between-study heterogeneity (τ²). Zero-event studies were retained without continuity correction using the PLO transformation implemented in the metafor package. For each meta-analysis we reported the pooled ORR, 95% CI, number of contributing studies, τ², and I². Summary statistics across all 540 meta-analyses are shown in Figure X, and full per-pair results are provided in Supplementary Table S1.

library(dplyr)
library(metafor)
library(purrr)
library(readr)

#------------------------------------------------------------
# 1. Load your data
#------------------------------------------------------------
dat <- read_csv("clinical_trials.csv")   # columns: cancer, drug, total_num, respond_num


#------------------------------------------------------------
# 2. Function to run one random-effects meta-analysis
#------------------------------------------------------------
run_meta <- function(df) {
  
  # metafor's escalc() computes the logit-transformed proportion + its sampling variance
  esc <- escalc(
    measure = "PLO",   # logit-transformed proportion (handles zero responders safely)
    xi      = respond_num,
    ni      = total_num,
    data    = df
  )
  
  # random-effects model (REML recommended)
  fit <- tryCatch(
    rma(yi, vi, data = esc, method = "REML"),
    error = function(e) return(NULL)
  )
  
  if (is.null(fit)) {
    return(tibble(
      pooled_logit = NA,
      pooled_orr   = NA,
      ci_lb        = NA,
      ci_ub        = NA,
      tau2         = NA,
      I2           = NA,
      k            = nrow(df)
    ))
  }
  
  tibble(
    pooled_logit = fit$b[1],
    pooled_orr   = transf.ilogit(fit$b),                    # back-transform to ORR
    ci_lb        = transf.ilogit(fit$ci.lb),
    ci_ub        = transf.ilogit(fit$ci.ub),
    tau2         = fit$tau2,
    I2           = fit$I2,
    k            = fit$k
  )
}

#------------------------------------------------------------
# 3. Apply meta-analysis to each cancer–drug pair
#------------------------------------------------------------
results <- dat %>%
  group_by(cancer, drug) %>%
  group_modify(~ run_meta(.x)) %>%
  ungroup()

#------------------------------------------------------------
# 4. Inspect summary results
#------------------------------------------------------------
print(results)

# Optional: summary distributions
results %>%
  summarise(
    median_ORR = median(pooled_orr, na.rm = TRUE),
    IQR_ORR_lb = quantile(pooled_orr, 0.25, na.rm = TRUE),
    IQR_ORR_ub = quantile(pooled_orr, 0.75, na.rm = TRUE),
    median_tau2 = median(tau2, na.rm = TRUE),
    median_I2 = median(I2, na.rm = TRUE)
  )

#------------------------------------------------------------
# 5. Save Supplementary Table
#------------------------------------------------------------
write.csv(results, "Supplementary_Table_S1_pooled_ORR_meta_analysis.csv")


library(ggplot2)

ggplot(results, aes(x = pooled_orr)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_bw() +
  labs(
    title = "Distribution of Pooled Clinical Response Rates",
    x = "Pooled ORR",
    y = "Number of drug–cancer pairs"
  )


library(ggridges)

ggplot(results, aes(x = pooled_orr, y = cancer, fill = cancer)) +
  geom_density_ridges(alpha = 0.7) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(
    title = "Pooled ORR by Cancer Type",
    x = "Pooled ORR",
    y = "Cancer"
  )


ggplot(results, aes(x = pooled_orr, y = I2)) +
  geom_point(alpha = 0.7, color = "darkred") +
  theme_bw() +
  labs(
    title = "Heterogeneity vs. Clinical Response Rate",
    x = "Pooled ORR",
    y = "I² (%)"
  )


ggplot(results, aes(x = drug, y = cancer, fill = pooled_orr)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(
    title = "Heatmap of Pooled ORR Across Drug–Cancer Pairs",
    x = "Drug",
    y = "Cancer",
    fill = "Pooled ORR"
  )









