#Breast Cancer Detection Project
library(tidyverse)

##### DATA CLEANING #####

#read in data set 
setwd("C:/Users/17739/Desktop/ST 537 Applied Multivariate and Longitudinal Data Analysis/R Project/Breast Cancer Git Hub/Data/")
bc.df = read.table('breast-cancer-wisconsin.data', sep=',')
colnames(bc.df) = c('Id', 'Clump Thickness', 'Cell Size Uniformity', 'Cell Shape Uniformity', 
                    'Marginal Adhesion', 'Single Epithelial Cell Size', 'Bare Nuclei', 'Bland Chromatin', 
                    'Normal Nucleoli', 'Mitosis', 'Class')

head(bc.df) #view first 6 rows of data
dim(bc.df) #check dimensions 
str(bc.df) #Bare Nuclei is only feature stored as chr instead of int
table(bc.df$Class) #frequency for each class

#write dataframe to csv 
#write.csv(bc.df,"breast-cancer-wisconsin.csv", row.names = FALSE)

#identify missing data
bc.df%>%filter(!`Bare Nuclei` %in% c("1","2","3","4","5","6","7","8","9","10"))

#calculte mean Bare Nuclei by Class 
mean.BareNuclei.2<-bc.df%>%group_by(Class)%>%summarize(Mean = mean(as.numeric(`Bare Nuclei`), na.rm=TRUE))%>%filter(Class == 2)
mean.BareNuclei.4<-bc.df%>%group_by(Class)%>%summarize(Mean = mean(as.numeric(`Bare Nuclei`), na.rm=TRUE))%>%filter(Class == 4)

#impute missing with mean
#Bare Nuclei (class 2) with mean Bare Nuclei for class 2
#Bare Nuclei (class 4) with mean Bare Nuclei for class 4
bc.df<-bc.df%>%mutate(`Bare Nuclei Revised` = ifelse(`Bare Nuclei` == "?" & Class == 2, round(mean.BareNuclei.2$Mean,2), 
                                                    ifelse(`Bare Nuclei` == "?" & Class == 4, round(mean.BareNuclei.4$Mean,2), `Bare Nuclei`))) #%>%filter(`Bare Nuclei` == "?")

#convert Bare Nuclei Revised to int
bc.df$`Bare Nuclei Revised`<-as.integer(bc.df$`Bare Nuclei Revised`)

str(bc.df)

#write revised dataframe to csv 
#write.csv(bc.df,"breast-cancer-wisconsin_revised.csv", row.names = FALSE)


##### PCA (Principle Components Analysis #####
#Scientific Question 2
#Can we use PCA to reduce the dimensionality of the data and to identify/summarize crucial variables?

library(knitr)
library(GGally)


bc.df.class<-bc.df$Class
bc.df.X<-select(bc.df, -c("Id","Bare Nuclei","Class"))

#correlation matrix
ggcorr(bc.df.X, label =  T, label_size = 3, label_round = 2)

cor(bc.df.X)
#shows that Cell Shape Uniformity and Cell Size Uniformity are highly correlated (0.91)

#check sd for each variable 
round(apply(bc.df.X, 2, sd),3) # need to standardize variables 

#standardize variables 
bc.df.X.std<-scale(bc.df.X, center = T, scale = T)
#check if scaled correctly 
apply(bc.df.X.std, 2,sd)

#perform PCA
data.pca <- prcomp(bc.df.X.std)
#Extract the importance of each component 
summary(data.pca)

#how many PCs to retain? 
screeplot(data.pca, type = "line", main = "Variance Explained by each PC")
#the elbow in the graph is at PC2
#PC 1 and 2 explain 74% of the variability

#review loadings for PC1 and PC2
round(data.pca$rotation[,1:2],3)
#PC1 is roughly equally weighted 
#PC2 is essentially the effect of Mitosis
PC1<-data.pca$rotation[,1]
PC2<-data.pca$rotation[,2]

PC.scores.1<-bc.df.X.std%*%PC1
PC.scores.2<-bc.df.X.std%*%PC2

#plot PC score 1 vs PC score 2 colored by class
plot(PC.scores.1,PC.scores.2, pch = 19, col = bc.df.class)