###########################################################################
############### PROJECT STATISTICAL LEARNING APPENDIX 2 ##################
##########################################################################


######################## DATA CLEANSING ##############################

#datasets downloaded from https://www.kaggle.com/nxpnsv/country-health-indicators
#import dataset
country_health = read.csv('country_health.csv')
dim(country_health) #180 70

#select only a subset of the most relevant countries
subrow = c(7,9,10,17,21,24,31,32,34,37,39,40,44,45,47,
           51,61,62,66,68,74,78,79,84,85,87,90,91,99,100,
           105,108,111,116,120,121,130,138,
           145,153,154,155,159,163,168,169,172,173,177)
country_health = country_health[subrow,]

#delete the column of the countries' names and save it
country_name = country_health$Country_Region
country_health = country_health[,-1]

#delete the variables for which no information is provided 
subcol = 8:61
country_health = country_health[,subcol]

#fix the names of the columns 
colnames(country_health)[1] = c('Cardiovascular.diseases %')
colnames(country_health)[2] = c('Cancers %')
colnames(country_health)[3] = c('Diabetes, blood and endocrine.diseases %')
colnames(country_health)[4] = c('Respiratory.diseases %')
colnames(country_health)[5] = c('Liver.disease %')
colnames(country_health)[6] = c('Diarrhea and common.infectious.disease %')
colnames(country_health)[7] = c('Musculoskeletal.disorders %')
colnames(country_health)[8] = c('HIV.AIDS.and.tuberculosis %')
colnames(country_health)[9] = c('Malaria and neglected.tropical.diseases %')
colnames(country_health)[10] = c('Nutritional.deficiencies %')
colnames(country_health)[12] = c('Share.of.deaths.from.smoking %')
colnames(country_health)[54] = c('obesity - adult.prevalence.rate')

#fix the data type
country_health = apply(country_health, 2, as.numeric)

#rename rows with countries' names
rownames(country_health) = country_name

#look at the number of null values for each column
for(i in seq(1:49)) {
  s = sum(is.na(country_health[,i]))/length(country_health[,i])
  cat('Percentage of Nan values in', colnames(country_health)[i], ':', s*100,'\n')}

#delete the columns containing too many null values
no = c(23,36,37,38,39,41,42,43)
country_health = country_health[,-no]

#remove observations with at least a null value
country_health_clean = country_health[complete.cases(country_health),] 
dim(country_health_clean) #44 46


#check for outliers
summary(country_health_clean)
boxplot(country_health_clean, names=1:46, col=rainbow(10,alpha=0.5),
        outcol='red', pch=8)
#outliers are present but robust methods in order to avoid 
#outliers to affect PCA instead of deleting observations 


########################### PCA ###########################

#standardize the dataset by subtracting the mean and dividing
#by the standard deviation
country_health_clean = scale(country_health_clean)

library(MASS)
#use a random subsample of 15 features

set.seed(123456789, kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
s=sample(1:46, 15)
country_health_clean2 = country_health_clean[,s]
colnames(country_health_clean2)

#assign row and column number
n = nrow(country_health_clean2)
p = ncol(country_health_clean2)

#calculate eigenvectors and eigenvalues starting from the
#robust correrlation matrix
rho_rob = cov.rob(country_health_clean2, cor=T)
rho = rho_rob$cor
autoval = eigen(rho)$values
autovec = eigen(rho)$vectors

#look at the amount of explained variance
autoval
#look at the percentage and cumulating percentage of 
#explained variance 
expl_pvar = autoval/p
expl_pcum = cumsum(expl_pvar)
expl_var = cbind(expl_pvar, expl_pcum)
expl_var


#use the screeplot to see graphically the explained variance 
#and select the optimal number of components which are above 1
plot(x=c(1:15) ,y=autoval, type='o', main = 'Scree diagram', 
     xlab = 'Principal components', ylab='eigenvalues', col='blue', lwd=2)
abline(h=1, col='red',lty=2,lwd=2)
axis(1, seq(0,15,1))

#plot the cumulating variance
plot(cumsum(expl_var[,1]*100), type='o', ylab='Cumulating PVE', xlab = 'Principal components', col='blue', lwd=2)
axis(1, seq(0,15,1))
axis(2, seq(0,100,5))
abline(h=74, col='red',lty=2,lwd=2)
#select the first 3 components


#calculate the loadings and build the components matrix
#for every principal component
comp_matrix = round(cbind(-eigen(rho)$vectors[,1]*sqrt(autoval[1]),
                          -eigen(rho)$vectors[,2]*sqrt(autoval[2]),
                          -eigen(rho)$vectors[,4]*sqrt(autoval[3])),
                    3)
rownames(comp_matrix) = colnames(country_health_clean2)
colnames(comp_matrix) = c("PC1","PC2","PC3")
comp_matrix

#compute communality and bind it to components matrix
communality = comp_matrix[,1]^2 + comp_matrix[,2]^2 + comp_matrix[,3]^2
comp = cbind(comp_matrix,communality)
comp 


#plot the LOADINGS PLOTS, which show the relation between
#features and PC
library(ggplot2)
library(ggrepel)

c_names = colnames(country_health_clean2)

#loading plot for the 1st and 2nd components
ggplot(data.frame(cbind(comp[,1],comp[,2])), aes(x = comp[,1], y = comp[,2], label = c_names)) +
  geom_point(col='blue')+
  geom_text_repel(size=5)+
  labs(title='LOADING PLOT',subtitle='PC1 vs.PC2', x='PC1', y='PC2')+
  theme_gray(base_size = 10)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), plot.subtitle= element_text(size = 10, face = "bold", hjust = 0.5))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")+
  xlim(-1,1)+
  ylim(-1,1)

#loading plot for the 2nd and 3rd components
ggplot(data.frame(cbind(comp[,2],comp[,3])), aes(x = comp[,2], y = comp[,3], label = c_names)) +
  geom_point(col='purple')+
  geom_text_repel(size=5)+
  labs(title='LOADING PLOT',subtitle='PC2 vs.PC3', x='PC2', y='PC3')+
  theme_gray(base_size = 10)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), plot.subtitle= element_text(size = 10, face = "bold", hjust = 0.5))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")+
  xlim(-1,1)+
  ylim(-1,1)


#calculate the scores
z_country_health_clean2 = scale(country_health_clean2, center=T, scale=T)
zscore = z_country_health_clean2 %*% autovec [,1:3] 

scores = round(cbind(-zscore[,1]/sqrt(autoval[1]),
                     -zscore[,2]/sqrt(autoval[2]),
                     -zscore[,3]/sqrt(autoval[3])
),2)

r_names = rownames(country_health_clean2)

#plot the SCORES PLOT, which shows the relation between 
#observations and principal components

#score plot of the 1st and 2nd components
ggplot(data.frame(cbind(scores[,1],scores[,2])), aes(x = scores[,1], y = scores[,2], label = r_names)) +
  geom_point(col='blue')+
  geom_text_repel(size=5)+
  labs(title='SCORES PLOT',subtitle='PC1 vs.PC2', x='PC1', y='PC2')+
  theme_gray(base_size = 10)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), plot.subtitle= element_text(size = 15, face = "bold", hjust = 0.5))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")

#score plot of the 2nd and 3rd components
ggplot(data.frame(cbind(scores[,2],scores[,3])), aes(x = scores[,2], y = scores[,3], label = r_names)) +
  geom_point(col='purple')+
  geom_text_repel(size=5)+
  labs(title='SCORES PLOT',subtitle='PC2 vs.PC3', x='PC2', y='PC3')+
  theme_gray(base_size = 10)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), plot.subtitle= element_text(size = 15, face = "bold", hjust = 0.5))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")


#create a princomp object using a robust correlation matrix
#in order to plot a biplot
pr = princomp(country_health_clean2, cor=T, covmat = rho_rob)

library(ggfortify)

#biplot of 1st and 2nd components
autoplot(pr, x=1, y=2, label = TRUE, label.size = 3, loadings=TRUE,
         loadings.label = TRUE, loadings.label.size = 3,
         shape = FALSE, loadings.label.vjust = 1.2, label.repel=T)+
  labs(title='BIPLOT',subtitle='PC1 vs.PC2')+
  theme_gray(base_size = 10)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), plot.subtitle= element_text(size = 15, face = "bold", hjust = 0.5))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")

#biplot of 2nd and 3rd components
autoplot(pr, x=2, y=3, label = TRUE, label.size = 3, loadings=TRUE,
         loadings.label = TRUE, loadings.label.size = 3,
         shape = FALSE, loadings.label.vjust = 1.2, label.repel=T)+
  labs(title='BIPLOT',subtitle='PC2 vs.PC3')+
  theme_gray(base_size = 10)+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), plot.subtitle= element_text(size = 15, face = "bold", hjust = 0.5))+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")


###################### K-MEANS CLUSTERING ######################

library(rgl)

colnames(scores) = c('PC1','PC2','PC3')
scores = data.frame(scores)
loadings = data.frame(comp)

#look at how within-variance decreases as the number of clusters
#increase
set.seed(1,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding'))
wss = (nrow(scores))*sum(apply(scores,2,var))
for (i in 2:10) wss[i] = sum(kmeans(scores, 
                                    centers=i)$withinss)
plot(1:10, wss, type="b", xlab="Number of clusters",
     ylab="Within Deviance", col='blue', lwd=2)
axis(1, seq(0,10,1))
axis(2, seq(0,120,10))
abline(h=41, col='red',lty=2,lwd=2)
abline(h=22, col='red',lty=2,lwd=2)


set.seed(98,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
#fit the k-means with 4 clusters starting from 50 different
#starting points
knn4 = kmeans(scores, 4, nstart=50) # 4 cluster solution

#mean of each cluster
mean4 = aggregate(scores,by=list(knn4$cluster),FUN=mean)
#matrix containing the observations and their corresponding
#cluster 
data_knn_4 = data.frame(scores, knn4$cluster)
#3-dimensional plot of 4 clusters k-means
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(4)[knn4$cluster], size=1.5, type='s', main='k-means clustering (4 clusters)') 


set.seed(123,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
#fit the k-means with 7 clusters starting from 50 different
#starting points
knn7 = kmeans(scores, 7, nstart=50) 

#matrix containing the observations and their corresponding
#cluster 
data_knn_7 = data.frame(scores, knn7$cluster)
#3-dimensional plot of 7 clusters k-means
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(7)[knn7$cluster], size=1.5, type='s',main='k-means clustering (7 clusters)') 


#means of 4 clusters
medie1 = aggregate(scores, list(knn4$cluster), mean)

#plot of the means of clusters for every principal component
library(ggplot2)
library(gridExtra)

m1_pc1 = data.frame(mean=medie1$PC1,cluster=1:4)
m1_pc2 = data.frame(mean=medie1$PC2,cluster=1:4)
m1_pc3 = data.frame(mean=medie1$PC3,cluster=1:4)


p1_1 = ggplot(m1_pc1, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#85C1E9')+
  labs(title='Principal Component 1')+
  theme_minimal()

p2_1 = ggplot(m1_pc2, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#F1948A')+
  labs(title='Principal Component 2')+
  theme_minimal()

p3_1 = ggplot(m1_pc3, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#D2B4DE', size=1)+
  labs(title='Principal Component 3')+
  theme_minimal()

grid.arrange(p1_1, p2_1, p3_1, nrow=1)


################## HIERARCHICAL CLUSTERING ##################

#compute the euclidian distance between scores
d1 = dist(scores, method='euclidian')


# fit hierarchical clustering with SINGLE LINKAGE
set.seed(123,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
h5 = hclust(d1, method='single')

#plot dendogram
plot(h5, main="Single (euclidian)")
abline(h=1.1, col='red',lty=2,lwd=2) 

#divide in 5 clusters
cluster5=cutree(h5, k=5)
#3-dimensional plot
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(5)[cluster5], size=1.5, type='s', main='SINGLE (EUCLIDIAN)')


# fit hierarchical clustering with CENTROID LINKAGE
set.seed(123,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
h4 = hclust(d1, method='centroid')
#plot dendogram
plot(h4, main="Centroid (euclidian)")
abline(h=1.25, col='red',lty=2,lwd=2) 

#divide in 5 clusters
cluster4=cutree(h4, k=5)
#3-dimensional plot
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(5)[cluster4], size=1.5, type='s', main='CENTROID (EUCLIDIAN)')


# fit hierarchical clustering with AVERAGE LINKAGE
set.seed(123,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
h1 = hclust(d1, method='average')

#plot the dendogram
plot(h1, main="Average (euclidian)")
abline(h=2, col='red',lty=2,lwd=2) #4 clusters
abline(h=2.75, col='blue',lty=2,lwd=2) #7 clusters

#divide in 4 clusters
cluster1_2=cutree(h1, k=4)
#3-dimensional plot
library(rgl)
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(4)[cluster1_2], size=1.5, type='s', main='AVERAGE (EUCLIDIAN)')

#divide in 7 clusters
cluster1=cutree(h1, k=7)
#3-dimensional plot
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(7)[cluster1], size=1.5, type='s', main='AVERAGE (EUCLIDIAN)')


# fit hierarchical clustering with COMPLETE LINKAGE 
set.seed(123,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
h2 = hclust(d1, method='complete')

#plot dendogram
plot(h2, main="Complete (euclidian)")
abline(h=3.7, col='red',lty=2,lwd=2) #4 clusters

#divide in 4 clusters
cluster2=cutree(h2, k=4)
#3-dimensional plot
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(4)[cluster2], size=1.5, type='s', main='COMPLETE (EUCLIDIAN)')

#means of clusters
medie2 = aggregate(scores, list(cluster2), mean)

#plot of the means of clsusters for every principal component
library(ggplot2)
library(gridExtra)

m2_pc1 = data.frame(mean=medie2$PC1,cluster=1:4)
m2_pc2 = data.frame(mean=medie2$PC2,cluster=1:4)
m2_pc3 = data.frame(mean=medie2$PC3,cluster=1:4)


p1_2 = ggplot(m2_pc1, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#85C1E9')+
  labs(title='Principal Component 1')+
  theme_minimal()

p2_2 = ggplot(m2_pc2, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#F1948A')+
  labs(title='Principal Component 2')+
  theme_minimal()

p3_2 = ggplot(m2_pc3, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#D2B4DE', size=1)+
  labs(title='Principal Component 3')+
  theme_minimal()

grid.arrange(p1_2, p2_2, p3_2, nrow=1)



# fit hierarchical clustering with WARD'S LINKAGE
set.seed(123,kind='Mersenne-Twister',
         normal.kind='Inversion', sample.kind='Rounding')
h3 = hclust(d1, method="ward.D2")

#plot dendogram
plot(h3, main="Ward linkage (euclidian)")
abline(h=3.5, col='red',lty=2,lwd=2) #5 clusters
abline(h=5, col='blue',lty=2,lwd=2) #4 clusters

#divide in 5 clusters
cluster3=cutree(h3, k=5) 
#3-dimensional plot
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(5)[cluster3], size=1.5, type='s', main='WARD LINKAGE (EUCLIDIAN)')

#divide in 4 clusters
cluster3_1=cutree(h3, k=4) 
#3-dimensional plot
plot3d(scores[,1], scores[,2], scores[,3], col = rainbow(4)[cluster3_1], size=1.5, type='s', main='WARD LINKAGE (EUCLIDIAN)')


#means of 4 clusters
medie3 = aggregate(scores, list(cluster3_1), mean)

#plot of the means of clusters for every principal component
library(ggplot2)
library(gridExtra)

m3_pc1 = data.frame(mean=medie3$PC1,cluster=1:4)
m3_pc2 = data.frame(mean=medie3$PC2,cluster=1:4)
m3_pc3 = data.frame(mean=medie3$PC3,cluster=1:4)


p1_3 = ggplot(m3_pc1, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#85C1E9')+
  labs(title='Principal Component 1')+
  theme_minimal()

p2_3 = ggplot(m3_pc2, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#F1948A')+
  labs(title='Principal Component 2')+
  theme_minimal()

p3_3 = ggplot(m3_pc3, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#D2B4DE', size=1)+
  labs(title='Principal Component 3')+
  theme_minimal()

grid.arrange(p1_3, p2_3, p3_3, nrow=1)


#means of 5 clusters
medie4 = aggregate(scores, list(cluster3), mean)

#plot of the means of clusters for every principal component
library(ggplot2)
library(gridExtra)

m4_pc1 = data.frame(mean=medie4$PC1,cluster=1:5)
m4_pc2 = data.frame(mean=medie4$PC2,cluster=1:5)
m4_pc3 = data.frame(mean=medie4$PC3,cluster=1:5)


p1_4 = ggplot(m4_pc1, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#85C1E9')+
  labs(title='Principal Component 1')+
  theme_minimal()

p2_4 = ggplot(m4_pc2, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#F1948A')+
  labs(title='Principal Component 2')+
  theme_minimal()

p3_4 = ggplot(m4_pc3, aes(x=cluster, y=mean)) +
  geom_bar(stat='identity', fill='#D2B4DE', size=1)+
  labs(title='Principal Component 3')+
  theme_minimal()

grid.arrange(p1_4, p2_4, p3_4, nrow=1)


#look at how the clusters have been aggregated at each step
agglo = function(hc){
  data.frame(row.names=paste0("Cluster",seq_along(hc$height)),
             height=hc$height,
             components=ifelse(hc$merge<0, 
                               hc$labels[abs(hc$merge)], paste0("Cluster",hc$merge)),
             stringsAsFactors=FALSE) }

agglo(h1)
agglo(h2)
agglo(h3)
agglo(h4)
agglo(h5)

