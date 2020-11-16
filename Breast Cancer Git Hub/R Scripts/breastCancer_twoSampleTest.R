## Author: Rong Huang
## Date: 11/15/2020
## Purpose: st537 final project-Breast Cancer Analysis, two sample test


## read in data set, data processing 
data = read.csv("breast-cancer-wisconsin_revised.csv", header=TRUE)[,-7]
data$Class = as.factor(data$Class)
str(data)

## data exploratory
# number of records for each class: Pie Chart with Percentages
slices = c(table(data$Class)[1], table(data$Class)[2])
pct = round(slices/sum(slices)*100)
lbls = paste(table(data$Class), c("Benign", "Malignant"), paste("(", pct, "%)", sep=""))
pie(slices,labels = lbls, col=c(4, 2), main="Pie Chart of Classes")

# pair plot
pairs(data[,-c(1,10)],
      col = c(4, 2)[as.numeric(data$Class)],   # Change color by group
      pch = c(8, 1)[as.numeric(data$Class)],   # Change points by group
      main = "Pair Plot of Benign and Malignant for All Variables", oma=c(2,2,2,16))
par(xpd = TRUE)
legend("bottomright", col=c(4, 2),pch = c(8, 1),
       legend = c("Benign", "Malignant"))

# boxplot
par(mfrow=c(3,3))
for (name in names(data)[-c(1,10)]){
  boxplot(data[,name]~ Class, data = data, col=c(4,2),
          xlab='', ylab = '', main=name)
}

## data analysis
library(ICSNP)
X = as.matrix(data[data$Class==2,-c(1,10)]) # class 2
Y = as.matrix(data[data$Class==4,-c(1,10)]) # class 4

#normality check
################## chisquare.plot ##############
chisquare.plot <- function(x, mark, title) {
  # x= a n x p data matrix, mark: number of
  # extreme points to mark
  # p=number of variables, n=sample size
  p <- ncol(x)
  n <- nrow(x)
  # xbar and s
  xbar <- colMeans(x)
  s <- cov(x)
  # Mahalanobis dist, sorted
  x.cen <- scale(x, center = T, scale = F)
  d2 <- diag(x.cen %*% solve(s) %*% t(x.cen))
  sortd <- sort(d2)
  # chi-sq quantiles
  qchi <- qchisq((1:n - 0.5)/n, df = p)
  # plot, mark points with heighest distance
  plot(qchi, sortd, pch = 19, xlab = "Chi-square quantiles",
       ylab = "Mahalanobis squared distances",
       main = title)
  points(qchi[(n - mark + 1):n], sortd[(n - mark + 1):n], cex = 3, col = "#990000")
}
###########################################
par(mfrow=c(1,2))
chisquare.plot(X,2, "Chi-square Q-Q Plot of Benign")
chisquare.plot(Y,2, "Chi-square Q-Q Plot of Malignant")

## Two-sample Hotelling's T2 test
HotellingsT2(X,Y)

## Bonferroni intervals
alpha = 0.05 # old significance level
p = 2 # number of intervals/variables
# mean vectors
xbar = colMeans(X)
ybar = colMeans(Y)
difference = xbar - ybar
# covariances of each group
S.x = cov(X)
S.y = cov(Y)
# data statistics summary
stats = round(cbind(xbar, ybar, sqrt(diag(S.x)), sqrt(diag(S.y))),3)
colnames(stats) = c('Benign Mean', 'Malignant Mean', 'Benign Sd', 'Malignant Sd')
stats
# sample sizes
m = nrow(X)
n = nrow(Y)
# pooled covariance matrix
S.pool = ((m-1)*S.x + (n-1)*S.y) / (m + n - 2)
# critical value
crit = qt(alpha/(2*p), df = m+n-2, lower.tail = F)
# Bonferroni intervals
half.width = crit*sqrt(diag(S.pool)*(1/m + 1/n))
lower = difference - half.width
upper = difference + half.width
int.bonf = cbind(lower, upper)
int.bonf


