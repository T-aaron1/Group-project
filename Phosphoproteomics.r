#R version 3.6.0
install.packages("ggplot2")
library(ggplot2)

sampledata<-read.table("az20.txt", sep="\t", header = TRUE)

#initial volcano plot with the sample data
#note alot of insignificant p values
#added x=10 straight line 
plot( y=-log(sampledata$AZ20_p.value),x=log(sampledata$AZ20_fold_change))
abline(h=10)

#assesing the level of variance between the control and treatment groups
#histograms show a similar positive skew of the variance for samples
hist(sampledata$AZ20_ctrlCV)
hist(sampledata$AZ20_treatCV)

?p.adjust
#intial adjustment of p values
p.adjust(sampledata$AZ20_p.value, method= "fdr" , n = length(sampledata$AZ20_p.value))

#adjustment of pvalues not needed as each phosphorylation test is not needed 
#each test is independent despite them coming from the same sample

#each substrate is a phosphosite and the first column is area under the ms1 peak for that phosposite
#the second column is that are when the inhibitor az20 is introduced
#fold change is the caluclation ratio value from the initial value to the final value(second)
#when the fold change is less than 1 the area under the peak has been reduced 
#the activity could be indicated by the fold change, reasearch whether an increase or decrease of fold change is relative to kinase activity and how

