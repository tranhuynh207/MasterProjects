rm(list = ls())
#Survey Methods Project Work
#Data on Migration
#Source OECD

#Variables
#Y=Nbr of Work Migrant
#X1=Average Annual Wage 
#X2=Poverty Rate Working Pop
#X3=Unemployement Rate Foreign
#X4=Working Pop
#X5=Income Inequality (GINI)

#Read the dataset

setwd("C:/Users/Martial Houessou/OneDrive - UGent/Academic/IMRD_2018_2020/Semester2/Survey_Methods/Project_Work/Migration/Data")

ItMig<-read.csv("Dataset_Mig1.csv", header=T, sep=",")
ItMig

#Descriptive Statiistics
dim(ItMig)
str(ItMig)
summary(ItMig)
plot(Work_Migrt_Nbr~Year,data=ItMig)

# X variables correlation with Y
plot(Work_Migrt_Nbr~Av_An_Wage_Euro + Pov_Rt_WP + Unemp_Rt + Wking_Pop + Income_Ineq, data = ItMig)
Corr1<-cor(ItMig$Work_Migrt_Nbr,ItMig$Av_An_Wage_Euro)
Corr2<-cor(ItMig$Work_Migrt_Nbr,ItMig$Pov_Rt_WP)
Corr3<-cor(ItMig$Work_Migrt_Nbr,ItMig$Unemp_Rt)
Corr4<-cor(ItMig$Work_Migrt_Nbr,ItMig$Wking_Pop)
Corr5<-cor(ItMig$Work_Migrt_Nbr,ItMig$Income_Ineq)

#Simple Linear Regression Model
MigLM<-lm(Work_Migrt_Nbr~Av_An_Wage_Euro + Unemp_Rt + Income_Ineq, data = ItMig)
MigLM
summary(MigLM)

#Evaluating the model validity
par(mfrow=c(2,2))
plot(MigLM)
par(mfrow=c(1,1))


#Inference
#It can be infered from this analysis that the variables "Average Annual Wage in Euro" and "Income Inequality" satistically affect significantly 
#the variable Y which is number of work migrant in Italy. Variables such as Poverty rate of working population, Working pop and Unemployement rate
#seem to not have any significant effect on the number of migrant to Italy. 
#An increase in the annual wage significantly increase the migration decision of work migrant from abroas hence an increase of work migrant in Italy
#Furthermore, it is noticed a reduction of income inequality also result an increase of work migrant. An equal income economy seems to motive migrant 
#to move from their country of origin to Italy.
#to be statistically significant