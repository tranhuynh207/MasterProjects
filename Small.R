samp<-read.csv("sample_group1.csv", header=T, sep=";")
cens<-read.csv("CENSUS_data_R.csv", header=T, sep=";")
str(samp)
str(cens)

#match cens and samp

N.cens <- subset(cens, COD.REGIO.SILC %in% samp$REGIO) 

#create variables for linear regression
samp$AGE_class1<-ifelse(samp$age<25, 1, 0)
samp$PROPERTY.HOUSE<-ifelse(samp$HH020==1,1,0 )
samp$THREE.MEMBERS<-ifelse(samp$hhsize==3,1,0)
samp$MALE<-ifelse(samp$RB090==1,1,0)

#Ordering sample and pop

samp<-samp[order(samp [,"PROV"]),]
N.cens<-N.cens[order(N.cens[,"COD.REGIO.SILC"]),]

View(table(samp$PROV))

library(sae)

#Direct estimation

HT <- direct(eqhhincome, dom=PROV,
             sweight=DB090,
             domsize=aggregate(samp$DB090, list(samp$PROV), sum),
             data=samp)
summary(HT$CV)

# categorize areas on the basis of CV
install.packages("carData")
library(car)
cate_CV<-recode(HT$CV, "0:16.5= '< 16.5%'; 16.51:33.3='between 16.6% and 33.3% '; else='>33.3%'")
table(cate_CV)

# map the results
install.packages("sp")
library(maptools)
library(RColorBrewer)
library(classInt)
library(sp)
library(lattice)

list.files()
ita=readShapePoly("Prov2011_g_WGS84.shp")

plot(ita)

# merge results and shape file
ita.merge<- merge(ita, HT, by.x= "COD_PROV", by.y="Domain")

spplot(ita.merge, "Direct", main = "Direct Estimate NUTS 3")
# number of classes
nclr <- 4
class <- classIntervals(ita.merge$Direct, nclr, style="quantile")

categ = class[["brks"]]

display.brewer.all()

plotclr <-brewer.pal(nclr,"Dark2")

plotclr2 <- c("gray72","gray40", "gray17", "gray7")

ita.merge$cutDIR <- cut(ita.merge$Direct, breaks = categ, 
                  right = F, include.lowest=T)  

# plot results with spplot
spplot(ita.merge, "cutDIR", col.regions=plotclr, 
       main = "Direct Estimates (NUTS 3)",
       colorkey = list(space="bottom",height = 0.9, 
                       labels = list(at = seq(0.5, length(categ) -0.5), 
                                     labels = round(categ,2))))


# unit level model
# linear model to check the significance of covariate
lm01 <- lm(eqhhincome~ AGE_class1 + MALE + THREE.MEMBERS + PROPERTY.HOUSE, data=samp)
summary(lm01)


# look for significant covariates with a stepwise method
slm01 <- step(lm01, direction="both",trace=0)
summary(slm01)


# population size
censN<-aggregate(samp$DB090, list(samp$PROV), sum)
N.cens1<-N.cens[,c(4,5,13,17,22)]

# EBLUP
unit.eblup<-pbmseBHF(eqhhincome~ AGE_class1 + MALE + THREE.MEMBERS + PROPERTY.HOUSE,
                      dom=PROV, meanxpop=N.cens1, 
                      popnsize=censN, B=500, data=samp) 
# 31 areas with zero sample size
N.cens2<- cens[,c(4,5,13,17,22)]
View(N.cens2)
setd<-setdiff(N.cens2$COD.PROV.SILC,samp$PROV)
setd
length(setd)

# select popolation mean for out of sample areas
X.m_tmp <- subset(N.cens2, COD.PROV.SILC %in% setd) 
X.m_o<-X.m_tmp[,-1]

# compute synthethic estimator
beta.hat<-unit.eblup$est$fit$summary$coefficients[,1]
betaest<-matrix(beta.hat,nrow=5,ncol=1) #nrow= number of covariates + Intercept
meanxpop <- cbind(1,X.m_o)
eblup_out<- as.matrix(meanxpop)%*%betaest    
RESULTS_out<-data.frame(Domain=X.m_tmp$COD.PROV.SILC,
                        Eblup=eblup_out)

# variance can be obtain with smoothing method 
# (see Wolter, Introduction to Variance Estimation)

# calculate the CV %
CV.eblup<-100*sqrt(unit.eblup$mse$mse)/abs(unit.eblup$est$eblup$eblup)
CV.eblup1_t<- recode(CV.eblup, "0:16.5= '< 16.5%'; 16.51:33.3='between 16.6% and 33.3% '; else='>33.3%'")
table(CV.eblup1_t)


# dataframe with results for sampled areas
RESULTS<-data.frame(Domain=HT$Domain, n=HT$SampSize, 
                    Dir=HT$Direct, SDdir=HT$SD, CVdir=HT$CV,
                    Eblup=unit.eblup$est$eblup$eblup, RMSEeblup=sqrt(unit.eblup$mse$mse),
                    CVeblup=CV.eblup)

# merge results and shape file
ita_sae<- merge(ita, RESULTS, by.x="COD_PROV", by.y="Domain")

spplot(ita.merge, "Direct", main = "Direct Estimate NUTS 3")
# number of calsses
nclr <- 4
class <- classIntervals(ita_sae$Eblup, nclr, style="quantile")
range(RESULTS$Dir)
range(RESULTS$Eblup)
quantile(RESULTS$Dir)
class<- classIntervals(RESULTS$Eblup, n=nclr, style="fixed",
                       fixedBreaks=c(min(RESULTS$Dir), 14808.460,
                                     17708.132, 19483.948,
                                     max(RESULTS$Dir))) 
#create the same scale with the previous plot
#breakpoints
categ = class[["brks"]]

plotclr <- brewer.pal(nclr,"YlOrRd")

ita_sae$cutBHF <- cut(ita_sae$Eblup, breaks = categ, 
                      right = F, include.lowest=T) 

ita_sae$cutDir <- cut(ita_sae$Dir, breaks = categ, 
                      right = F, include.lowest=T) 

spplot(ita_sae, "cutBHF", col.regions=plotclr, 
       main = "Eblup Estimates (NUTS 3)",
       colorkey = list(space="bottom",height = 0.9, 
                       labels = list(at = seq(0.5, length(categ) -0.5), 
                                     labels = round(categ,2))))

spplot(ita_sae, "cutDir", col.regions=plotclr, 
       main = "Direct Estimates (NUTS 3)",
       colorkey = list(space="bottom",height = 0.9, 
                       labels = list(at = seq(0.5, length(categ) -0.5), 
                                     labels = round(categ,2))))

#plot all sae estimates 
tmp<-RESULTS[,c(1,6)]
head(tmp)


RESULTS_ALL<- rbind(RESULTS_out,tmp)

# merge results and shape file
ita_sae2<- merge(ita, RESULTS_ALL, by.x="COD_PROV", by.y="Domain")

ita_sae2$cutBHF <- cut(ita_sae2$Eblup, breaks = categ, 
                       right = F, include.lowest=T) 
spplot(ita_sae2, "cutBHF", col.regions=plotclr, 
       main = "Eblup Estimates (NUTS 3)",
       colorkey = list(space="bottom",height = 0.9, 
                       labels = list(at = seq(0.5, length(categ) -0.5), 
                                     labels = round(categ,2))))


# compare estimates
cor(RESULTS$Dir, RESULTS$Eblup) #the higher -> move into the same direction -> the better 

#bias diagnostic plot > show whether estimates good or not
plot(RESULTS$Dir,RESULTS$Eblup, col="black",lwd=1,lty=1,pch=1, 
     main="Bias diagnostic plot", 
     xlab="EBLUP estimate", ylab="Direct survey estimate", ylim = c(9000, 45000))
abline(lm(RESULTS$Dir~RESULTS$Eblup), col="blue",lwd=1,lty=1) 
abline(0,1, col="red",lwd=1,lty=2) #line with 100% corelation 
legend("topleft",legend= c("Fitted line", "45 degree line"),
       lty=c(1,2),col=c("blue","red"),lwd=c(1,2), cex=0.75)

# if the estimates follow the 45? line there are no differences in the estimates. 
# Differently from this plot you can see over-under estimates


#CVs of direct estimates and Eblup estimates
RESULTS<- RESULTS[order(RESULTS$n, na.last = FALSE, decreasing = FALSE), ]
plot(RESULTS$n, RESULTS$CVdir, type = "p",col = "red", 
     lwd = 1, pch = 16, lty = 1, ylab = "CV", 
     xlab = "sample size", cex.axis = 1,
     cex.lab = 1, ylim=c(1,30))
points(RESULTS$n,RESULTS$CVeblup, type = "p", 
       col = "green", lwd = 1, pch = 25, lty = 1)


# CI
#rank districts
RESULTS <- RESULTS[order(RESULTS$Eblup), ] #
RESULTS <- cbind(RESULTS, c(1:dim(RESULTS)[1]))
colnames(RESULTS)[9] <- "rankEBLUP"

library(plotrix)
plotCI(RESULTS$rankEBLUP, RESULTS$Eblup, 
       ui=RESULTS$Eblup+ 1.96*RESULTS$RMSEeblup, 
       li=RESULTS$Eblup- 1.96*RESULTS$RMSEeblup,
       xlab= "Districts (sorted by mean of income)",
       ylab = "Mean of income", 
       main="Confidence Interval Plot", ylim=c(5000,60000),xaxt="n")

plotCI(RESULTS$rankEBLUP+0.5, RESULTS$Dir, 
       ui=RESULTS$Dir+ 1.96*RESULTS$SDdir, 
       li=RESULTS$Dir- 1.96*RESULTS$SDdir, col="red" ,add=T)


text(1:70,RESULTS$Eblup + 2.1*RESULTS$RMSEeblup,RESULTS$Domain, 
     cex=0.6,srt = 90, adj=c(-0.1,0.3))

#CI overlap => cannot say statistical differences


# FH model
# 

# FH model the mean of eqhhincome
X.m1 <- subset(N.cens2, COD.PROV.SILC %in% samp$PROV) 

FH.data<-data.frame(Domain=HT$Domain, y.dir=HT$Direct, SD=HT$SD,X.m1[,-1])

area.eblup1<-mseFH(formula=FH.data$y.dir ~ FH.data$AGE_class1+ FH.data$MALE+ 
                     FH.data$THREE.MEMBERS+ FH.data$PROPERTY.HOUSE, 
                   vardir=FH.data$SD^2)
CV.areaeblup1<-100*sqrt(area.eblup1$mse)/abs(area.eblup1$est$eblup)
summary(CV.areaeblup1)
cor(FH.data$y.dir, area.eblup1$est$eblup)

RESULTS_FH<-data.frame(Domain=HT$Domain, n=HT$SampSize, 
                       Dir=HT$Direct, SDdir=HT$SD, CVdir=HT$CV,
                       Eblup=area.eblup1$est$eblup, 
                       RMSEeblup=sqrt(area.eblup1$mse),
                       CVeblup=CV.areaeblup1)


#bias diagnostic plot
plot(RESULTS_FH$Dir,RESULTS_FH$Eblup, col="black",lwd=1,lty=1,pch=1, 
     main="Bias diagnostic plot", 
     xlab="FH estimate", ylab="Direct survey estimate", ylim = c(9000, 45000))
abline(lm(RESULTS_FH$Dir~RESULTS_FH$Eblup), col="blue",lwd=1,lty=1) 
abline(0,1, col="red",lwd=1,lty=2) 
legend("topleft",legend= c("Fitted line", "45 degree line"),
       lty=c(1,2),col=c("blue","red"),lwd=c(1,2), cex=0.75)


#CVs of direct estimates and Eblup estimates
RESULTS_FH<- RESULTS_FH[order(RESULTS_FH$n, na.last = FALSE, decreasing = FALSE), ]
plot(RESULTS_FH$n, RESULTS_FH$CVdir, type = "p",col = "red", 
     lwd = 1, pch = 16, lty = 1, ylab = "CV", 
     xlab = "sample size", cex.axis = 1,
     cex.lab = 1, ylim=c(1,30))
points(RESULTS_FH$n,RESULTS_FH$CVeblup, type = "p", 
       col = "green", lwd = 1, pch = 25, lty = 1)
