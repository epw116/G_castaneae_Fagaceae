setwd("~/R files")
### Nut Rot Analysis

## Oak

oak <- read.csv("OakNutRotBothTests.csv")
oak$Rank <- factor(oak$Rank,order=TRUE,levels=c(0,1,2,3,4))
oak$Rank

# With only two levels of test, can't use random effect
# So treating as fixed, as a block

fm2 <- clm(Rank~Species*Treatment+Test,data=oak,HESS=TRUE)
anova(fm2,type="III")

pairs(emmeans(fm2,~Treatment))
pairs(emmeans(fm2,~Species))
pairs(emmeans(fm2,~Treatment | Species))

pairs(emmeans(fm1,~Treatment | Species))
pairs(emmeans(fm2,~Treatment | Species))

## Chesnut

cnut <- read.csv("ChestnutNutRotBothTests.csv")
cnut$Rank <- factor(cnut$Rank,order=TRUE,levels=c(0,1,2,3,4))
cnut$Rank

# With only two levels of test, can't use random effect
# So treating as fixed, as a block

fmc2 <- clm(Rank~Species*Treatment+Test,data=cnut,HESS=TRUE)
anova(fm2c,type="III")

pairs(emmeans(fmc2,~Treatment))
pairs(emmeans(fmc2,~Species))
pairs(emmeans(fmc2,~Treatment | Species))

pairs(emmeans(fmc1,~Treatment | Species))
pairs(emmeans(fmc2,~Treatment | Species))

### Canker

library(lme4)
library(car)
library(tidyverse)

canker <- read.csv("BranchCankerAllTests.csv")

## November

nov2 <- canker %>% filter(Test=="Nov") %>%
  mutate(Trt2 = case_when(grepl("Control",Treatment)~"Control",
                          TRUE~"Treatment")) %>%
  separate(Treatment,c("Location","Other"))

fm2bc <- lm(bcPower(Canker.Area,-2/3)~Trt2*Species+Location,data=nov2)
Anova(fm2bc,type="III")

pairs(emmeans(fm2bc,~Trt2))
pairs(emmeans(fm2bc,~Species))
pairs(emmeans(fm2bc,~Trt2 | Species))

bctranendo1 <- make.tran("boxcox",-2/3)
res.tbendo <- with(bctranendo1,
                   lm(linkfun(Canker.Area) ~Trt2*Species+Location,data=nov2,
                      contrasts=list(Trt2=contr.sum, Species=contr.sum)))
int.emmendo1<-emmeans(res.tbendo,~Trt2:Species)
summary(int.emmendo1,type="response")

## February

feb2 <- canker %>% filter(Test=="Feb") %>%
  mutate(Trt2 = case_when(grepl("Control",Treatment)~"Control",
                          TRUE~"Treatment")) %>%
  separate(Treatment,c("Location","Other"))

fm2bc <- lm(bcPower(Canker.Area,-2/3)~Trt2*Species+Location,data=feb2)
Anova(fm2bc,type="III")

pairs(emmeans(fm2bc,~Trt2))
pairs(emmeans(fm2bc,~Species))
pairs(emmeans(fm2bc,~Trt2 | Species))

bctraneco1 <- make.tran("boxcox",-2/3)
res.tbeco <- with(bctraneco1,
                   lm(linkfun(Canker.Area) ~Trt2*Species+Location,data=feb2,
                      contrasts=list(Trt2=contr.sum, Species=contr.sum)))
int.emmeco1<-emmeans(res.tbeco,~Trt2:Species)
summary(int.emmeco1,type="response")


## April

apr2 <- canker %>% filter(Test=="Apr") %>%
  mutate(Trt2 = case_when(grepl("Control",Treatment)~"Control",
                          TRUE~"Treatment")) %>%
  separate(Treatment,c("Location","Other"))

fm2bc <- lm(bcPower(Canker.Area,-1/10)~Trt2*Species+Location,data=apr2)
Anova(fm2bc,type="III")

pairs(emmeans(fm2bc,~Trt2))
pairs(emmeans(fm2bc,~Species))
pairs(emmeans(fm2bc,~Trt2 | Species))

bctranact1 <- make.tran("boxcox",-1/10)
res.tbact <- with(bctranact1,
                  lm(linkfun(Canker.Area) ~Trt2*Species+Location,data=apr2,
                     contrasts=list(Trt2=contr.sum, Species=contr.sum)))
int.emmact1<-emmeans(res.tbact,~Trt2:Species)
summary(int.emmact1,type="response")

### Stem

stem1 <- read.csv("StemCankerTest.csv")

library(lme4)
library(car)

## Stem Canker 

fm1bc <- lm(bcPower(Canker.Area,-.5)~Treatment*Species,data=stem1)
Anova(fm1bc,type="III")

pairs(emmeans(fm1bc,~Treatment))
pairs(emmeans(fm1bc,~Species))
pairs(emmeans(fm1bc,~Treatment | Species))

bctranstem1 <- make.tran("boxcox",-.5)
res.tbstem <- with(bctranstem1,
                  lm(linkfun(Canker.Area) ~Treatment*Species,data=stem1,
                     contrasts=list(Treatment=contr.sum, Species=contr.sum)))
int.emmstem1<-emmeans(res.tbstem,~Treatment:Species)
summary(int.emmstem1,type="response")

## Stem Girdling

stem2 <- read.csv("StemTestGirdle.csv",stringsAsFactors = TRUE)
stem2$Rank <- factor(stem2$Rank,order=TRUE,levels=c(0,1,2,3,4))

library(rstanarm)

sglma <- stan_polr(Rank~Species*Treatment,data=stem2,
                   prior = R2(0.25), prior_counts = dirichlet(1))
summary(sglma)

sims <- as.matrix(sglma)

newdata <- data.frame(Species=c(stem2$Species[1],stem2$Species[1]),
                      Treatment=c(stem2$Treatment[1],stem2$Treatment[2]))
newdata
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.00625,0.99375))
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.025,0.975))

newdata <- data.frame(Species=c(stem2$Species[3],stem2$Species[3]),
                      Treatment=c(stem2$Treatment[1],stem2$Treatment[2]))
newdata
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.00625,0.99375))
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.025,0.975))

newdata <- data.frame(Species=c(stem2$Species[5],stem2$Species[5]),
                      Treatment=c(stem2$Treatment[1],stem2$Treatment[2]))
newdata
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.00625,0.99375))
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.025,0.975))

newdata <- data.frame(Species=c(stem2$Species[6],stem2$Species[6]),
                      Treatment=c(stem2$Treatment[1],stem2$Treatment[2]))
newdata
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.00625,0.99375))
quantile(posterior_linpred(sglma,newdata=newdata)[,1]-
           posterior_linpred(sglma,newdata=newdata)[,2],
         probs = c(0.025,0.975))