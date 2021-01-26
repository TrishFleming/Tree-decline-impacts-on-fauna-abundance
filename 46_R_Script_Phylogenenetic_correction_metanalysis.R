# Fleming et al 2021 Global meta-analysis of tree decline impacts on fauna. Biological Reviews

# SETUP LIBRARIES AND WORKING DIRECTORY -----
R.Version()
library(lattice)
library(nlme)
library(mgcv)
library(lme4)
library(MASS)
library(ggplot2)
library(plyr)
library(Rmisc)
library(psych) # the command describe carries out summary statistics
library(lm.beta)
library(cowplot) #to arrange graphs on a grid
library(lmerTest)
library(MuMIn) 
library(DHARMa)
library(ggeffects) #for prediction
library(forcats) #order by effect size
library(metafor)
library(ape) #(Paradis et al 2004)
library(phytools)
library(MCMCglmm) # for meta-analysis
library(dplyr) # for data manipulation
library(MCMCglmm) #(Hadfield, 2010a)    https://ourcodingclub.github.io/tutorials/mcmcglmm/ 
library(buildr) #to build formulae
library(robumeta)  #https://www.youtube.com/watch?v=d1pYHfCKhyA
library(metafor)   #https://www.youtube.com/watch?v=d1pYHfCKhyA
library(dplyr)     #https://www.youtube.com/watch?v=d1pYHfCKhyA

# (1) DATA processing -----
library(dmetar)  # https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/ 
library(esc)
library(metafor)  # tutorial https://stats.stackexchange.com/questions/116659/mixed-effects-meta-regression-with-nested-random-effects-in-metafor-vs-mixed-mod 

#__(1a) Calculate effect sizes and their variance from raw data using Metafor  -----
setwd("C:/Users/20043940/Dropbox/TF data/1. Tree declines")

#Hedges g from means and sd https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/f.html 
yellow<- read.csv("43 calc g.csv", header = TRUE) 
g_calcMEANS <- escalc(n1i = n1, n2i = n2, m1i = Decline.mean1, m2i = Intact.mean2, 
                      sd1i = sd1, sd2i = sd2, data = yellow, measure = "SMD",  append = TRUE)
write.csv(g_calcMEANS, "means.csv")

yellowR<- read.csv("43 calc g R.csv", header = TRUE) 
#Hedges g from r https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/f.html 
g_calcR <- escalc(measure="COR",ri=r, ni=n.captured, data=yellowR,  append = TRUE)
write.csv(g_calcR, "corrR.csv")

yellowT<- read.csv("T_Drever.csv", header = TRUE) 
g_calcT <- escalc(measure="COR",ti=t, ri=r, ni=n.captured, data=yellowT,  append = TRUE)
write.csv(g_calcT, "t_tests.csv")

yellowX<- read.csv("Adams chi.csv", header = TRUE) 
g_calcX <- escalc(measure = "SMD",  chisq=chi, n1i = n1, n2i = n2, totaln=n.captured, data=yellowX,  append = TRUE)
esc_chisq(chisq=7.488095238,totaln=14.97619048,es.type="cox.or")

yellowX<- read.csv("calc g from Chi.csv", header = TRUE) 
g_calcX <- escalc(measure="SMD",chisq=chi, ni = n.captured, data=yellowX,  append = TRUE)   
esc_chisq(chisq=7.488095238,totaln=14.97619048,es.type="cox.or")
esc_chisq(chisq=9.9,totaln=100,es.type="cox.or")
escal

write.csv(g_calcT, "t_tests.csv")
yellowCounts <- escalc(measure="COR",ti=t, ri=r, ni=n.captured, data=yellowT,  append = TRUE)
write.csv(yellowCounts, "AlsopCounts.csv")

#__(1b) Test normality -----
setwd("C:/Users/20043940/Dropbox/TF data/1. Tree declines/Phylogenetic corrections/ALL birds")
black<- read.csv("allDATA.csv")
#or plot by group, e.g. for birds:
blackBIRDS <- black[which(black$Class=='Bird'),]

qqnorm(black$yi)       #this does not follow the line        
qqline(black$yi)
qqnorm(black$log_yi)   #this follows the line
qqline(black$log_yi) 


#__(1c) structure of random moderators -----
#__based on Metafor tutorials: https://ourcodingclub.github.io/tutorials/mcmcglmm/ and https://rstudio-pubs-static.s3.amazonaws.com/10913_5858762ec84b458d89b0f4a4e6dd5e81.html
plot(black$yi, (1/black$vi)) # draw funnel plot of slope and precision (1/SE)
#NOTE the data derived from chi tests had very small vi, produced large values for the inverse function
randomtest <- MCMCglmm(log_yi ~ 1, random = ~Publication + SpeciesSynonyms, data = black)
summary(randomtest)
par(mfrow = c(1,2))
hist(mcmc(randomtest$VCV)[,"Publication"])
hist(mcmc(randomtest$VCV)[,"SpeciesSynonyms"])
par(mfrow=c(1,1)) # Reset the plot panel back to single plots
plot(randomtest$Sol) #tests for model convergence.  The trace can be used to assess mixing (or convergence), while the density is like a smoothed histogram of the estimates of the posterior distribution that the model produced for every iteration of the model. To make sure your model has converged, the trace plot should look like a fuzzy caterpillar.  
plot(randomtest$VCV) #tests relationship for random factors.  
a <- 1000
prior1 <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
randomerror <- MCMCglmm(Log.effect.size.g ~ 1, random = ~Publication + SpeciesSynonyms, 
                        data = black, prior = prior1, nitt = 50000)
summary(randomerror)
plot(randomerror$VCV) #there is now a larger effective sample size.  The models look to have mixed much better too. This is also good.

#__(1d) Heterogeneity analysis -----
# Based on https://www.youtube.com/watch?v=d1pYHfCKhyA How to perform a transparent meta-analysis; Daniel Quintan. Written by Daniel S. Quintana.  A non-technical primer for conducting a meta-analysis to synthesize correlational data  If you have any questions or comments, contact me over via email daniel.quintana@medisin.uio.no or Twitter @dsquintana. Any updates to the script will be posted to http://github.com/dsquintana/corr_meta

setwd("C:/Users/20043940/Dropbox/TF data/1. Tree declines")
black<- read.csv("allDATA.csv")
#NOTE: data are sorted from largest negative to largest positive effects
blackBIRDS <- black[which(black$Class=='Bird'),]
blackMAMMS <- black[which(black$Class=='Mammal'),]
blackREPS<- black[which(black$Class=='Reptile'),]
blackINVERTS<- black[which(black$Taxon=='Arthropod'),]

resBIRDS <- rma(yi, vi, data=blackBIRDS) 
resMAMMS <- rma(yi, vi, data=blackMAMMS) 
resREPS <- rma(yi, vi, data=blackREPS) 
resINVERTS <- rma(yi, vi, data=blackINVERTS)
resBIRDS
resMAMMS
resREPS
resINVERTS

resBIRDS_log <- rma(log_yi, vi, data=blackBIRDS) 
resMAMMS_log <- rma(log_yi, vi, data=blackMAMMS) 
resREPS_log <- rma(log_yi, vi, data=blackREPS) 
resINVERTS_log <- rma(log_yi, vi, data=blackINVERTS)
resBIRDS_log
resMAMMS_log 
resREPS_log
resINVERTS_log

# example outputs for 'resBIRDS'
# "Random-Effects Model (k = 446; tau^2 estimator: REML)" <-  used a random-effects model with k=446 cases and that the degree of heterogeneity (tau^2) was calculated using a restricted maximum-likelihood estimator.
# "tau^2 (estimated amount of total heterogeneity): 0.9584 (SE = 0.0792)"
# "I^2 (total heterogeneity / total variability):   100.00%"; i.e. 100% of variation reflected actual differences in the population mean.
# "Test for Heterogeneity: Q(df = 445) = 47084097974.6844, p-val < .0001" <- significant level of variation
# Significant p-value for model results indicates that the included studies do NOT share a common effect size.
# "estimate" = estimated model coefficient (i.e., the summary effect size) 

# Bajaut plot can inidicate which studies influencing overall heterogeneity. 
# Studies that fall to the top right quadrant of the Baujat plot contribute most to both these factors. 
b_resBIRDS <- rma(yi, vi, data=blackBIRDS, slab=RecordID)    
baujat(b_resBIRDS)
b_resMAMMS <- rma(yi, vi, data=blackMAMMS, slab=RecordID)    
baujat(b_resMAMMS)
b_resREPS <- rma(yi, vi, data=blackREPS, slab=RecordID)   
baujat(b_resREPS)
b_resINVERTS <- rma(yi, vi, data=blackINVERTS, slab=RecordID)    
baujat(b_resINVERTS)

# __(1e) Diagnostics to identify potential outliers and influential cases -----
# Datapoints marked with an asterisk fulfil the criteria as an 'influential study'. 
inf <- influence(resBIRDS)
inf_log <- influence(resBIRDS_log)
inf <- influence(resMAMMS)
inf_log <- influence(resMAMMS_log)
inf <- influence(resREPS)
inf_log <- influence(resREPS_log)
inf <- influence(resINVERTS)
inf_log <- influence(resINVERTS_log)
print(inf)
plot(inf)
print(inf_log)
plot(inf_log)

# Forest plots with effect size and 95% CIs.  Studies with larger squares contributed more to the summary effect size. Summary effect size represented as the polygon at the bottom
fpBIRDS <- forest(resBIRDS, xlim=c(-12,12), digits=c(1,1), cex=.8)
# too large for plotting
# print to file ...
fpMAMMS <- forest(resMAMMS, xlim=c(-10,10), digits=c(2,1), cex=.8)
pdf("fig fpMAMMS.jpg", width = 6, height = 10) # Open a new pdf file
plot_grid(fpMAMMS, align = "v", nrow = 1, rel_heights = c(1))
dev.off() # Close the file
# ... or print on screen
fpREPS <- forest(resREPS, xlim=c(-8,10), digits=c(2,1), cex=.8)
fpREPS
fpINVERTS <- forest(resINVERTS, xlim=c(-10,10), digits=c(2,1), cex=.8)
fpINVERTS

### funnel plot
funnel(resBIRDS, xlab = "Effect size (g)")
funnel(resMAMMS, xlab = "Effect size (g)")
funnel(resREPS, xlab = "Effect size (g)")
funnel(resINVERTS, xlab = "Effect size (g)")

# __(1f) Diagnostics to identify potential publication bias -----
# Egger's regression and Rank correlation tests for evidence of publication bias (asymmetry in effect sizes).
regtest(resBIRDS)
regtest(resMAMMS)
regtest(resREPS)
regtest(resINVERTS)
ranktest(resBIRDS)
ranktest(resMAMMS)
ranktest(resREPS)
ranktest(resINVERTS)


# (2) MCMCglmm including phylogenetic correction -----
# Based on 
#   Fayard et al https://onlinelibrary.wiley.com/doi/abs/10.1111/brv.12606 
#   http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm
# models need the species listed as 'animal', and a column called 'phylo' <- 'name of the matrix should represent the column in data associated with random term'
# Make sure there are no spaces between genus and speceies - use an underscore

#    split.direct.sum <- deals with error due to splitting long formulae -----
split.direct.sum<-function(x){  
  
  if(is.na(x)){
    return(NULL)
  }else{
    x <- gsub("\n", "", x, fixed=TRUE)  # required if the random formula is very long and broken over lines
    openB<-gregexpr("\\(", x)[[1]]
    closeB<-gregexpr("\\)", x)[[1]]
    true_openB<-openB
    true_closeB<-closeB
    
    for(i in 1:length(openB)){
      dist<-outer(openB, closeB, function(x,y){y-x})
      dist[which(dist<0)]<-Inf
      dist<-which(dist==min(dist), arr.ind=T)[1,] 
      true_openB[i]<-openB[dist[1]]
      true_closeB[i]<-closeB[dist[2]]
      openB<-openB[-dist[1]]
      closeB<-closeB[-dist[2]]
    }
    
    plus<-gregexpr("\\+", x)[[1]]
    internals<-matrix(mapply(function(x,y){(plus>x & plus<y)}, x=true_openB, y=true_closeB), length(plus), length(true_openB))
    rterms<-strsplit(x, "")[[1]]
    rterms[plus[which(rowSums(internals)!=0)]]<-"leaveMCMCleave"
    rterms<-paste(rterms, collapse="")
    rterms<-strsplit(rterms, " *\\+ *")[[1]]
    rterms<-gsub("leaveMCMCleave", "+", rterms)
    return(rterms)
  }
}





# (b) BIRDS including SD for ALL data n=446 effect size----
# notes (1) S=0, T=1, (2) V included in primary diet, (3) make sure that species names have an underscore, not space
# (4) ordered from ground to trees - ground is baseline
#__ (b1) Bird data -----
setwd("C:/Users/20043940/Dropbox/TF data/Tree declines/Phylogenetic corrections/ALL birds")
# DATA
black<- read.csv("allDATA.csv")
blackBIRDS <- black[which(black$Class=='Bird'),]
blackBIRDS$animal <- as.factor(blackBIRDS$animal) # Pedigree function needs column called "animal"
blackBIRDS$Nesting.guild <- factor(x=blackBIRDS$Nesting.guild, levels=c("GRO", "LOW", "CAV", "TREE")) #to set order - comparison with GRO
blackBIRDS$Foraging.guild <- factor(x=blackBIRDS$Foraging.guild, levels=c("GRO", "LOW", "BAR", "CAN", "AER")) #to set order - comparison 
blackBIRDS$Primary.diet <- factor(x=blackBIRDS$Primary.diet, levels=c("O", "V", "I", "S", "N")) #to set order - comparison excludes "V" = combined with O 
blackBIRDS$Decline.cause<- as.factor(blackBIRDS$Decline.cause)
blackBIRDS$Temporal.or.Spatial<- scale(blackBIRDS$Temporal.or.Spatial) # S=0, T=1
blackBIRDS$Decline.permanence<- scale(blackBIRDS$Decline.permanence) #use as continuous, not categorical
blackBIRDS$Decline.scale<- scale(blackBIRDS$Decline.scale) #use as continuous, not categorical
blackBIRDS$Decline.speed<- scale(blackBIRDS$Decline.speed)
# PHYLOGENY
phylo<-read.nexus("phylo.nex")
summary.phylo(phylo) # checks if tree has imported correctly (number of tips)
phylo$tip.label
phylo$node.label <- NULL
#phylo2=drop.tip(phylo, tip=c(2:3)) # TRIMMING OPTION to keep 1 rotifera (to root) and other acantho
is.rooted(phylo) # checks if the tree is rooted (necessary)
phylo=force.ultrametric(phylo, method="extend") # to make the tree ultrametric
is.ultrametric(phylo) # checks if the tree is ultrametric (necessary)    TRUE
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)
plot(phylo)  #too large to plot

#__ (b2) Sensitivity analysis on global models testing effect of (1) raw or log-transformed g, (2) the inclusion of species as a random factor, (3) the inclusion of phylogenetic correction (PC) as a random factor-----
priorSD=list(R=list(R1=list(V=1,nu=0.002),R2=list(V=1,nu=0.002),R3=list(V=1,nu=0.002),
                    R4=list(V=1,nu=0.002),R5=list(V=1,nu=0.002),R6=list(V=1,nu=0.002)), 
             G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))  
prior_noSD=list(R=list(R1=list(V=1,nu=0.002),R2=list(V=1,nu=0.002),R3=list(V=1,nu=0.002),
                       R4=list(V=1,nu=0.002),R5=list(V=1,nu=0.002),R6=list(V=1,nu=0.002)), 
                G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))


Mlmer_blackBIRD <-lmer(yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed +  Primary.diet * Foraging.guild  +
                   Nesting.guild + (1|Publication), data=blackBIRDS)
Mglmm_blackBIRD <- MCMCglmm(yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed +  Primary.diet * Foraging.guild  +
                        Nesting.guild,  random=~Publication + animal, singular.ok=TRUE,prior=NULL, nitt=50000,
                      thin=250,burnin=1000, data=blackBIRDS)# Bird species as a random factor, but no phylogeny
Mpc_blackBIRD<-MCMCglmm(yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed +  
                          Primary.diet * Foraging.guild  + Nesting.guild,random=~Publication + animal, 
                  family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior,
                  data=blackBIRDS,nitt=50000,burnin=1000,thin=500, singular.ok=TRUE)   
MlmerLOGblackBIRDS_noSD <-lmer(log_yi ~ Temporal.or.Spatial + Decline.speed + Decline.permanence	+ Decline.scale	+ 
                              Primary.diet * Foraging.guild  + Nesting.guild + (1|Publication), data=blackBIRDS)  
MglmmLOGblackBIRDS_noSD <- MCMCglmm(log_yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+  Decline.speed + 
                                   Primary.diet * Foraging.guild  + Nesting.guild, random=~Publication + Species, 
                                 singular.ok=TRUE,  prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackBIRDS) # Bird species as a random factor, but no phylogeny
MpcLOGblackBIRDS_noSD <-MCMCglmm(log_yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+  Decline.speed + 
                                   Primary.diet * Foraging.guild  +  Nesting.guild, 
                                 random=~Publication + Species, family="gaussian",
                              ginverse=list(phylo=inv.phylo$Ainv),
                              prior=NULL, data=blackBIRDS,nitt=50000,burnin=1000,thin=500)
Mglmm_yi_blackBIRDS_SD <- MCMCglmm(yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed + 
                                    Primary.diet * Foraging.guild  + Nesting.guild, random=~Publication + animal + idh(vi):units, 
                                  singular.ok=TRUE,  prior=priorSD, nitt=50000,thin=250,burnin=1000, data=blackBIRDS) # Bird species as a random factor, but no phylogeny
Mpc_yi_blackBIRD_SD <-MCMCglmm(yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed +
                                 Primary.diet * Foraging.guild  +  Nesting.guild, 
                              random=~Publication + animal + idh(vi):units, family="gaussian",
                              ginverse=list(phylo=inv.phylo$Ainv),
                              prior=priorSD, data=blackBIRDS,nitt=50000,burnin=1000,thin=500)    # top model with phylogeny and SD

MglmmLOGblackBIRDS_SD <- MCMCglmm(log_yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed + 
                                    Primary.diet * Foraging.guild  + Nesting.guild, 
                                  random=~Publication + animal + idh(vi):units, 
                                  singular.ok=TRUE,  prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackBIRDS) # Bird species as a random factor, but no phylogeny
# top model:
MpcLOGblackBIRD_SD <-MCMCglmm(log_yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed + 
                                Primary.diet * Foraging.guild  +  Nesting.guild, 
                              random=~Publication + animal + idh(vi):units, family="gaussian",
                              ginverse=list(phylo=inv.phylo$Ainv),
                              prior=priorSD, data=blackBIRDS,nitt=50000,burnin=1000,thin=500)    # top model with phylogeny and SD
BIC(Mlmer_blackBIRD , Mglmm_blackBIRD , Mglmm_yi_blackBIRDS_SD , Mpc_blackBIRD, Mpc_yi_blackBIRD_SD ,MlmerLOGblackBIRDS_noSD , MglmmLOGblackBIRDS_noSD , MglmmLOGblackBIRDS_SD , MpcLOGblackBIRDS_noSD , MpcLOGblackBIRD_SD )
DIC(Mlmer_blackBIRD , Mglmm_blackBIRD , Mglmm_yi_blackBIRDS_SD , Mpc_blackBIRD, Mpc_yi_blackBIRD_SD ,MlmerLOGblackBIRDS_noSD , MglmmLOGblackBIRDS_noSD , MglmmLOGblackBIRDS_SD , MpcLOGblackBIRDS_noSD , MpcLOGblackBIRD_SD )
AIC(Mlmer_blackBIRD , Mglmm_blackBIRD , Mglmm_yi_blackBIRDS_SD , Mpc_blackBIRD, Mpc_yi_blackBIRD_SD ,MlmerLOGblackBIRDS_noSD , MglmmLOGblackBIRDS_noSD , MglmmLOGblackBIRDS_SD , MpcLOGblackBIRDS_noSD , MpcLOGblackBIRD_SD )

#__ (b2) Dredge on pc-corrected model -----
gm_blackBIRDS <-MCMCglmm.updateable(log_yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed +
                                      Primary.diet * Foraging.guild  +  Nesting.guild, 
                                       random=~Publication + animal + idh(vi):units, family="gaussian",
                                       ginverse=list(phylo=inv.phylo$Ainv),
                                       prior=priorSD, data=blackBIRDS,nitt=50000,burnin=1000,thin=500)
options(na.action = "na.fail"); dredge.gm_blackBIRDS<- dredge(gm_blackBIRDS, rank="DIC", m.lim = c(0, 5),evaluate=T); options(na.action = "na.omit")
dredge.gm_blackBIRDS # this final version has vertebrate-eaters grouped with omnivores
write.csv(dredge.gm_blackBIRDS, "dredge_blackBIRDS.csv")

retained_blackBIRDS <-MCMCglmm(log_yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed + 
                                      Primary.diet  +  Nesting.guild, 
                                    random=~Publication + animal + idh(vi):units, family="gaussian",
                                    ginverse=list(phylo=inv.phylo$Ainv),
                                    prior=priorSD, data=blackBIRDS,nitt=50000,burnin=1000,thin=500)
summary(retained_blackBIRDS)

top.models <- get.models(dredge.gm_blackBIRDS, subset = delta<2)
model.avg <- summary(model.avg(top.models))
MA.est.table <- model.avg$coefmat.full
CI <- as.data.frame(confint(model.avg)) ## calculate confidence intervals 
df <- as.data.frame(MA.est.table)
wCI <-cbind(df, CI)
col1 <- rep("Birds", length(df$Estimate))
df2 <- cbind(col1, wCI)
sw <- as.data.frame(model.avg$sw)
col1 <- rep("Birds", length(model.avg$sw))
sw2 <- cbind(col1, sw); names(sw2) <- c("Type", "SW")
summary(model.avg)

#__ (b4) calculate error and distribution of data -----
# funnel plots



# (m) MAMMALS -----
# __(m1) Mammal DATA -----
setwd("C:/Users/20043940/Dropbox/TF data/Tree declines/Phylogenetic corrections/Mammals with SD")
#NOTES: (1) S=0, T=1, (2) make sure that species names have an underscore, not space
# (3) ordered from ground to trees - ground is baseline, (4) Smithopsis/Antechinus renamed Smithopsis
# DATA
black<- read.csv("allDATA.csv")
blackMAMMS <- black[which(black$Class=='Mammal'),]
blackMAMMS$Nesting.guild <- factor(x=blackMAMMS$Nesting.guild, levels=c("GRO", "TREE")) #to set order - comparison with GRO
blackMAMMS$Foraging.guild <- factor(x=blackMAMMS$Foraging.guild, levels=c("GRO", "CAN", "AER")) #to set order - comparison 
blackMAMMS$Primary.diet <- factor(x=blackMAMMS$Primary.diet, levels=c("I", "L", "N", "S")) #to set order - comparison 
blackMAMMS$animal <- as.factor(blackMAMMS$animal) # Pedigree function needs column called "animal"
blackMAMMS$Log.Mb <- scale(blackMAMMS$Log.Mb)
blackMAMMS$Decline.cause<- as.factor(blackMAMMS$Decline.cause)
blackMAMMS$Temporal.or.Spatial<- scale(blackMAMMS$Temporal.or.Spatial) # S=0, T=1
blackMAMMS$Decline.permanence<- scale(blackMAMMS$Decline.permanence) #use as continuous, not categorical
blackMAMMS$Decline.scale<- scale(blackMAMMS$Decline.scale) #use as continuous, not categorical
blackMAMMS$Decline.speed<- scale(blackMAMMS$Decline.speed) #use as continuous, not categorical
# PHYLOGENY
phylo<-read.nexus("phylo.nex")# MAMMALS   WARNING -> name of the matrice should represent the column in data associated with random term
phylo=read.tree("phylo.csv")  # BIRDS this version works
summary.phylo(phylo) # checks if tree has imported correctly (number of tips)
phylo$tip.label
phylo$node.label <- NULL
#phylo2=drop.tip(phylo, tip=c(2:3)) # TRIMMING OPTION to keep 1 rotifera (to root) and other acantho
is.rooted(phylo) # checks if the tree is rooted (necessary)
phylo=force.ultrametric(phylo, method="extend") # to make the tree ultrametric
is.ultrametric(phylo) # checks if the tree is ultrametric (necessary)    TRUE
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)
plot(phylo)
# __(m2) Sensitivity analysis on global models testing effect of (1) raw or log-transformed g, (2) the inclusion of species as a random factor, (3) the inclusion of phylogenetic correction (PC) as a random factor -----
Mlmer_MAMMAL <-lmer(yi ~ Temporal.or.Spatial + method.of.capture + 
                      Decline.permanence	+ Decline.scale	+ Decline.speed + 
                      Primary.diet * Foraging.guild  + Nesting.guild + Log.Mb  +
                      (1|Publication), data=blackMAMMS)
Mglmm_MAMMAL <- MCMCglmm(yi ~ Temporal.or.Spatial + method.of.capture + 
                           Decline.permanence	+ Decline.scale	+ Decline.speed + 
                           Primary.diet * Foraging.guild  + Nesting.guild+ Log.Mb,
                         random=~Publication + animal, singular.ok=TRUE,
                         prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackMAMMS)       # Bird species as a random factor, but no phylogeny
Mpc_MAMMAL <-MCMCglmm(yi ~ Temporal.or.Spatial + method.of.capture + 
                        Decline.permanence	+ Decline.scale	+ Decline.speed +  
                        Primary.diet * Foraging.guild  + Nesting.guild+ Log.Mb,
                      random=~Publication + animal, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                      prior=prior,              data=blackMAMMS, nitt=50000,burnin=1000,thin=500)  #with phylogeny    
Mglmm_MAMMAL_SD <- MCMCglmm(yi ~ Temporal.or.Spatial + method.of.capture + 
                           Decline.permanence	+ Decline.scale	+ Decline.speed + 
                           Primary.diet * Foraging.guild  + Nesting.guild+ Log.Mb,
                         random=~Publication + animal + idh(vi):units, singular.ok=TRUE,
                         prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackMAMMS)       # Bird species as a random factor, but no phylogeny
Mpc_MAMMAL_SD <-MCMCglmm(yi ~ Temporal.or.Spatial + method.of.capture + 
                        Decline.permanence	+ Decline.scale	+ Decline.speed +  
                        Primary.diet * Foraging.guild  + Nesting.guild+ Log.Mb,
                      random=~Publication + animal + idh(vi):units, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                      prior=NULL,  data=blackMAMMS, nitt=50000,burnin=1000,thin=500)  #with phylogeny

# USING HEDGES G TRANSFORMED Log.effect.size.g
#Top model:
MlmerLOG_MAMMAL <-lmer(log_yi ~ Temporal.or.Spatial + method.of.capture + 
                         Decline.permanence	+ Decline.scale	+ Decline.speed +  
                         Primary.diet * Foraging.guild  + Nesting.guild + Log.Mb + 
                         (1|Publication), data=blackMAMMS)

MglmmLOG_MAMMAL <- MCMCglmm(log_yi~ Temporal.or.Spatial + method.of.capture + 
                              Decline.permanence	+ Decline.scale	+ Decline.speed +  
                              Primary.diet * Foraging.guild  + Nesting.guild + Log.Mb,
                            random=~Publication + animal, singular.ok=TRUE,
                            prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackMAMMS) # Bird species as a random factor, but no phylogeny

MpcLOG_MAMMAL<-MCMCglmm(log_yi ~ Temporal.or.Spatial + method.of.capture + 
                          Decline.permanence	+ Decline.scale	+ Decline.speed + 
                          Primary.diet * Foraging.guild  + Nesting.guild + Log.Mb,
                        random=~Publication + animal, 
                        family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                        prior=prior, data=blackMAMMS, nitt=50000,burnin=1000,thin=500)  #with phylogeny

MglmmLOG_MAMMAL_SD <- MCMCglmm(log_yi ~ Temporal.or.Spatial + method.of.capture + 
                              Decline.permanence	+ Decline.scale	+ Decline.speed +  
                              Primary.diet * Foraging.guild  + Nesting.guild + Log.Mb,
                            random=~Publication + animal + idh(vi):units, singular.ok=TRUE,
                            prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackMAMMS) # Bird species as a random factor, but no phylogeny

MpcLOG_MAMMAL_SD <-MCMCglmm(log_yi ~ Temporal.or.Spatial + method.of.capture + 
                          Decline.permanence	+ Decline.scale	+ Decline.speed + 
                          Primary.diet * Foraging.guild  + Nesting.guild + Log.Mb,
                        random=~Publication + animal + idh(vi):units, 
                        family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                        prior=NULL, data=blackMAMMS, nitt=50000,burnin=1000,thin=500)  #with phylogeny
#or try prior=priorSD  

AIC(Mlmer_MAMMAL,Mglmm_MAMMAL,Mglmm_MAMMAL_SD, Mpc_MAMMAL,Mpc_MAMMAL_SD,   MlmerLOG_MAMMAL,MglmmLOG_MAMMAL,MglmmLOG_MAMMAL_SD, MpcLOG_MAMMAL, MpcLOG_MAMMAL_SD)
BIC(Mlmer_MAMMAL,Mglmm_MAMMAL,Mglmm_MAMMAL_SD, Mpc_MAMMAL,Mpc_MAMMAL_SD,   MlmerLOG_MAMMAL,MglmmLOG_MAMMAL,MglmmLOG_MAMMAL_SD, MpcLOG_MAMMAL, MpcLOG_MAMMAL_SD)
DIC(Mlmer_MAMMAL,Mglmm_MAMMAL,Mglmm_MAMMAL_SD, Mpc_MAMMAL,Mpc_MAMMAL_SD,   MlmerLOG_MAMMAL,MglmmLOG_MAMMAL,MglmmLOG_MAMMAL_SD, MpcLOG_MAMMAL, MpcLOG_MAMMAL_SD)

# __(m3) Dredge on pc-corrected model -----
# global model (all fixed effects) 
prior=list(R=list(R1=list(V=1,nu=0.002),R2=list(V=1,nu=0.002),R3=list(V=1,nu=0.002),R4=list(V=1,nu=0.002),R5=list(V=1,nu=0.002),R6=list(V=1,nu=0.002)),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))

#gm = global model     species are listed in the database as animal	phylo	SpeciesSynonyms
MCMCglmm.updateable<- updateable(MCMCglmm) # updates dredge (MuMin package) function to do models selection with bayesian models 
gm_blackMAMMAL<- MCMCglmm.updateable(log_yi ~ Temporal.or.Spatial + method.of.capture + 
                                  Decline.permanence	+ Decline.scale	+ Decline.speed + 
                                  Primary.diet * Foraging.guild  + Nesting.guild + Log.Mb,
                                random=~Publication + animal + idh(vi):units, 
                                family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                                prior=NULL, data=blackMAMMS, nitt=50000,burnin=1000,thin=500)
summary(gm_blackMAMMAL)
options(na.action = "na.fail"); dredge.gm_blackMAMMAL<- dredge(gm_blackMAMMAL, rank="DIC",evaluate=T); options(na.action = "na.omit")  # DIC is the bayesian equivalent for AIC
dredge.gm_blackMAMMAL  # table with models ranked by values of DIC (best model=lowest DIC)
write.csv(dredge.gm_blackMAMMAL, "dredge_blackMAMMS.csv")

retained_blackMAMM<- MCMCglmm(log_yi ~ Temporal.or.Spatial + Decline.permanence	+ Decline.scale	+ Decline.speed + 
                                       Foraging.guild  + Nesting.guild,
                                     random=~Publication + animal + idh(vi):units, 
                                     family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                                     prior=NULL, data=blackMAMMS, nitt=50000,burnin=1000,thin=500)
summary(retained_blackMAMM)

best.model_MAMMAL <- get.models(dredge.gm_blackMAMMAL, subset = delta==0) 
top.models_MAMMAL <- get.models(dredge.gm_blackMAMMAL, subset = delta<2)
model.avg_MAMMAL <- summary(model.avg(top.models_MAMMAL))
MA.est.table_MAMMAL <- model.avg_MAMMAL$coefmat.full
CI <- as.data.frame(confint(model.avg_MAMMAL)) ## calc confidence intervals 
df <- as.data.frame(MA.est.table_MAMMAL)
wCI <-cbind(df, CI)
col1 <- rep("Mammals", length(df$Estimate))
df_MAMMAL <- cbind(col1, wCI)
sw <- as.data.frame(model.avg_MAMMAL$sw)
col1 <- rep("Mammals", length(model.avg_MAMMAL$sw))
sw_MAMMAL <- cbind(col1, sw); names(sw_MAMMAL) <- c("Type", "SW")

summary(best.model_MAMMAL)

MpcLOG_MAMMAL_retained<-MCMCglmm(Log.effect.size.g ~ Decline.speed + Decline.scale + Foraging.guild ,
                                 random=~Publication + animal, 
                                 family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                                 prior=prior, data=a, nitt=500000,burnin=1000,thin=500)  #with phylogeny
plot(MpcLOG_MAMMAL_retained)
qqplot(MpcLOG_MAMMAL_retained)
summary(MpcLOG_MAMMAL_retained)
summary(top.models_MAMMAL)

#__(m4) calculate error and distribution of data -----
#__ funnel plots 
plot(aqua$Log.effect.size.g, I(1/aqua$sd)) # this makes the funnel plot of slope (rate of change in days/year) and precision (1/SE)

#__ error calcs 
randomtest_MAMMAL <-MCMCglmm(Log.effect.size.g ~ 1, 
                             random=~Publication + animal + idh(sd):units, family="gaussian",
                             ginverse=list(phylo=inv.phylo$Ainv),
                             prior=priorSD, data=aqua,nitt=50000,burnin=1000,thin=500)
summary(randomtest_MAMMAL)
par(mfrow = c(1,2)); hist(mcmc(randomtest_MpcLOGaquaMAMMAL_SD$VCV)[,"Publication"]); hist(mcmc(randomtest_MpcLOGaquaMAMMAL_SD$VCV)[,"animal"]); par(mfrow=c(1,1)) # Reset the plot panel back to single plots
plot(randomtest_MpcLOGaquaMAMMAL_SD$Sol) #tests for model convergence.  The trace is like a time series of what your model did while it was running and can be used to assess mixing (or convergence), while the density is like a smoothed histogram of the estimates of the posterior distribution that the model produced for every iteration of the model. To make sure your model has converged, the trace plot should look like a fuzzy caterpillar.  
plot(randomtest_MpcLOGaquaMAMMAL_SD$VCV) #tests relationship for random factors.  It looks like some of the variances of the random effects haven’t mixed very well at all. The effective sample size is also very small. Maybe we could improve this by increasing the number of iterations, but because the chain seems to be stuck around zero, it looks like we’ll need to use a stronger prior than the default.
xsim <- simulate(randomtest_MpcLOGaquaMAMMAL_SD) 
plot(aqua$Log.effect.size.g, I(1/aqua$sd))
points(xsim, I(1/aqua$sd), col = "red") 


# (r) REPTILES -----
# __(r1) Reptile DATA -----
# notes (1) S=0, T=1 < ALL REPTILE DATA SPATIAL, (2) reptile data lacked variation for decline.permanence (3) make sure that species names have an underscore, not space
# (4) ordered from ground to trees - ground is baseline
setwd("C:/Users/20043940/Dropbox/TF data/Tree declines/Phylogenetic corrections/ALL reptiles")
# DATA
black<- read.csv("allDATA.csv")
blackREPS <- black[which(black$Class=='Reptile'),]
blackREPS$animal <- as.factor(blackREPS$animal) # Pedigree function needs column called "animal"
blackREPS$Foraging.guild <- factor(x=blackREPS$Foraging.guild, levels=c("GRO", "BAR")) #to set order - comparison 
blackREPS$Primary.diet <- factor(x=blackREPS$Primary.diet, levels=c("O", "I")) #to set order - comparison excludes "V" = combined with O 
blackREPS$Decline.cause<- as.factor(blackREPS$Decline.cause)
blackREPS$Temporal.or.Spatial<- scale(blackREPS$Temporal.or.Spatial) # S=0, T=1
blackREPS$Decline.permanence<- scale(blackREPS$Decline.permanence) #use as continuous, not categorical
blackREPS$Decline.scale<- scale(blackREPS$Decline.scale) #use as continuous, not categorical
blackREPS$Decline.speed<- scale(blackREPS$Decline.speed)
# PHYLOGENY
phylo<-read.nexus("phylo.trees.nex")
summary.phylo(phylo) # checks if tree has imported correctly (number of tips)
phylo$tip.label
phylo$node.label <- NULL
#phylo2=drop.tip(phylo, tip=c(2:3)) # TRIMMING OPTION to keep 1 rotifera (to root) and other acantho
is.rooted(phylo) # checks if the tree is rooted (necessary)
phylo=force.ultrametric(phylo, method="extend") # to make the tree ultrametric
is.ultrametric(phylo) # checks if the tree is ultrametric (necessary)    TRUE
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)
plot(phylo)

# __(r2) Sensitivity analysis on global models testing effect of (1) raw or log-transformed g, (2) the inclusion of species as a random factor, (3) the inclusion of phylogenetic correction (PC) as a random factor-----
priorSD=list(R=list(R1=list(V=1,nu=0.002),R2=list(V=1,nu=0.002),R3=list(V=1,nu=0.002),
                    R4=list(V=1,nu=0.002),R5=list(V=1,nu=0.002),R6=list(V=1,nu=0.002)), 
             G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))  
prior_noSD=list(R=list(R1=list(V=1,nu=0.002),R2=list(V=1,nu=0.002),R3=list(V=1,nu=0.002),
                       R4=list(V=1,nu=0.002),R5=list(V=1,nu=0.002),R6=list(V=1,nu=0.002)), 
                G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))
Mlmer_REPTILE <-lmer(yi ~  Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild  +
                       (1|Publication), data=blackREPS)
Mglmm_REPTILE <- MCMCglmm(yi ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                          random=~Publication + animal, singular.ok=TRUE,
                          prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackREPS)       #  species as a random factor, but no phylogeny
Mpc_REPTILE <-MCMCglmm(yi ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild ,
                       random=~Publication + animal, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                       prior=prior,data=blackREPS, nitt=50000,burnin=1000,thin=500)  #with phylogeny    
Mglmm_REPTILE_SD <- MCMCglmm(yi ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                          random=~Publication + animal + idh(vi):units, singular.ok=TRUE,
                          prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackREPS)       #  species as a random factor, but no phylogeny
Mpc_REPTILE_SD <-MCMCglmm(yi ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild ,
                       random=~Publication + animal + idh(vi):units, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                       prior=NULL,data=blackREPS, nitt=50000,burnin=1000,thin=500)  #with phylogeny    

# USING HEDGES G TRANSFORMED Log.effect.size.g
MlmerLOG_REPTILE <-lmer(log_yi ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild + 
                          (1|Publication), data=blackREPS)

MglmmLOG_REPTILE <- MCMCglmm(log_yi~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                             random=~Publication + animal, singular.ok=TRUE,
                             prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackREPS) #  species as a random factor, but no phylogeny

MpcLOG_REPTILE<-MCMCglmm(log_yi ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                         random=~Publication + animal, 
                         family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                         prior=prior, data=blackREPS, nitt=50000,burnin=1000,thin=500)  #with phylogeny
MglmmLOG_REPTILE_SD <- MCMCglmm(log_yi~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                             random=~Publication + animal + idh(vi):units, singular.ok=TRUE,
                             prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackREPS) #  species as a random factor, but no phylogeny

MpcLOG_REPTILE_SD<-MCMCglmm(log_yi ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                         random=~Publication + animal + idh(vi):units, 
                         family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                         prior=NULL, data=blackREPS, nitt=50000,burnin=1000,thin=500)  #with phylogeny

DIC(Mlmer_REPTILE,Mglmm_REPTILE,Mpc_REPTILE,Mglmm_REPTILE_SD,Mpc_REPTILE_SD,    MlmerLOG_REPTILE,MglmmLOG_REPTILE,MpcLOG_REPTILE,     MglmmLOG_REPTILE_SD, MpcLOG_REPTILE_SD)
AIC(Mlmer_REPTILE,Mglmm_REPTILE,Mpc_REPTILE,Mglmm_REPTILE_SD,Mpc_REPTILE_SD,    MlmerLOG_REPTILE,MglmmLOG_REPTILE,MpcLOG_REPTILE,     MglmmLOG_REPTILE_SD, MpcLOG_REPTILE_SD)
BIC(Mlmer_REPTILE,Mglmm_REPTILE,Mpc_REPTILE,Mglmm_REPTILE_SD,Mpc_REPTILE_SD,    MlmerLOG_REPTILE,MglmmLOG_REPTILE,MpcLOG_REPTILE,     MglmmLOG_REPTILE_SD, MpcLOG_REPTILE_SD)

# __(r3) Dredge on pc-corrected model -----
# NOTE: the model without phylogenetic corection (PC) was marginally stronger than the one with PC.  Carried out sensitivity analyses comparing with and without PC.
MCMCglmm.updateable<- updateable(MCMCglmm)  
gm_REPTILE_noPC<- MCMCglmm.updateable(log_yi ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                                      random=~Publication + animal + idh(vi):units, singular.ok=TRUE,
                                      prior=NULL, nitt=50000,thin=250,burnin=1000, data=blackREPS)
options(na.action = "na.fail"); dredge.gm_REPTILE_noPC<- dredge(gm_REPTILE_noPC, rank="DIC", evaluate=T); options(na.action = "na.omit")
dredge.gm_REPTILE_noPC   # note - had to remove decline.permanence to get the model to run
write.csv(dredge.gm_REPTILE_noPC , "dredge_blackREPS.csv")


gm_REPTILE<- MCMCglmm.updateable(log_yi ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                                      random=~Publication + animal + idh(vi):units, 
                                      family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                                      prior=NULL, data=blackREPS, nitt=50000,burnin=1000,thin=500)
options(na.action = "na.fail"); dredge.gm_REPTILE<- dredge(gm_REPTILE, rank="DIC", evaluate=T); options(na.action = "na.omit")
dredge.gm_REPTILE   # note - had to remove decline.permanence to get the model to run
write.csv(dredge.gm_REPTILE, "dredge_blackREPS_mod_v.csv")

retained_REPTILE<- MCMCglmm(log_yi ~ Decline.speed,
                            random=~Publication + animal + idh(vi):units, 
                            family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                            prior=NULL, data=blackREPS, nitt=50000,burnin=1000,thin=500)
MpcLOG_REPTILE_retained<-MCMCglmm(Log.effect.size.g_REPTILE ~ Decline.permanence	 + Foraging.guild + Decline.scale, 
                                  random=~Publication + animal, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                                  prior=prior, data=a, nitt=50000,burnin=1000,thin=500)  

summary(retained_REPTILE)
summary(MpcLOG_REPTILE_retained)

top.models_REPTILE <- get.models(dredge.gm_REPTILE, subset = delta<2)
model.avg_REPTILE <- summary(model.avg(top.models_REPTILE))
MA.est.table_REPTILE <- model.avg_REPTILE$coefmat.full
CI <- as.data.frame(confint(model.avg_REPTILE)) ## calc confidence intervals
df <- as.data.frame(MA.est.table_REPTILE)
wCI <-cbind(df, CI)
col1 <- rep("REPTILEs", length(df$Estimate))
df_REPTILE <- cbind(col1, wCI)
swREPTILE <- as.data.frame(model.avg_REPTILE$sw)
col1 <- rep("REPTILEs", length(model.avg_REPTILE$sw))
sw_REPTILE <- cbind(col1, sw); names(sw_REPTILE) <- c("Type", "SW")

#__(r4) calculate error and distribution of data -----
# funnel plots
plot(blackREPS$yi, I(1/blackREPS$vi)) # this makes the funnel plot of slope (rate of change in days/year) and precision (1/SE)

# error calcs 
randomtest_REPTILEerror <-MCMCglmm(Log.effect.size.g ~ 1, 
                                          random=~Publication + animal + idh(sd):units, family="gaussian",
                                          ginverse=list(phylo=inv.phylo$Ainv),
                                          prior=priorSD, data=a,nitt=50000,burnin=1000,thin=500)
summary(randomtest_REPTILEerror)
par(mfrow = c(1,2)); hist(mcmc(randomtest_REPTILEerror$VCV)[,"Publication"]); hist(mcmc(randomtest_REPTILEerror$VCV)[,"animal"]); par(mfrow=c(1,1)) # Reset the plot panel back to single plots
plot(randomtest_REPTILEerror$Sol) #tests for model convergence.  
plot(randomtest_REPTILEerror$VCV) #tests relationship for random factors.  
xsim <- simulate(randomtest_REPTILEerror) 
plot(a$Log.effect.size.g, I(1/a$sd))
points(xsim, I(1/a$sd), col = "red") 


# (i) INVERTS -----
# __(i1) Invert DATA -----
setwd("C:/Users/20043940/Dropbox/TF data/Tree declines/Phylogenetic corrections/ALL inverts")
black<- read.csv("allDATA.csv")
blackINVERTS <- black[which(black$Taxon=='Arthropod'),]
# notes (1) S=0, T=1 <- only S data for invert, (2) make sure that species names have an underscore, not space
# (3) ordered from ground to trees - ground is baseline
# DATA
blackINVERTS$Foraging.guild <- factor(x=blackINVERTS$Foraging.guild, levels=c("GRO", "LOW","CAN")) #to set order - comparison 
blackINVERTS$animal <- as.factor(blackINVERTS$animal) # Pedigree function needs column called "animal"
blackINVERTS$Decline.cause<- as.factor(blackINVERTS$Decline.cause)
blackINVERTS$Decline.permanence<- scale(blackINVERTS$Decline.permanence) #use as continuous, not categorical
blackINVERTS$Decline.scale<- scale(blackINVERTS$Decline.scale) #use as continuous, not categorical
blackINVERTS$Decline.speed<- scale(blackINVERTS$Decline.speed)
# PHYLOGENY
phylo<-read.nexus("phylo.trees.nex") #INVERTS
summary.phylo(phylo) # checks if tree has imported correctly (number of tips)
phylo$tip.label
phylo$node.label <- NULL
#phylo2=drop.tip(phylo, tip=c(2:3)) # TRIMMING OPTION to keep 1 rotifera (to root) and other acantho
is.rooted(phylo) # checks if the tree is rooted (necessary)
phylo=force.ultrametric(phylo, method="extend") # to make the tree ultrametric
is.ultrametric(phylo) # checks if the tree is ultrametric (necessary)    TRUE
inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)
plot(phylo)

# __(i2) Sensitivity analysis on global models testing effect of (1) raw or log-transformed g, (2) the inclusion of species as a random factor, (3) the inclusion of phylogenetic correction (PC) as a random factor-----
Mlmer_INVERT <-lmer(g ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild  +
                      (1|Publication), data=a)
Mglmm_INVERT <- MCMCglmm(g ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                         random=~Publication + animal, singular.ok=TRUE,
                         prior=NULL, nitt=50000,thin=250,burnin=1000, data=a)       # Bird species as a random factor, but no phylogeny
Mpc_INVERT <-MCMCglmm(g ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                      random=~Publication + animal, family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                      prior=prior,              data=a, nitt=50000,burnin=1000,thin=500)  #with phylogeny    
MlmerLOG_INVERT <-lmer(Log.effect.size.g ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild  +
                         (1|Publication), data=a)
MglmmLOG_INVERT <- MCMCglmm(Log.effect.size.g ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                            random=~Publication + animal, singular.ok=TRUE,
                            prior=NULL, nitt=50000,thin=250,burnin=1000, data=a) # Bird species as a random factor, but no phylogeny
MpcLOG_INVERT<-MCMCglmm(Log.effect.size.g ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                        random=~Publication + animal, 
                        family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                        prior=prior, data=a, nitt=50000,burnin=1000,thin=500)  #with phylogeny
# with standard deviation
MglmmLOG_INVERT_wSD <- MCMCglmm(Log.effect.size.g ~ Decline.permanence	+ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                                random=~Publication + animal  + idh(sd):units, singular.ok=TRUE,
                                prior=priorSD, nitt=50000,thin=250,burnin=1000, data=a) #  species as a random factor, but no phylogeny

MpcLOG_INVERT_wSD <-MCMCglmm(Log.effect.size.g ~ Decline.scale	+ Decline.speed + Primary.diet + Foraging.guild,
                             random=~Publication + animal + idh(sd):units, 
                             family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                             prior=priorSD, data=a, nitt=50000,burnin=1000,thin=500)  #with phylogeny

DIC(Mlmer_INVERT,Mglmm_INVERT,Mpc_INVERT,MlmerLOG_INVERT,MglmmLOG_INVERT,MpcLOG_INVERT,MglmmLOG_INVERT_wSD, MpcLOG_INVERT_wSD)


# __(i3) Dredge on pc-corrected model -----
MCMCglmm.updateable<- updateable(MCMCglmm)  
gm_INVERT_pc<- MCMCglmm.updateable(Log.effect.size.g ~ Decline.scale	+ Decline.permanence + Decline.speed + Primary.diet + Foraging.guild,
                                random=~Publication + animal + idh(sd):units, 
                                family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),singular.ok=TRUE,
                                prior=priorSD, data=a, nitt=50000,burnin=1000,thin=500)
options(na.action = "na.fail"); dredge.gm_INVERT_pc <- dredge(gm_INVERT_pc, rank="DIC"); options(na.action = "na.omit")
dredge.gm_INVERT_pc

best.models_INVERT <- get.models(dredge.gm_INVERT, subset = delta==0) 
top.models_INVERT <- get.models(dredge.gm_INVERT, subset = delta<2)
model.avg_INVERT <- summary(model.avg(top.models_INVERT))
MA.est.table_INVERT <- model.avg_INVERT$coefmat.full
CI <- as.data.frame(confint(model.avg_INVERT)) ## calc confidence intervals 
df <- as.data.frame(MA.est.table_INVERT)
wCI <-cbind(df, CI)
col1 <- rep("INVERTs", length(df$Estimate))
df_INVERT <- cbind(col1, wCI)
sw <- as.data.frame(model.avg_INVERT$sw)
col1 <- rep("INVERTs", length(model.avg_INVERT$sw))
sw_INVERT <- cbind(col1, sw); names(sw_INVERT) <- c("Type", "SW")

#__(i4) calculate error and distribution of data -----
# funnel plots
plot(a$Log.effect.size.g, I(1/a$sd)) # this makes the funnel plot of slope (rate of change in days/year) and precision (1/SE)

#__ error calcs 
randomtest_INVERT <-MCMCglmm(Log.effect.size.g ~ 1, 
                             random=~Publication + animal + idh(sd):units, family="gaussian",
                             ginverse=list(phylo=inv.phylo$Ainv),
                             prior=priorSD, data=a,nitt=50000,burnin=1000,thin=500)
summary(randomtest_INVERT)

# THE END -----