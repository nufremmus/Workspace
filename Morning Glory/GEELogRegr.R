# we add more specifications in gee(), such as link function = "logit" or family = "binomialbinomial"
library(MASS)

#library for estimating equations
library("gee")

#library for importing a newick tree
library("ape")
MyTree <- read.tree("~/Google Drive/correlated association/ITS_w_outgrps_cured.newick")

#library for comparative phylogeny
library("phytools")

#library for read_excel
library("readxl")
mydata <- read_excel("~/Google Drive/correlated association/Alkaloid Seed Mass Table for Doug Aug 31 2016_clean.xlsx", 
                                                                            na = "null")

#lib for inv.logit
library("boot")
############################################################
#        clean data: since there more data in excel than 
#        in newick, we need to remove the extra data in excel
############################################################

I = match(MyTree$tip.label,mydata$Species)
I_remain_data = I[which(!is.na(I))]
#remove extra data
mydata <-mydata[I_remain_data,]

############################################################
#        pruning MyTree
############################################################
Index_extra_tree = which(is.na(I))

# Species_extra = MyTree$tip.label[Index_extra]
pruned.tree<-drop.tip(MyTree,MyTree$tip.label[Index_extra_tree])

# verify
match(pruned.tree$tip.label,mydata$Species)

colnames(mydata) = c("Species" , "Num_Acc","Num_Alk","SeedMass")

#add a new feature: the chance of having Ergot
mydata["Alk_P"] <- mydata["Num_Alk"]/mydata["Num_Acc"]

#sort mydata 
mydata<-mydata[match(pruned.tree$tip.label,mydata$Species),]

############################################################
#       generalized linear models
############################################################

library(nlme)

#phylogenetic correlation assuming the continuous phenotypes evolved according to Brownian motion
sigma<- vcv.phylo(pruned.tree, cor=TRUE)

gee_model1 = gee(Alk_P ~ SeedMass,id = c(rep(0,82)), data = mydata, corstr = "fixed",family = binomial, R = sigma )
beta0 = gee_model1$coefficients[1]
beta1 = gee_model1$coefficients[2]

# assume beta_1 is zero
gee_model1_null = gee(Alk_P ~ 1,id = c(rep(0,82)), data = mydata, corstr = "fixed",family = binomial, R = sigma )
# get the chance of having Ergot for morning glory family under the null hypothesis that beta_1 is zero
beta0_null <- gee_model1_null$coefficients[1]
p_null <- as.numeric(inv.logit(beta0_null))

############ Simulation to get distribution of beta under the null hypothesis ##############
data_sim <- mydata
mu <- c(rep(0,82))
# sigma<- vcv.phylo(pruned.tree, cor=TRUE)
Nsim = 2000

# the normalized numer of Ergot for each taxa follows normal distribution (0, sigma) 
#DOUBLE CHECK, right now the error has unit variances for each species
noise <- mvrnorm(Nsim, mu = mu, Sigma = sigma )
# sim_Alk is the simulated number of Ergot each taxa should have under the null hypothesis
sim_Alk_P =   matrix( c(rep(0,Nsim*82)),nrow=Nsim,ncol=length(mydata$Species)) 
sim_Alk_num =   matrix( c(rep(0,Nsim*82)),nrow=Nsim,ncol=length(mydata$Species)) 
# data_sim$Alk =   matrix( c(rep(0,Nsim*82)),nrow=Nsim,ncol=length(mydata$Species)) 
for(i in 1:82){
#### add noise to beta_null then use inv.logit link function to transform that into simulated presence of alkaloids
  #simulated chance of having alkaloids
  sim_Alk_P[,i] = inv.logit(noise[,i]+beta0_null)
  #simulate the number of sample that have alkaloids according to binomial distribution
  sim_Alk_num[,i] = rbinom(Nsim,mydata$Num_Acc[i],sim_Alk_P[,i])
}

# generate beta1 distribution under the null hypothesis
beta_sim<-c(rep(0,Nsim))
count_exception = 0
for(j in 1:Nsim){
  # When divergences happen in gee(), that loop is skipped
  error_beta_sim<-tryCatch({
  #get the simulated chance of having Ergot for each species
  data_sim$Alk_P = sim_Alk_num[j,]/data_sim$Num_Acc 
  # get the simulated beta_1 from the simulated data
  beta_sim[j] = gee(Alk_P ~ SeedMass,id = c(rep(0,82)), data = data_sim, family = binomial, corstr = "fixed", R = sigma )$coefficient[2]},
  error=function(e){e}
)
  if(inherits(error_beta_sim,"error")){
    count_exception = count_exception+1
  }
}

#get rid of the values of beta_sim that diverged
# beta_sim = beta_sim[beta_sim!=0]
# get p-value
mean(abs(beta_sim) > abs(beta1))
