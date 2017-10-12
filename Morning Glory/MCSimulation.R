#changed gee parameters: family = binomial and made sure sim_alk are nonnegative

#library for estimating equations
library("gee")

#for generating random multivariate normal r.v.'s
library(mvtnorm)

#library for importing a newick tree
library("ape")
MyTree <- read.tree("~/Google Drive/correlated association/ITS_w_outgrps_cured.newick")

#library for comparative phylogeny
library("phytools")

#library for read_excel
library("readxl")

#library for inv.logit
library("boot")

#library for nonlinear models
library(nlme)

#lib for excel
library(readxl)
mydata <- read_excel("~/Google Drive/correlated association/Alkaloid Seed Mass Table for Doug Aug 31 2016_clean.xlsx",
                     na = "null")
############################################################
#        clean data: more data in excel than in newick
############################################################
I = match(MyTree$tip.label,mydata$Species)
I_remain_data = I[which(!is.na(I))]
mydata <-mydata[I_remain_data,]

# mydata <- mydata[match(MyTree$tip.label,mydata$Species),]
colnames(mydata) = c("Species" , "Num_Acc","Num_Alk","SeedMass")
mydata["Alk_P"] <- mydata["Num_Alk"]/mydata["Num_Acc"]

############################################################
#        pruning MyTree
############################################################
Index_extra_tree = which(is.na(I))
pruned.tree<-drop.tip(MyTree,MyTree$tip.label[Index_extra_tree])

#sort mydata so the species match the MyTree tips' labels
mydata<-mydata[match(pruned.tree$tip.label,mydata$Species),]
# get maximum phylogeny for the given topology using Alk tips

MP.tree<-pruned.tree
MP.tree$Alk_P<- mydata$Alk_P[match(pruned.tree$tip.label,mydata$Species)]
MP.tree$Alk_P[which(MP.tree$Alk_P>0)] = 1

#approximate ancestral nodes with maximum parsimony
library("mvSLOUCH")
MP.ouchtree<- ape2ouch(MP.tree)
#correct the names of the internal nodes
MP.ouchtree@nodelabels[1:(MP.ouchtree@nnodes-MP.ouchtree@nterm)]<-as.character(
  1:(MP.ouchtree@nnodes-MP.ouchtree@nterm))

# root = 2 gives us root 1 ; 2nd state in level (0 1) is 1
# acctran accelerates transformation along evolution, thus realizing our assumption that if a child has Ergot, their parent are more likely to have Ergot as well
# The experiment is different from our assumption that most ancestral nodes will be 1.....
# The reality: most ancestral nodes are 0 and setting root = 1 does not affect the ancestral nodes close to root being 0
# NOTE: acctran or deltran does not affect how we calculate count_flip
# root = 2, set jump to a negative number
# root = 1, set jump to a positive number
MP.fitch<- fitch.mvsl(MP.ouchtree, MP.tree$Alk_P, deltran = FALSE, acctran = TRUE, root = 1)
MP1.ouchtree<-MP.ouchtree
plot(MP.ouchtree)
for(i in 1:length(MP1.ouchtree@nodelabels)){
  MP1.ouchtree@nodelabels[i] = as.character(MP.fitch[i])
}
MP1.ouchtree@nodelabels[1] = "1"

#simulate the SeedMass
library(MASS)
mu_SM = c(rep(mean(mydata$SeedMass),82))
#should I add sig_SM in the simulation
sig_SM = var(mydata$SeedMass)
sigma<- vcv.phylo(pruned.tree, cor=TRUE) #phylogenetic correlation according to brownian motion
Nsim_SM = 10
#assumed that each species has the same variance
sim_SM <- mvrnorm(Nsim_SM, mu = mu_SM, Sigma = sigma*sig_SM )

#calculate the number of jumps for each species
count_flip = c(rep(0,163)) #only uses the latter 82 entries so the index is consistent
for(k in 82:163){ #1:81 are for internal nodes
  ancestor_curr = k
  while(ancestor_curr!=1){
    ancestor_prev = ancestor_curr
    ancestor_curr = as.numeric(MP.ouchtree@ancestors[ancestor_prev])
    fitch_curr = as.numeric(MP.fitch[ancestor_curr])
    fitch_prev = as.numeric(MP.fitch[ancestor_prev])
    # ergot is gained during evolution: 1 if value = 0; 2 if value = 1.
    if((fitch_curr==1) & (fitch_prev==2) ){
      count_flip[k] = count_flip[k]+1;
    }
    # ergot is lost during evolution
    if((fitch_curr==2) & (fitch_prev==1) ){
      count_flip[k] = count_flip[k]-1;
    }
    }}

# the constant we add to the seed mass if that species gain Ergot in their evolution
jump = 60

# if there are even number of jumps make it zero jump; odd then one jump
sim_SM[,which(count_flip[82:163]==1)] = sim_SM[,which(count_flip[82:163]==1)]+jump # adding jumps to the species
#Exception handling
if(length(which((sim_SM<0) == TRUE)) > 0 ){
  cat("WRONG simulation")
  Sys.sleep(5)
}
###########################################################
#       generalized linear models
############################################################
#initialization
pvalue = rep(0,Nsim_SM) #stores pvalues for each round of simulation
length_sim = rep(0,Nsim_SM) #stores the number of simulated beta1 under H0 in each simulation
count_exception_beta1 = NULL #stores the round index when divergences of gee() for beta1 happens

#repeat simulations Nsim_SM times
for(round_SM in 1:Nsim_SM){
    #stores the loop index when divergences of gee() for beta_sim happens
    count_exception_beta_sim = NULL
    #ERROR HANDLING for getting beta1 from gee; if gee() diverges that loop is skipped
    Error_beta1 <- tryCatch(
      {  mydata["sim_SM"]<-sim_SM[round_SM,]
      beta1 = gee(Alk_P ~ sim_SM,id = c(rep(0,82)), data = mydata, corstr = "fixed",family=binomial, R = sigma )$coefficients[2]
      beta_null = gee(Alk_P ~ 1,id = c(rep(0,82)), data = mydata, corstr = "fixed",family=binomial, R = sigma )$coefficients[1]
      p_null <- inv.logit(beta_null)
      p_null = as.numeric(p_null)
      },
      error=function(e) e)

  # stores the round index for beta1 when gee() diverges
  if(inherits(Error_beta1, "error")){
    count_exception_beta1 = c(count_exception_beta1,round_SM)
    next}
  
  #continue simulation for beta1 under null hypothesis
  ############ Simulation to get distribution of beta ##############
  data_sim <- mydata
  mu <- c(rep(0,82))
  Nsim = 500
  noise <- mvrnorm(Nsim, mu = mu, Sigma = sigma )
  sim_Alk_num =   matrix( c(rep(0,Nsim*82)),nrow=Nsim,ncol=length(mydata$Species))
  sim_Alk_p   =   matrix( c(rep(0,Nsim*82)),nrow=Nsim,ncol=length(mydata$Species))
  
  #simulate number of sample with ergot alkaloids for each species
  for(i in 1:82){
    sim_Alk_p[,i]   = inv.logit(noise[,i]+beta_null)
    sim_Alk_num[,i] = rbinom(Nsim,data_sim$Num_Acc[i],sim_Alk_p[,i]) 
  }

  # generate beta distribution under the null hypothesis
  beta_sim<-c(rep(0,Nsim))
  for(j in 1:Nsim){
    # if gee() diverges that loop is skipped
    Error_beta_sim <- tryCatch(
    {data_sim$Alk_P = sim_Alk_num[j,]/data_sim$Num_Acc
    beta_sim[j] = gee(Alk_P ~ sim_SM,id = c(rep(0,82)), data = data_sim, family=binomial,corstr = "fixed", R = sigma )$coefficient[2]},
    error=function(e) e
    )
    
    #stores the loop index for beta_sim when gee() diverges
    if (inherits(Error_beta_sim, "error")){
      count_exception_beta_sim = c(count_exception_beta_sim,j)
      next}
  }#end for loop for j
  
  #clean beta_sim to get rid of ones diverged
  if (!is.null(count_exception_beta_sim)){
     beta_sim = beta_sim[-count_exception_beta_sim]
  }
  
  # get p-value
  pvalue[round_SM] = mean(abs(beta_sim) > abs(beta1))
  length_sim[round_SM] = length(beta_sim)
}#end round_SM

 pvalue_effective = pvalue[-count_exception_beta1]
 save(pvalue_effective,count_exception_beta1,pvalue,file = "pvalue.RData")
