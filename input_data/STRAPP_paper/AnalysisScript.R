library(BAMMtools);
source("traitDependentBAMM.R")


#fish data analysis for figure 3a

#read in BAMM diversification result
v <- read.tree('fish_naturecomm.tre');
events_cr <- 'fish2cr_event_data_250m.csv'
ed_cr <- getEventData(v, events_vr, burnin=0.1)

#read in BAMM morphology--body size result and get the rate of body size evolution for each species 

vmorph <- read.tree('fishtree_morph.txt');
eventsmorph <- 'fishsize_event_data.txt';
edmorph <- getEventData(vmorph, eventsmorph, burnin=0.1, nsamples=1000, type ='trait');
traits <- getTipRates(edmorph)$beta.avg;

#Permutation analysis
fish.permu<-traitDependentBAMM(ed_cr,traits,10000,return.full = T,logrates = T,two.tailed = T,method='s')

#check the difference between observed correlation and permuted correlation
fish_diff<-abs(fish.permu$obs.corr)-abs(fish.permu$null)
hist(fish_diff)

#Type I error rate on the fish pylogeny, figure 3b

#simulate traits

simulateTraitsMatrix <- function(phy, sigma2 = 1, rootstate=0, reps=100, verbose=F){  
  states <- matrix(0, nrow=nrow(phy$edge), ncol=reps);  
  root <- length(phy$tip.label) + 1;
  isRootEdge <- phy$edge[,1] == root;
  rootstate <- rep(rootstate, reps);  
  states[isRootEdge, ][1,] <- rootstate + rnorm(reps, sd=sqrt(phy$edge.length[isRootEdge][1] * sigma2));
  states[isRootEdge, ][2,] <- rootstate + rnorm(reps, sd=sqrt(phy$edge.length[isRootEdge][2] * sigma2));  
  for (i in 1:nrow(phy$edge)){
    if (verbose){
      cat(i, '\n');
    }
    node <- phy$edge[i,2];
    par <- phy$edge[i,1];
    if (par != root){
      par.state <- states[phy$edge[,2] == par, ];
      states[i,] <- par.state + rnorm(reps, sd=sqrt(phy$edge.length[i] * sigma2));
    }
    
  }
  rownames(states) <- phy$edge[,2];
  phy$states <- states[as.character(1:length(phy$tip.label)), ];  
  return(phy)
}
ttree<-drop.tip(v,tip = v$tip.label[! v$tip.label %in% names(traits)])
ttree<-simulateTraitsMatrix(ttree,reps = 1000)

#perform permutation on each of the simulated trait data
p1<-sapply(1:dim(ttree$states)[2],function(i){
  lt<-ttree$states[,i]
  names(lt)<-ttree$tip.label
  l<-traitDependentBAMM(ed_cr,lt,1000,return.full = F,logrates = T,two.tailed = T)
  l$p.value
})

#calculate the type I error rate.
sum(p1<=0.05)/length(p1)


#subset analysis on fish morpholoy data Figure 4a
subsample<-rev(c(25,50,100,200,400,800,1600,3200 )) #the number of taxa to sub sample
rep<-1:100 # for each subset size, do 100 replicates
p<-matrix(0,nrow = length(subsample),ncol = length(rep))
for (i in 1:length(rep)){
  t<-traits
  for (s in 1:length(subsample)){
    t<-sample(t,size = subsample[s],replace = F)
    l<-traitDependentBAMM(ed_cr,t,1000,return.full = F,logrates = T,two.tailed = T)
    p[s,i]<-l$p.value
  }
  cat(i, format(Sys.time(), "%a %b %d %X %Y"),"\n")
}
row.names(p)<-as.character(subsample)
p5<-apply(X = p,MARGIN = 1,function(i){sum(i<=0.05)}) # the number of replicates out 100 detected a significant correlation at 0.05 level
p1<-apply(X = p,MARGIN = 1,function(i){sum(i<=0.01)}) # the number of replicates out 100 detected a significant correlation at 0.01 level

#subset analysis with simulated Brownian traits on fish morphology

p<-matrix(0,nrow = length(subsample),ncol = length(rep))
for (i in 1:length(rep)){
  t<-ttree$states[,i]
  names(t)<-ttree$tip.label
  for (s in 1:length(subsample)){
    t<-sample(t,size = subsample[s],replace = F)
    l<-traitDependentBAMM(ed_cr,t,1000,return.full = F,logrates = T,two.tailed = T)
    p[s,i]<-l$p.value
  }
  cat(i, Sys.time(),"\n")
}
p5<-apply(X = p,MARGIN = 1,function(i){sum(i<=0.05)}) # the number of replicates out 100 detected a significant correlation at 0.05 level-- false positives
p1<-apply(X = p,MARGIN = 1,function(i){sum(i<=0.01)}) # the number of replicates out 100 detected a significant correlation at 0.01 level-- false positives

#analysis on bird dichromatism

#dichromatism as a continous trait
dm<-read.table("Amenta_data.txt",header=T,sep=" ",stringsAsFactors = F)
trait<-dm$pca
names(trait)<-dm$species
trait<-trait[! is.na(trait)]

bird_event<-"hack2vr_event_data_250m.csv"
bird_treefile<-"birds_hackett.tre"
bird_tree<-read.tree(bird_treefile)
ed<-getEventData(bird_tree,eventdata = bird_event,burnin = 0.1)

bird_dm<-traitDependentBAMM(ephy = ed,traits = trait,reps = 100,return.full = T,logrates = T,two.tailed = T)
#check the difference in absolute spearman's correlation coefficients
bird_dm_diff<-abs(bird_dm$obs.corr)-abs(bird_dm$null)

#dichromatism as binary trait Figure 5a
bitrait<-sapply(trait, function(x){if(x>median(trait)){1}else{0}}) #use 50% cutoff,  i.e. half of the species are dichromatic
bird_bidm<-traitDependentBAMM(ephy = ed,traits = bitrait,reps = 10000,return.full = T,logrates = T,method="m", two.tailed = T)
bitrait<-sapply(trait, function(x){if(x>quantile(trait,probs = 0.80)){1}else{0}}) #use 80% cutoff, i.e., 20% of the species are dichromatic
bird_bidm2<-traitDependentBAMM(ephy = ed,traits = bitrait,reps = 10000,return.full = T,logrates = T,method="m", two.tailed = T)

#type I error rate with binary traits simulated on the bird phylogeny Figure 5b

library(phylolm)
library(phyclust)

birdtree<-as.phylo(ed)
birdtree<-drop.tip(birdtree,tip = birdtree$tip.label[! birdtree$tip.label %in% names(trait)])
birdtree<-rescale.rooted.tree(birdtree, scale.height = 1/get.rooted.tree.height(birdtree)) #scale the root height of the tree for simulation
trait<-trait[birdtree$tip.label]

rep=100;
beta<-c(0,log(0.25)) #two different stationary frequencies, 0 correspond to 50:50, log(0.25) corresponds to 20:80
alpha<-c( exp(-1),1,exp(1),9) # four different phylogenetic signal (1, 0, -1, and -log(9)), from strong to weak

bi_p1<-data.frame(beta=numeric(0),alpha=numeric(0),typeI=numeric(0),stringsAsFactors = F)
p<-matrix(0, ncol = rep,nrow = 0)

for (b in 1:length(beta)){
  for(a in 1:length(alpha)){
    cat(beta[b],alpha[a],"\n")
    sim_bin<-rbinTrait(n=rep,phy = birdtree,beta = beta[b],alpha = alpha[a]) # simulate binary trait
    
    #get rid of sets of simulations where no variation is observed between traits
    x<-apply(sim_bin,2,function(x) sum(x)/dim(sim_bin)[1])
    sim_bin<-sim_bin[,x>0 & x<1,drop = FALSE]
    while (dim(sim_bin)[2]<rep){
      y<-rbinTrait(n=rep,phy = birdtree,beta = beta[b],alpha = alpha[a])
      
      x<-apply(y,2,function(x) sum(x)/dim(y)[1])
      y<-y[,x>0 & x<1,drop = FALSE]
      sim_bin<-cbind(sim_bin,y)
    }
    
    pp<-sapply(1:rep,function(i){
      lt<-sim_bin[,i]
      names(lt)<-row.names(sim_bin)
      l<-traitDependentBAMM(ed,lt,1000,return.full = F,logrates = T,two.tailed = T,method = "m")
      l$p.value
    })
    p<-rbind(p,pp)
    bi_p1<-rbind(bi_p1,data.frame(beta=beta[b],alpha=alpha[a],typeI=sum(pp<=0.05)/rep,stringsAsFactors = F))
  }
}










