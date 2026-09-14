#From Beaulieu & O'Meara SI

GetBatchAsr <- function(file.name, n.cores){
	models <- system(paste("ls -1 ", file.name, "*.Rsave", sep=""), intern=TRUE)
	for(i in 1:length(models)){
		load(models[i])
		if(is.null(hisse.fit$hidden.states)){
			hisse.fit.recon <- MarginRecon(hisse.fit$phy, hisse.fit$data, f=hisse.fit$f, pars=hisse.fit$solution, four.state.null=TRUE, aic=hisse.fit$AIC, n.cores=n.cores)
			save(hisse.fit.recon, file=paste(strsplit(models[i],"Rsave")[[1]],"recon.Rsave", sep=""))
		}else{
			hisse.fit.recon <- MarginRecon(hisse.fit$phy, hisse.fit$data, f=hisse.fit$f, pars=hisse.fit$solution, hidden.states=hisse.fit$hidden.states, aic=hisse.fit$AIC, n.cores=n.cores)
			save(hisse.fit.recon, file=paste(strsplit(models[i],"Rsave")[[1]],"recon.Rsave", sep=""))
		}
	}
}



#skip this, as all show CID or BiSSE
#dataset<-c("macro","primtreb","primtrent","primcyan","ceph","anycyano","lich")
#filepath<-"/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/"
#n.cores=8

#for(pq in 1:length(dataset)){
#	file.name=paste("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/",dataset[pq],".bounded.ex.model",sep="")
#	GetBatchAsr(file.name=file.name, n.cores=n.cores)
#}

#skip this - all CID
#dataset<-c("macro_hisse_fractions_macro_1_micro_0.5","macro_hisse_fractions_macro_0.9_micro_0.4","macro_hisse_fractions_macro_0.7_micro_0.25","macro_hisse_fractions_macro_0.5_micro_0.1","macro_hisse_fractions_macro_0.9_micro_1","macro_hisse_fractions_macro_0.8_micro_1","macro_hisse_fractions_macro_0.7_micro_1","macro_hisse_fractions_macro_0.6_micro_1","macro_hisse_fractions_macro_0.5_micro_1")
#filepath<-"/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/"
#n.cores=8

#for(pq in 1:length(dataset)){
#	file.name=paste("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/",dataset[pq],".bounded.ex.model",sep="")
#	GetBatchAsr(file.name=file.name, n.cores=n.cores)
#}


dataset<-c("peltigerales_macro","peltigerales_ceph","caliciales_macro","cladoniineae_macro","cladoniineae_trip","lecanorales_macro","teloschistales_macro")
filepath<-"/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/"
n.cores=8

for(pq in 1:length(dataset)){
	file.name=paste("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/",dataset[pq],".bounded.ex.model",sep="")
	GetBatchAsr(file.name=file.name, n.cores=n.cores)
}


require(hisse)
dataset<-c("primtrent","ceph","anycyano")
filepath<-"/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/"
n.cores=8
for(pq in 1:length(dataset)){
	file.name=paste("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/",dataset[pq],".bounded.ex.model",sep="")
	GetBatchAsr(file.name=file.name, n.cores=n.cores)
}
