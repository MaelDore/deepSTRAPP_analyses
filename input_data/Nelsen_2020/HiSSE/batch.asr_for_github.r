 ########Start of the model average tip and node algorithm#########
GetBatchAsr <- function(file.name, n.cores){
	models <- system(paste("ls -1 ", file.name, "*.Rsave", sep=""), intern=TRUE)
	for(i in 1:length(models)){
		load(models[i])
		if(is.null(hisse.fit$hidden.states)){
			pp <- MarginRecon(hisse.fit$phy, hisse.fit$data, f=hisse.fit$f, pars=hisse.fit$solution, four.state.null=TRUE, aic=hisse.fit$AIC, n.cores=n.cores)
			save(pp, file=paste(file.name, i, "recon", "Rsave", sep="."))
		}else{
			pp <- MarginRecon(hisse.fit$phy, hisse.fit$data, f=hisse.fit$f, pars=hisse.fit$solution, hidden.states=hisse.fit$hidden.states, aic=hisse.fit$AIC, n.cores=n.cores)
			save(pp, file=paste(file.name, i, "recon", "Rsave", sep="."))
		}
	}
}
#GetBatchAsr("achene.ex.model", n.cores=24)


my.analyses<-c("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/cladoniineae_ceph.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/cladoniineae_macro.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/caliciales_macro.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/lecanorales_macro.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/peltigerales_ceph.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/peltigerales_macro.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/Lecanoromycetidae_macro.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/macro.bounded.ex.model","/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/primcyan.bounded.ex.model"
for(p in 1:length(my.analyses)){
	GetBatchAsr(my.analyses[p],n.cores=12)
}