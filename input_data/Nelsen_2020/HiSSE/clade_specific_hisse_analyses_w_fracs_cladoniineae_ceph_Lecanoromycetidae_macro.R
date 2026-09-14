
require(ape)
require(phytools)
require(hisse)
require(parallel)


	
SimFunction <- function(phy=NULL,sim.dat=NULL,sampling.file=NULL,file.name=NULL){
	sampling.file[,2]<-sampling.file[,1]
	sampling.file.new <- sampling.file[phy$tip.label,]
	sampling.matrix <- as.matrix(sampling.file.new[,1])
	sampling.f <- as.vector(sampling.matrix)
	trans.rates.hisse <- TransMatMaker(hidden.states=TRUE)
	trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
	trans.rates.hisse.test <- TransMatMaker(hidden.states=TRUE)
	trans.rates.hisse.test <- ParDrop(trans.rates.hisse.test, c(3,5,8,9,10,12))
	trans.rates.hisse.red <- trans.rates.hisse
	trans.rates.hisse.red.test <- trans.rates.hisse.test
	trans.rates.hisse.red[!is.na(trans.rates.hisse.red) & !trans.rates.hisse.red == 0] = 1
	trans.rates.hisse.red.test[!is.na(trans.rates.hisse.red.test) & !trans.rates.hisse.red.test == 0] = 1
#	trans.rates.hisse.red.test[2,1]=2
#	trans.rates.hisse.red.test[1,2]=3
	trans.rates.bisse <- TransMatMaker(hidden.states=FALSE)
	trans.rates.bisse.red <- trans.rates.bisse
	trans.rates.bisse.red[!is.na(trans.rates.bisse.red)] = 1
	
	hisse.fit <- NA
	bisse.fit <- NA
	
	RunModel <- function(model.number){
		if(model.number==1){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse))	
		}
		if(model.number==2){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse))
		}
		if(model.number==3){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,2,0,0), trans.rate=trans.rates.bisse.red))	
		}
		if(model.number==4){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=FALSE, turnover.anc=c(1,2,0,0), eps.anc=c(1,1,0,0), trans.rate=trans.rates.bisse.red))
		}
		if(model.number==5){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,2,2), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==6){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==7){
			try(hisse.fit <- hisse.null4(phy, sim.dat, f=sampling.f, bounded.search=TRUE))
		}
		if(model.number==8){
			try(hisse.fit <- hisse.null4(phy, sim.dat, f=sampling.f, bounded.search=TRUE, eps.anc=rep(1,8)))
		}
		if(model.number==9){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==10){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==11){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==12){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))		
		}
		if(model.number==13){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,2,3,4), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==14){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==15){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,2), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==16){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,1,2), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))	
		}
		if(model.number==17){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==18){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))	
		}
		if(model.number==19){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,2,1,3), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==20){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,2,1,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))	
		}
		if(model.number==21){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red))
		}
		if(model.number==22){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red))	
		}
		if(model.number==23){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,2,3), trans.rate=trans.rates.hisse.red.test))
		}
		if(model.number==24){
			try(hisse.fit <- hisse(phy, sim.dat, f=sampling.f, bounded.search=TRUE, hidden.states=TRUE, turnover.anc=c(1,1,2,3), eps.anc=c(1,1,1,1), trans.rate=trans.rates.hisse.red.test))	
		}		
		save(phy, hisse.fit, file=paste(file.name, model.number, "bounded.Rsave", sep="."))
	}
	mclapply(1:24, RunModel, mc.cores=12)
}


#dat.nd<-read.csv(file="~/price_et_al_analyses/Dataset.modified.csv",stringsAsFactors=FALSE)
#head(dat.nd)
#tr.dr[[1]]
dataset<-"Lecanoromycetidae_macro"
filepath<-"/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/"
tax<-read.csv(file="ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_for_bamm_subclass.csv",header=TRUE,stringsAsFactors=FALSE)
tax$edited_name<-gsub("-","_",tax$edited_name)


for(pq in 1:length(dataset)){
	if(dataset[pq] %in% "Lecanoromycetidae_macro"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_macro_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$subclass %in% "Lecanoromycetidae"],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$subclass %in% "Lecanoromycetidae"]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$subclass %in% "Lecanoromycetidae"])),FALSE)
		write.tree(tr.dr,file="Lecanoromycetidae_macro_hisse.tre")
	}
	if(dataset[pq] %in% "peltigerales_macro"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_macro_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Peltigerales"],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Peltigerales"]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$order %in% "Peltigerales"])),FALSE)
		write.tree(tr.dr,file="peltigerales_macro_hisse.tre")
	}
	if(dataset[pq] %in% "peltigerales_ceph"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_ceph_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Peltigerales"],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Peltigerales"]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$order %in% "Peltigerales"])),FALSE)
		write.tree(tr.dr,file="peltigerales_ceph_hisse.tre")
	}
	if(dataset[pq] %in% "caliciales_macro"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_macro_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Caliciales"],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Caliciales"]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$order %in% "Caliciales"])),FALSE)
		write.tree(tr.dr,file="caliciales_macro_hisse.tre")
	}
	if(dataset[pq] %in% "cladoniineae_macro"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_macro_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")])),FALSE)
		write.tree(tr.dr,file="cladoniineae_macro_hisse.tre")
	}
	if(dataset[pq] %in% "cladoniineae_ceph"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_ceph_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")])),FALSE)
		write.tree(tr.dr,file="cladoniineae_ceph_hisse.tre")
	}
	if(dataset[pq] %in% "lecanorales_macro"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_macro_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Lecanorales"],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Lecanorales"]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$order %in% "Lecanorales"])),FALSE)
		write.tree(tr.dr,file="lecanorales_macro_hisse.tre")
	}
	if(dataset[pq] %in% "teloschistales_macro"){
		df <-read.csv(file=paste(filepath,"ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr_macro_forcorhmm.csv",sep=""),header=TRUE,stringsAsFactors=FALSE)
		sampling.fractions<-read.delim("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_fractions.txt")
		sampling.fractions<-data.frame(sampling.fractions[,2],row.names=sampling.fractions[,1],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		sampling.fractions<-data.frame(sampling.fractions[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Teloschistales"],],row.names=rownames(sampling.fractions)[rownames(sampling.fractions) %in% tax$edited_name[tax$order %in% "Teloschistales"]],stringsAsFactors=FALSE)
		colnames(sampling.fractions)<-"f"
		tr.dr<-read.tree("~/lichen_redo_hisse_analyses_w_fracs/ExaML_result.subclass_family_constraint_ML_TREE.trimmed_treepl_med_families_cv_LADDERIZED.tre")
		tr.dr<-ladderize(drop.tip(tr.dr,c(tax$edited_name[!tax$order %in% "Teloschistales"])),FALSE)
		write.tree(tr.dr,file="teloschistales_macro_hisse.tre")
	}		
	#df<-data.frame(dat.nd[,c("Taxon",dataset[pq])],stringsAsFactors=FALSE)
	colnames(df)<-c("Genus_species","states")
	rownames(df)<-df$Genus_species
	SimFunction(phy=tr.dr,sim.dat=df,sampling.file=sampling.fractions,file.name=paste("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/",dataset[pq],".bounded.ex.model",sep=""))
	EmpiricalRunSummary <- function(file.name){
		models <- system(paste("ls -1 ", file.name, "*.bounded.Rsave", sep=""), intern=TRUE)
		for(i in 1:length(models)){
			load(models[i])
			write.table(t(c(hisse.fit$loglik, hisse.fit$AIC, max(hisse.fit$index.par)-1, hisse.fit$solution)), file=paste("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/",dataset[pq],"_empirical.results.bounded.txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
		}
	}
	EmpiricalRunSummary(paste("/home/mpnelsen/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/",dataset[pq],".bounded.ex.model",sep=""))
}

