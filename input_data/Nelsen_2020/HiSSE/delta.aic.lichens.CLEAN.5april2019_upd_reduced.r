
#THESE IN THE delta.aic.lichens.CLEAN.30jan2019 were problematic: "macro_hisse_fractions_macro_1_micro_0.5","macro_hisse_fractions_macro_0.9_micro_0.4","macro_hisse_fractions_macro_0.7_micro_0.25","macro_hisse_fractions_macro_0.5_micro_0.1","macro_hisse_fractions_macro_0.9_micro_1","macro_hisse_fractions_macro_0.8_micro_1","macro_hisse_fractions_macro_0.7_micro_1","macro_hisse_fractions_macro_0.6_micro_1","macro_hisse_fractions_macro_0.5_micro_1","Lecanoromycetidae_macro"

#they are in same order as first batch like macro, primtreb, etc. dat$ModelNo<-c(10:19,1,20:24,2:9)
mods<-c("substrate")
path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
for(x in 1:length(mods)){
	dat<-read.delim(file=paste(path,mods[x],"_empirical.results.bounded.redo.txt",sep=""),stringsAsFactors=FALSE,header=FALSE)
	colnames(dat)<-c("loglik","AIC","index.par","solution")
	if(mods[x] %in% c("substrate")){
		#"macro_hisse_fractions_macro_0.5_micro_0.1" has problems with likelihood values - looks like a boundary hit
		dat$ModelNo<-c(10:19,1,20:24,2:9)
	}	
	dat$deltaAIC<-0
	dat$RelLik<-0
	dat$AICweight<-0
	best<-min(dat$AIC)
	for(p in 1:nrow(dat)){
		dat$deltaAIC[p]<-dat$AIC[p]-best
		dat$RelLik[p]<-exp(-0.5*dat$deltaAIC[p])
	}
	sum.rel.lik<-sum(dat$RelLik)
	for(p in 1:nrow(dat)){
		dat$AICweight[p]<-exp(-0.5*dat$deltaAIC[p])/sum.rel.lik
	}	
	dat<-dat[order(dat$ModelNo,decreasing=FALSE),]
	write.csv(dat,file=paste(path,mods[x],"_empirical.results.w.wts.MODS.CLEAN.5april2019.csv",sep=""))
}


df<-as.data.frame(matrix(nrow=24,ncol=(length(mods)*2)))
colnames(df)[seq(from=1,to=length(colnames(df)),by=2)]<-paste(mods,"AIC",sep=".")
colnames(df)[seq(from=2,to=length(colnames(df)),by=2)]<-paste(mods,"AICw",sep=".")
rownames(df)<-c("BiSSE: all free","BiSSE: extfrac0=extfrac1","BiSSE: q's equal","BiSSE: extfrac0=extfrac1, q's equal","CID-2: q's equal","CID-2: extfrac's, q's equal","CID-4: q's equal","CID-4: extfrac's equal, q's equal","HiSSE: q's equal","HiSSE: extfr's equal, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q's equal","HiSSE: q0B1B=0, q1B0B=0, all other q's equal","HiSSE: extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q0B1B=0,q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal")

for(x in 1:length(mods)){
	dat<-read.csv(file=paste(path,mods[x],"_empirical.results.w.wts.MODS.CLEAN.5april2019.csv",sep=""),stringsAsFactors=FALSE)
	#dat<-dat[order(dat$ModelNo),]
	colnames(dat)[5:12]<-c("turnover_0a","turnover_1a","turnover_0b","turnover_1b","extfrac_0a","extfrac_1a","extfrac_0b","extfrac_1b")
	for(z in 1:nrow(dat)){
		dat$a0_sp[z]<-dat$turnover_0a[z]/(1+dat$extfrac_0a[z])
		dat$a1_sp[z]<-dat$turnover_1a[z]/(1+dat$extfrac_1a[z])
		dat$b0_sp[z]<-dat$turnover_0b[z]/(1+dat$extfrac_0b[z])
		dat$b1_sp[z]<-dat$turnover_1b[z]/(1+dat$extfrac_1b[z])

		dat$a0_ex[z]<-(dat$turnover_0a[z]*dat$extfrac_0a[z])/(1+dat$extfrac_0a[z])
		dat$a1_ex[z]<-(dat$turnover_1a[z]*dat$extfrac_1a[z])/(1+dat$extfrac_1a[z])
		dat$b0_ex[z]<-(dat$turnover_0b[z]*dat$extfrac_0b[z])/(1+dat$extfrac_0b[z])
		dat$b1_ex[z]<-(dat$turnover_1b[z]*dat$extfrac_1b[z])/(1+dat$extfrac_1b[z])

		dat$a0_net[z]<-dat$a0_sp[z]-dat$a0_ex[z]
		dat$a1_net[z]<-dat$a1_sp[z]-dat$a1_ex[z]
		dat$b0_net[z]<-dat$b0_sp[z]-dat$b0_ex[z]
		dat$b1_net[z]<-dat$b1_sp[z]-dat$b1_ex[z]

		df[dat$ModelNo[z],((x*2)-1)]<-dat$AIC[z]
		df[dat$ModelNo[z],(x*2)]<-round(dat$AICweight[z],3)
		if(df[dat$ModelNo[z],(x*2)]<0.001){
			df[dat$ModelNo[z],(x*2)]<-"<0.001"
		}
	colnames(dat)[1]<-"Model"
	dat$Model<-c("BiSSE: all free","BiSSE: extfrac0=extfrac1","BiSSE: q's equal","BiSSE: extfrac0=extfrac1, q's equal","CID-2: q's equal","CID-2: extfrac's, q's equal","CID-4: q's equal","CID-4: extfrac's equal, q's equal","HiSSE: q's equal","HiSSE: extfr's equal, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q's equal","HiSSE: q0B1B=0, q1B0B=0, all other q's equal","HiSSE: extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q0B1B=0,q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal")
	write.csv(dat,file=paste(path,mods[x],"_empirical.results.w.wts.MODS.CLEAN.5april2019.csv",sep=""))
	}
}



#make summary file:
mods<-c("macro","primtreb","primtrent","primcyan","ceph","anycyano","lich","peltigerales_macro","peltigerales_ceph","caliciales_macro","cladoniineae_macro","cladoniineae_trip","lecanorales_macro","teloschistales_macro","macro_hisse_fractions_macro_1_micro_0.5","macro_hisse_fractions_macro_0.9_micro_0.4","macro_hisse_fractions_macro_0.7_micro_0.25","macro_hisse_fractions_macro_0.5_micro_0.1","macro_hisse_fractions_macro_0.9_micro_1","macro_hisse_fractions_macro_0.8_micro_1","macro_hisse_fractions_macro_0.7_micro_1","macro_hisse_fractions_macro_0.6_micro_1","macro_hisse_fractions_macro_0.5_micro_1","Lecanoromycetidae_macro")
path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
df<-as.data.frame(matrix(nrow=24,ncol=(length(mods)*2)))
colnames(df)[seq(from=1,to=length(colnames(df)),by=2)]<-paste(mods,"AIC",sep=".")
colnames(df)[seq(from=2,to=length(colnames(df)),by=2)]<-paste(mods,"AICw",sep=".")
rownames(df)<-c("BiSSE: all free","BiSSE: extfrac0=extfrac1","BiSSE: q's equal","BiSSE: extfrac0=extfrac1, q's equal","CID-2: q's equal","CID-2: extfrac's, q's equal","CID-4: q's equal","CID-4: extfrac's equal, q's equal","HiSSE: q's equal","HiSSE: extfr's equal, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q's equal","HiSSE: q0B1B=0, q1B0B=0, all other q's equal","HiSSE: extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q0B1B=0,q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal")

for(x in 1:length(mods)){
	dat<-read.csv(file=paste(path,mods[x],"_empirical.results.w.wts.MODS.CLEAN.30jan2019.csv",sep=""),stringsAsFactors=FALSE)
	for(z in 1:nrow(dat)){
		df[dat$ModelNo[z],((x*2)-1)]<-dat$AIC[z]
		df[dat$ModelNo[z],(x*2)]<-round(dat$AICweight[z],3)
		if(df[dat$ModelNo[z],(x*2)]<0.001){
			df[dat$ModelNo[z],(x*2)]<-"<0.001"
		}
	}
}

write.csv(df,file=paste(path,"summary.MODS.CLEAN.5april2019.csv",sep=""))
#and make summary for ALL modelfittings with speciation/extinction?


#NOTE that 30jan2019 read in for files that were corrected above
#these ended up in Tables S16 (Lecanoromycetidae macro) and S23 (undersampled diversity)
#CID 4 is still best for AIC macro in S16 and was manually pasted in corrected SI.
#Similarly, CID4 is still best model for S23

#make summary file:
mods<-c("macro_hisse_fractions_macro_1_micro_0.5","macro_hisse_fractions_macro_0.9_micro_0.4","macro_hisse_fractions_macro_0.7_micro_0.25","macro_hisse_fractions_macro_0.5_micro_0.1","macro_hisse_fractions_macro_0.9_micro_1","macro_hisse_fractions_macro_0.8_micro_1","macro_hisse_fractions_macro_0.7_micro_1","macro_hisse_fractions_macro_0.6_micro_1","macro_hisse_fractions_macro_0.5_micro_1")
path<-"~/Desktop/hisse_analyses_bounded/"
df<-as.data.frame(matrix(nrow=24,ncol=(length(mods)*2)))
colnames(df)[seq(from=1,to=length(colnames(df)),by=2)]<-paste(mods,"AIC",sep=".")
colnames(df)[seq(from=2,to=length(colnames(df)),by=2)]<-paste(mods,"AICw",sep=".")
rownames(df)<-c("BiSSE: all free","BiSSE: extfrac0=extfrac1","BiSSE: q's equal","BiSSE: extfrac0=extfrac1, q's equal","CID-2: q's equal","CID-2: extfrac's, q's equal","CID-4: q's equal","CID-4: extfrac's equal, q's equal","HiSSE: q's equal","HiSSE: extfr's equal, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q's equal","HiSSE: q0B1B=0, q1B0B=0, all other q's equal","HiSSE: extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac0A=extfrac1A=extfrac0B, q0B1B=0,q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q's equal","HiSSE: netturn0A=netturn0B, extfrac0A=extfrac0B, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn0B, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q's equal","HiSSE: netturn0A=netturn1A, extfrac0A=extfrac1A, q0B1B=0, q1B0B=0, all other q's equal","HiSSE: netturn0A=netturn1A, extfrac's equal, q0B1B=0, q1B0B=0, all other q's equal")

for(x in 1:length(mods)){
	dat<-read.csv(file=paste(path,mods[x],"_empirical.results.w.wts.MODS.CLEAN.5april2019.csv",sep=""),stringsAsFactors=FALSE)
	for(z in 1:nrow(dat)){
		df[dat$ModelNo[z],((x*2)-1)]<-dat$AIC[z]
		df[dat$ModelNo[z],(x*2)]<-round(dat$AICweight[z],3)
		if(df[dat$ModelNo[z],(x*2)]<0.001){
			df[dat$ModelNo[z],(x*2)]<-"<0.001"
		}
	}
}

write.csv(df,file=paste(path,"summary.MODS.CLEAN.S23.csv",sep=""))
#and make summary for ALL modelfittings with speciation/extinction?


