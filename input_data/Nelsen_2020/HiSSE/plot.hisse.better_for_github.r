plot.hisse.states.hack<-function (x, rate.param, arc.marg, min.clade.size, do.observed.only = TRUE, rate.colors = NULL, 
    state.colors = NULL, edge.width.rate = 5, edge.width.state = 2, 
    type = "fan", rate.range = NULL, show.tip.label = TRUE, fsize = 1, 
    lims.percentage.correction = 0.001, legend = "tips", legend.position = c(0, 
        0.2, 0, 0.2), legend.cex = 0.4, legend.kernel.rates = "auto", 
    legend.kernel.states = "auto", legend.bg = "cornsilk3",rank="family", ...) 
{
	p2p<-diversitree:::plot2.phylo
	t <- max(branching.times(x[[1]]$phy))
	obj<-p2p(x[[1]]$phy,type="fan",label.offset = t * 1/4, show.tip.label = FALSE)
    hisse.results <- x
    if (class(hisse.results) == "hisse.states") {
        if (is.null(hisse.results$aic)) {
            hisse.results$aic = 42
        }
        tmp.list <- list()
        tmp.list[[1]] <- hisse.results
        hisse.results <- tmp.list
    }
    par(fig = c(0, 1, 0, 1), new = FALSE)
    if (is.null(rate.colors)) {
        rate.colors <- c("blue", "red")
    }
    if (is.null(state.colors)) {
        state.colors <- c("white", "black")
    }
    rates.tips <- hisse:::ConvertManyToRate(hisse.results, rate.param, 
        "tip.mat")
    rates.internal <- hisse:::ConvertManyToRate(hisse.results, rate.param, 
        "node.mat")
    states.tips <- NA
    states.internal <- NA
    if (do.observed.only) {
        states.tips <- hisse:::ConvertManyToBinaryState(hisse.results, 
            "tip.mat")
        states.internal <- hisse:::ConvertManyToBinaryState(hisse.results, 
            "node.mat")
    }
    else {
        stop("So far we can easily plot just the binary observed state; if you want to plot the hidden states, use a different function")
    }
    tree.to.plot <- hisse.results[[1]]$phy
    if (!show.tip.label) {
        rep.rev <- function(x, y) {
            result <- paste(rep(y, x), collapse = "", sep = "")
            return(result)
        }
        tree.to.plot$tip.label <- sapply(sequence(length(tree.to.plot$tip.label)), 
            rep.rev, " ")
        fsize = 0
    }
    rate.lims <- range(c(rates.tips, rates.internal))
    if (!is.null(rate.range)) {
        if (min(rate.range) > min(rate.lims) | max(rate.range) < 
            max(rate.lims)) {
            warning(paste("Did not override rate.lims: the specified rate.range (", 
                rate.range[1], ", ", rate.range[2], ") did not contain all the values for the observed rates (", 
                rate.lims[1], ", ", rate.lims[2], ")"))
        }
        else {
            rate.lims <- rate.range
        }
    }
    rate.lims[1] <- rate.lims[1] - lims.percentage.correction * 
        abs(rate.lims[1])
    rate.lims[2] <- rate.lims[2] + lims.percentage.correction * 
        abs(rate.lims[2])
    rate.tree <- hisse:::contMapGivenAnc(tree = tree.to.plot, x = rates.tips, 
        plot = FALSE, anc.states = rates.internal, lims = rate.lims, 
        ...)
    rate.colors <- colorRampPalette(rate.colors, space = "Lab")(length(rate.tree$cols))
    rate.tree$cols[] <- rate.colors
    plot(rate.tree, outline = FALSE, lwd = edge.width.rate, legend = FALSE, 
        type = type, fsize = fsize, ...)
	require(circlize)
	#THIS NUMBER NEEDS TO BE ADJUSTED
	#Caliciales:
	#circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-2.71,2.71),canvas.ylim=c(-2.71,2.71))
	#Peltigerales:
	#circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1.74,1.74),canvas.ylim=c(-1.74,1.74))
	circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1*arc.marg,arc.marg),canvas.ylim=c(-1*arc.marg,arc.marg))
	circos.initialize("a", xlim = c(0, 1),sector.width=1) # 'a` just means there is one sector 

	timescale<-read.csv(file="/Volumes/welwitschia/li'l rascal iv/Users/matthewnelsen/Documents/projects/carboniferous_buildup/carboniferous_revisions/25nov2015_data_pull_sap/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

	for(t in 1:nrow(timescale)){
		timescale$RGB[t]<-rgb(timescale$Col_R[t]/255,timescale$Col_G[t]/255,timescale$Col_B[t]/255,alpha=0.2)
	}
	require(phytools)
	require(phytools)
	ht<-max(nodeHeights(x[[1]]$phy))
	timey<-timescale[timescale$Type %in% "Period",]
	timey<-timey[timey$End<ht,]

	for(q in 1:nrow(timey)){
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-timey$End[q])/ht,rou2=(ht-timey$Start[q])/ht,col=timey$RGB[q],border=NA)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-50)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-100)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-150)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-200)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-250)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-300)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-350)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-400)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
		draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-450)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
	}
	for(j in 3:(nrow(timey)-1)){
		#print(x)
		text(x=(ht-rev(timey$Midpoint[j]))/ht,y=-0.03,labels=timey$Abbrev[j],cex=0.75,col="gray50")
	}
	circos.clear()
    par(fig = c(0, 1, 0, 1), new = TRUE)
    plot(rate.tree, outline = FALSE, lwd = edge.width.rate, legend = FALSE, 
        type = type, fsize = fsize, ...)
    par(fig = c(0, 1, 0, 1), new = TRUE)
    state.lims <- range(c(states.tips, states.internal))
    state.lims[1] <- state.lims[1] - lims.percentage.correction * 
        abs(state.lims[1])
    state.lims[2] <- state.lims[2] + lims.percentage.correction * 
        abs(state.lims[2])
    state.tree <- hisse:::contMapGivenAnc(tree = tree.to.plot, x = states.tips, 
        plot = FALSE, anc.states = states.internal, lims = state.lims, 
        ...)
    state.colors <- colorRampPalette(state.colors, space = "Lab")(length(rate.tree$cols))
    state.tree$cols[] <- state.colors
    plot(state.tree, outline = FALSE, lwd = edge.width.state, 
        legend = FALSE, type = type, fsize = fsize, ...)
	w=1/50
	div.arcs<-diversitree:::arcs
	col.bar<-"black"
	n <- obj$n.spp
	offset.bar=w * (n + 3)
	offset.lab=w * (n + 4)
	dy <-1/6
	dy <- dy/n * 2 * pi
	yy <- obj$xy$theta[seq_len(obj$Ntip)]
	#y0 <- tapply(yy - dy, tax2$Tribe, min)
	#y1 <- tapply(yy + dy, tax2$Tribe, max)
	#print("above tapply")
	y0 <- tapply(yy - dy, dat[,rank], min)
	y1 <- tapply(yy + dy, dat[,rank], max)
	x.bar <- rep(max(obj$xx) + offset.bar, length(y0))
	x.lab <- rep(max(obj$xx) + offset.lab, length(y0))
	#print("above keepers")
	keepers<-names(table(dat[,rank]))[table(dat[,rank])>=min.clade.size]
	#print("pas keepers")
	#print("above keepers")
	y0<-y0[names(y0) %in% keepers]
	y1<-y1[names(y1) %in% keepers]
	x.bar<-x.bar[1:length(keepers)]
	x.lab<-x.lab[1:length(keepers)]
	#print("above div.arcs")
	div.arcs(y0, y1, x.bar-max(x.bar)+ht+5, col = col.bar, lwd = 2)
	div.rt<-diversitree:::radial.text
	font=3
	cex=1
	ym <- (y0 + y1)/2
	str <- names(y0)
	str[!str %in% keepers]<-NA
	str<-gsub("_"," ",str)
	col.lab="black"
	div.rt(x.lab-x.lab+ht+10, ym, str, col = col.lab, font = font,cex = cex) 
    if (legend != "none") {
        par(fig = legend.position, new = TRUE)
        plot(x = c(-0.1, 1.1), y = c(-1.5, 1.5), xlab = "", ylab = "", 
            bty = "n", type = "n", xaxt = "n", yaxt = "n")
        rect(-0.1, -1.1, 1.1, 1.1, border = NA, col = legend.bg)
        rates.to.plot <- c()
        states.to.plot <- c()
        if (legend == "all" | legend == "tips") {
            rates.to.plot <- append(rates.to.plot, rates.tips)
            states.to.plot <- append(states.to.plot, states.tips)
        }
        if (legend == "all" | legend == "internal") {
            rates.to.plot <- append(rates.to.plot, rates.internal)
            states.to.plot <- append(states.to.plot, states.internal)
        }
        if (legend.kernel.rates == "auto") {
            if (length(unique(rates.to.plot)) <= 4) {
                legend.kernel.rates <- "hist"
            }
            else {
                legend.kernel.rates <- "rectangular"
            }
        }
        if (legend.kernel.states == "auto") {
            if (length(unique(states.to.plot)) <= 4) {
                legend.kernel.states <- "hist"
            }
            else {
                legend.kernel.states <- "rectangular"
            }
        }
        rates.density <- hisse:::GetNormalizedDensityPlot(rates.to.plot, 
            rate.lims, legend.kernel.rates)
        states.density <- hisse:::GetNormalizedDensityPlot(states.to.plot, 
            state.lims, legend.kernel.states)
        states.density$y <- (-1) * states.density$y
        par(lend = 1)
        segments(x0 = rates.density$x, y0 = rep(0, length(rates.density$y)), 
            y1 = rates.density$y, col = rate.colors[1 + as.integer(round((length(rate.colors) - 
                1) * rates.density$x))], lwd = ifelse(legend.kernel.rates == 
                "hist", 4, 1))
        text(x = 0, y = 1.2, labels = format(rate.lims[1], digits = 2), 
            cex = legend.cex)
        text(x = 1, y = 1.2, labels = format(rate.lims[2], digits = 2), 
            cex = legend.cex)
        text(x = 0.5, y = 1.2, labels = rate.param, cex = legend.cex)
        segments(x0 = states.density$x, y0 = rep(0, length(states.density$y)), 
            y1 = states.density$y, col = state.colors[1 + as.integer(round((length(state.colors) - 
                1) * states.density$x))], lwd = ifelse(legend.kernel.states == 
                "hist", 4, 1))
        text(x = 0, y = -1.2, labels = "0", cex = legend.cex)
        text(x = 1, y = -1.2, labels = "1", cex = legend.cex)
        text(x = 0.5, y = -1.2, labels = "State", cex = legend.cex)
    }
    return(list(rate.tree = rate.tree, state.tree = state.tree))
	circos.clear()
}


#files = system("ls -1 | grep .Rsave", intern=TRUE)
my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
#files<-list.files(path=my.path,pattern="Lecanoromycetidae_macro.bounded.ex.model.*.recon.Rsave")
#files<-list.files(path=my.path,pattern="teloschistales_macro.bounded.ex.model.*.recon.Rsave")
#files<-list.files(path=my.path,pattern="lich.bounded.ex.model.*.recon.Rsave")
#files<-list.files(path=my.path,pattern="primtreb.bounded.ex.model.*.recon.Rsave")
#arc.marg<-1.74
min.clade.size=3
#files<-list.files(path=my.path,pattern="peltigerales_macro.bounded.ex.model.*.recon.Rsave")
#arc.marg<-1.74
#files<-list.files(path=my.path,pattern="peltigerales_ceph.bounded.ex.model.*.recon.Rsave")
#arc.marg<-2.71
#files<-list.files(path=my.path,pattern="caliciales_macro.bounded.ex.model.*.recon.Rsave")
#arc.marg<-2.8
#files<-list.files(path=my.path,pattern="teloschistales_macro.bounded.ex.model.*.recon.Rsave")
arc.marg<-2.8
files<-list.files(path=my.path,pattern="cladoniineae_macro.bounded.ex.model.*.recon.Rsave")


my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.74
min.clade.size=3
files<-list.files(path=my.path,pattern="peltigerales_macro.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
#dev.off()



my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.74
min.clade.size=3
files<-list.files(path=my.path,pattern="peltigerales_ceph.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-280,280),ylim=c(-280,280),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)





my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.75
min.clade.size=3
files<-list.files(path=my.path,pattern="caliciales_macro.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)









my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.82
min.clade.size=3
files<-list.files(path=my.path,pattern="teloschistales_macro.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-180,180),ylim=c(-180,180),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)










my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.7
min.clade.size=3
files<-list.files(path=my.path,pattern="cladoniineae_macro.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)


my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.7
min.clade.size=3
files<-list.files(path=my.path,pattern="cladoniineae_ceph.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-100,100),ylim=c(-100,100),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)







my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.7
min.clade.size=10
files<-list.files(path=my.path,pattern="lecanorales_macro.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-200,200),ylim=c(-200,200),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-200,200),ylim=c(-200,200),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-200,200),ylim=c(-200,200),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-200,200),ylim=c(-200,200),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-200,200),ylim=c(-200,200),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
#dev.off()






my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.7
min.clade.size=10
files<-list.files(path=my.path,pattern="Lecanoromycetidae_macro.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-350,350),ylim=c(-350,350),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-350,350),ylim=c(-350,350),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-350,350),ylim=c(-350,350),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-350,350),ylim=c(-350,350),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-350,350),ylim=c(-350,350),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=3,edge.width.state=0.5)
#dev.off()














my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.65
min.clade.size=20
files<-list.files(path=my.path,pattern="^macro.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
#dev.off()


my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.65
min.clade.size=20
files<-list.files(path=my.path,pattern="^lich.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
#dev.off()


my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.65
min.clade.size=20
files<-list.files(path=my.path,pattern="^primcyan.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
#dev.off()


my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.65
min.clade.size=20
files<-list.files(path=my.path,pattern="^primtreb.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = pp
  rm(pp)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
#dev.off()









my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.65
min.clade.size=20
files<-list.files(path=my.path,pattern="^primtrent.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = hisse.fit.recon
  rm(hisse.fit.recon)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
#dev.off()


my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.65
min.clade.size=20
files<-list.files(path=my.path,pattern="^ceph.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = hisse.fit.recon
  rm(hisse.fit.recon)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)

#to do: trent, anycyano, ceph...doing reconstruction on vortex now (29 jan 2019)

my.path<-"~/Desktop/other/lecanoromycetes_diversification_2017_analyses/hisse_analyses/lichen_redo_hisse_analyses_w_fracs/hisse_analyses_bounded/"
arc.marg<-1.65
min.clade.size=20
files<-list.files(path=my.path,pattern="^anycyano.bounded.ex.model.*.recon.Rsave")
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(paste(my.path,files[i],sep=""))
  hisse.results.list[[i]] = hisse.fit.recon
  rm(hisse.fit.recon)
}
dat<-read.csv(file="~/Desktop/other/lecanoromycetes_diversification_2017_analyses/lecanoromycetes_correlated_evolution/ALL_w_accs_UPDATED_next3d_updatedforcladeanalysis_nocatmod_updforconstr.csv",stringsAsFactors=FALSE)
dat$edited_name<-gsub("-","_",dat$edited_name)
rownames(dat)<-dat$edited_name
dat<-dat[dat$edited_name %in% hisse.results.list[[1]]$phy$tip.label,]
dat<-dat[match(hisse.results.list[[1]]$phy$tip.label,rownames(dat)),]
#pdf(file=paste(my.path,"peltigerales_macro_speciation_hisse.pdf"))
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="speciation",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="net.div",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="turnover",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
plot.hisse.states.hack(hisse.results.list,arc.marg=arc.marg,rate.param="extinction.fraction",show.tip.label=FALSE,rate.colors=c("plum4","orange","yellow"),xlim=c(-390,390),ylim=c(-390,390),rank="family",min.clade.size=min.clade.size,legend.cex=1,edge.width.rate=1.5,edge.width.state=0.2)
