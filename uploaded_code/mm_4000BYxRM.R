#4000 BYxRM data
load("/home/jbloom/Dropbox/4000BYxRM/uploaded_data/pheno_raw.RData")
#pheno_raw
load("/home/jbloom/Dropbox/4000BYxRM/uploaded_data/cross.RData")

#mcor=cor(cross$pheno, use='pairwise.complete.obs')
#png(file='~/Dropbox/4000BYxRM/trait_correlations.png', width=1080, height=1080)
#par(oma=c(2,2,2,2))
#corrplot(mcor, tl.pos='td', type='upper', order='hclust', method='shade', diag=F, addCoef.col='black')
#dev.off()

#x=read.xls('~/Dropbox/4000BYxRM/working_tables/1000v4000_BYxRM.xlsx')
#plot(x[,4], x$broad4000, xlab='1000', ylab='4000', main='repeatibility')
#plot(x[,6], x$add4000, xlab='1000', ylab='4000', main='additive')
#plot(x$broad4000-x$add4000, (x[,4]-x[,6]))
#cor.test(x[,6], x$add4000)
#cor.test(x[,4], x$broad4000)

#4000 BYxRM data
#load('~/Dropbox/cross.RData')
library(regress)
library(rrBLUP)
library(qtl)
library(fields)
library(gputools)
library(Matrix)
source('/home/jbloom/Dropbox/4000BYxRM/uploaded_code/calcMM.R')

#load tetrad.chart
load('/home/jbloom/Dropbox/4000BYxRM/uploaded_data/tetrad.chart')

# some useful functions
extractScaledPhenotype=function(impcross,scaleVar=FALSE){apply(impcross$pheno, 2, scale, scale=scaleVar)}
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$data }))*2)-3 }
extractGenotype.argmax=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$argmax }))*2)-3 }
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}
get.LOD.by.COR = function(n.pheno, pheno, gdata, doGPU=T) {
    if(doGPU) {
    # Lynch and Walsh p. 454
    return( (-n.pheno*log(1-gpuCor(pheno, gdata, use='pairwise.complete.obs')$coefficients^2))/(2*log(10)) )  }
    else {
    return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10)) )  }
}
# needs negative log likelihoods
calc.pval.LRtest=function(null,full) { pchisq(-2*(null-full),1, lower.tail=FALSE) }
calc.BLUPS= function(G,Z,Vinv,y,X,B ){    G%*%t(Z)%*%Vinv%*%(y- X%*%B)     }
extractVarCompResults = function(r) {list(sigma=r$sigma, sigma.cov=r$sigma.cov, llik=r$llik) }
extractVarCompResultsJB = function(r) {list(sigma=r$Var, sigma.cov=r$invI, Bhat=as.vector(r$Bhat), llik=as.vector(r$llik)) }

#extract genotype data
gdata     = extractGenotype(cross)

gcoord=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,1])
# chr name
chr=as.character(do.call('rbind', strsplit(colnames(gdata), '_'))[,2])

#for 4000 BYxRM-------------------------------------------------------------
chr = gsub("chr", "", chr)
chr.lab = match(chr, as.character(as.roman(1:16)))

pos=as.numeric(do.call('rbind', strsplit(colnames(gdata), '_'))[,3])

# fix this
all.strain.names=names(pheno_raw[[1]])
#good.strain.names=all.strain.names[-bad.segs]

#nsample  =nrow(gdata)
#nmarker  =ncol(gdata)
newMM.results=list()
newMM.cblups =list()
newMM.peaks =list()
newMM.blupResids = list()

for(phenotype in names(pheno_raw) ) {
    s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'

    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------

    y=ny[!is.na(ny)]

    strain.names=(names(y))

    ##### extra filter############################
    # test if 
    #bsfilt=strain.names %in% good.strain.names
    #y=y[bsfilt]
    #strain.names=strain.names[bsfilt]
    ##############################################
    
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 

    strains.with.phenos=match(unique.sn, all.strain.names)

    n.strains=length(strains.with.phenos)
   
    # calculate relatedness matrix from genomewide marker data
    # also, in this case, this is equivalent to (gdata %*% t(gdata))
    A  = A.mat(gdata[strains.with.phenos,], shrink=FALSE)/2
    
    # added tet effect
    Tet = tetrad.chart[strains.with.phenos,]
    #Tet[Tet==0]=-1
    Tet = Tet %*% t(Tet)
    Tet = Tet/ncol(tetrad.chart)
    
    AA = A*A
    AAA = A*A*A

    y.avg=as.vector(by(y, names(y), mean, na.rm=T))
    # Calculate LODS, avg data
    LODS.avg=get.LOD.by.COR(n.strains, 
                            y.avg, gdata[strains.with.phenos,])
    options(warn = -1)

    mm.broad = calcMM(y, alg='fs',reps=TRUE)
    mm.broad = extractVarCompResultsJB(mm.broad)

    mm.A = calcMM(y, B=list(A=A), alg='fs',reps=TRUE)
    newMM.blupResids[[phenotype]] = calc.BLUPS(mm.A$Var[1]*A +  mm.A$Var[2]* Strain+mm.A$Var[3]*diag(n.strains), Z, mm.A$W,  y, rep(1,length(y)), mm.A$Bhat) - 
                 calc.BLUPS(mm.A$Var[1]*A , Z, mm.A$W,  y, rep(1,length(y)), mm.A$Bhat)  
    mm.A = extractVarCompResultsJB(mm.A)
   
    #LODS=get.LOD.by.COR(n.strains, blupResids, gdata[strains.with.phenos,])
    #bminusp= residuals(lm(as.vector(newMM.blupResids[[phenotype]])~gdata[strains.with.phenos,peaks]))
    #LODS=get.LOD.by.COR(n.strains, bminusp, gdata[strains.with.phenos,])

    mm.AA = calcMM(y, B=list(A=A,AA=AA), alg='fs',reps=TRUE)

    mm.tet = calcMM(y, B=list(A=A,AA=AA,TT=Tet), alg='ai',reps=TRUE)
    mm.tet = calcMM(y, B=list(TT=Tet), alg='ai',reps=TRUE)

    mm.AA = extractVarCompResultsJB(mm.AA)

    # equivalent, maybe even faster and more memory efficient
    print(phenotype)
    print(mm.AA)

    #list to store relatedness matrices built per chromosome
    #A.chr=list()
    # cov matrices for genome minus a chromosome
    A.chr.m = list()
    for(i in unique(chr.lab)){
        print(i)
        A.chr.m[[i]]=(A.mat(gdata[strains.with.phenos,chr.lab!=i], shrink=FALSE)/2)
    }
    names(A.chr.m)=paste('m_', names(A.chr.m), sep='')
    
    m.chr.blups.wo = list()
    mm.null.lliks = rep(NA,16)

    for ( i in 1:16) {
     print(i)
       mm.null = calcMM(y, B=list(A.m=A.chr.m[[i]], AA=AA), conv.val=1e-4, alg='ai',reps=T, Var=mm.AA$sigma)
       b = calc.BLUPS(mm.null$Var[1]*A.chr.m[[i]]+ mm.null$Var[2]*AA, Z, mm.null$W,  y, rep(1,length(y)),mm.null$Bhat)         
       y.all = calc.BLUPS(mm.null$Var[1]*A.chr.m[[i]]+ mm.null$Var[2]*AA + mm.null$Var[3]*Strain+mm.null$Var[4]*diag(n.strains), Z, mm.null$W,  y, rep(1,length(y)),mm.null$Bhat)         
       m.chr.blups.wo[[i]]=y.all-b
       mm.null.lliks[i]=as.vector(mm.null$llik)
    }

    c.blups=do.call('cbind', m.chr.blups.wo)
    c.blups = sapply(m.chr.blups.wo, as.vector)
    newMM.cblups[[phenotype]]=c.blups
    #mm.full.all.chr=regress(y~1, ~ZA[[1]]+ZA[[2]]+ZA[[3]]+ZA[[4]]+ZA[[5]]+ZA[[6]]+ZA[[7]]+ZA[[8]]+ZA[[9]]+ZA[[10]]+ZA[[1]]+ZA[[12]]+ZA[[13]]+ZA[[14]]+ZA[[15]]+ZA[[16]]+ZStraintZ, pos=rep(TRUE,18), verbose=9)
   
    LODS=get.LOD.by.COR(n.strains, c.blups, gdata[strains.with.phenos,])
    #plot(LODS[1,],ylim=c(0, max(LODS)) )
    #for(i in 2:16) { points(LODS[i,], col=i) }

    uc=unique(chr.lab)
    peaks=c()  

    #h2 nqtlh2 nqtlint epph2
    pdfout=paste('~/Dropbox/new_mm/4000BYxRM/', phenotype, '.pdf', sep='')
    pdf(file=pdfout, width=11, height=8)
    for ( i in 1:16){
        print(i)
        chr.ind=seq(1,ncol(gdata))[chr.lab==uc[i] ]
        g.s = gdata[strains.with.phenos,chr.ind]
        r.vec= c.blups[,i]
        plot.max=max(LODS[i,][which(chr.lab==uc[i])])
        plot(LODS[i,][which(chr.lab==uc[i])], ylim=c(0, plot.max+2), type='l',lwd=3, main=paste('Chr', i), ylab='LOD')
        mp  = which.max(LODS[i,][chr.ind])
        m.l = max(LODS[i,][chr.ind])
        perm.thresh=quantile(apply(get.LOD.by.COR(n.strains, replicate(1000,sample(r.vec)), g.s),1,max), .95)
        peaks.picked=0
        while(m.l>perm.thresh) {
            abline(v=mp, col=peaks.picked+1)
            peaks.picked=peaks.picked+1
            peaks=c(peaks, chr.ind[mp])
            r.vec=residuals(lm(r.vec~g.s[,mp]-1))
            lod.resid = as.vector(get.LOD.by.COR(n.strains, r.vec, g.s) )
            points(lod.resid, type='l', col=peaks.picked+1)
            mp  = which.max(lod.resid)
            m.l = max(lod.resid)
            perm.thresh=quantile(apply(get.LOD.by.COR(n.strains, replicate(1000,sample(r.vec)), g.s),1,max), .95)
        }
       # naive
       points(LODS.avg[1,][which(chr.lab==uc[i])], type='l',lwd=2, lty=3)
    }
    dev.off()
    newMM.peaks[[phenotype]]=peaks
    
    peak.A = A.mat(gdata[strains.with.phenos,peaks], shrink=FALSE)/2
    pAA    = peak.A*peak.A
    pAG    = peak.A *A 
    
    mm.peaks.AA = calcMM(y, B=list(A=peak.A,AA=pAA), alg='fs',reps=TRUE)
    mm.peaks.AA = extractVarCompResultsJB(mm.peaks.AA)

    mm.peaks.AAg = calcMM(y, B=list(A=peak.A,AA=pAG), alg='fs',reps=TRUE)
    mm.peaks.AAg = extractVarCompResultsJB(mm.peaks.AAg)

    
    newMM.results[[phenotype]]=list(mm.broad=mm.broad,
                                    mm.A=mm.A, 
                                    mm.AA=mm.AA, 
                                    mm.peaks.AA=mm.peaks.AA, 
                                    mm.peaks.AAg =mm.peaks.AAg)

    results=list(newMM.peaks, newMM.results, newMM.cblups, newMM.blupResids )
    #save(file='~/Dropbox/new_mm/4000BYxRM/Results.RData',results )
}
load('/home/jbloom/Dropbox/4000BYxRM/Results.RData')

newMM.peaks=results[[1]]
newMM.peaks=lapply(newMM.peaks, function(x) {unique(x)})
newMM.results=results[[2]]
newMM.cblups=results[[3]]
newMM.blupResids=results[[4]]

s.cross=sim.geno(cross)
peak.ldrop=list()
# add QTL peak detection 
for(i in 1:20) {
    print(i)
    aqmarkers=colnames(gdata)[newMM.peaks[[i]]]
    find.markerpos.obj=find.markerpos(cross, aqmarkers)
    makeqtl.obj=makeqtl(s.cross,find.markerpos.obj$chr, find.markerpos.obj$pos )
    rqtl.obj = refineqtl(s.cross, pheno.col=i, makeqtl.obj, method='imp', model='normal', verbose=T, maxit=1, maxit.fitqtl=1)
    peak.ldrop[[i]]=lapply(1:length(aqmarkers), function(x) lodint(rqtl.obj, qtl.index=x) )
}

peak.widths=lapply(peak.ldrop, function(x) {
      zname= sapply(x,  function (y) {
            mname.peak.ldrop =    as.numeric( do.call('rbind', strsplit(rownames(y), '_'))[,3])
             #return(mname.peak.ldrop)
            z=mname.peak.ldrop[c(1,length(mname.peak.ldrop))]
            return(z[2]-z[1])
})
      x2=x
      names(x2)=zname
      return(x2)})
peak.widths.vec=sapply(peak.widths, function(x) {as.numeric(names(x))})
names(pheno_raw)

pwsort=peak.widths
for(i in 1:20) {
pwsort[[i]]=peak.widths[[i]][as.character(sort(as.numeric(names(peak.widths[[i]])),decreasing=T))]
}
#names(pheno_raw)
# [1] "1_CobaltChloride_1"    "1_CopperSulfate_1"     "1_Diamide_1"           "1_E6-Berbamine_1"      "1_Ethanol_1"           "1_Formamide_1"         "1_Hydroxyurea_1"       "1_IndolaceticAcid_1"  
# [9] "1_Lactate_1"           "1_Lactose_1"           "1_MagnesiumChloride_1" "1_ManganeseSulfate_1"  "1_Menadione_1"         "1_Neomycin_1"          "1_Raffinose_1"         "1_Trehalose_1"        
#[17] "1_Xylose_1"            "1_YNB_1"               "1_YPD_1"               "1_Zeocin_1"           
#save(peak.widths, file='~/Dropbox/peak_widths.RData')

#vc3
load('/home/jbloom/Dropbox/4000BYxRM/uploaded_data/vc3.RData')

sigmas=lapply(newMM.results, function(x) lapply(x, function(y) y$sigma ))
sigma.cov=lapply(newMM.results, function(x) lapply(x, function(y) sqrt(diag(y$sigma.cov)) ))
sigmas3 = lapply(vc3, function(x) x$sigma)
sigmas3.cov = lapply(vc3, function(x) sqrt(diag(x$sigma.cov)))
for(i in 1:20) {
    sigmas[[i]]$mm.AAA = sigmas3[[i]]
    sigma.cov[[i]]$mm.AAA = sigmas3.cov[[i]]

}
nf=lapply(sigmas, function(x) {sapply(x, sum)})

for(i in 1:20){
    for(j in 1:6) {
    sigmas[[i]][[j]]=sigmas[[i]][[j]]/nf[[i]][j]
    sigma.cov[[i]][[j]]=sigma.cov[[i]][[j]]/nf[[i]][j]
    }
}

# Supp Table 1 
#avgMM.table
repMM.table=do.call('cbind', lapply(sigmas, function(x) {do.call('rbind', x) }))
colnames(repMM.table) = names(sigmas)
colnames(repMM.table)=gsub('1_', '', colnames(repMM.table))
colnames(repMM.table)=gsub('_1', '', colnames(repMM.table))
#"mm.broad"     "mm.A"         "mm.AA"        "mm.peaks.AA"  "mm.peaks.AAg" "mm.AAA"
write.table(repMM.table, file='~/Dropbox/4000BYxRM/repMM_table.txt', sep='\t', quote=F,row.names=F)
repMM.se.table=do.call('cbind', lapply(sigma.cov, function(x) do.call('c', x)))
write.table(repMM.se.table, file='~/Dropbox/4000BYxRM/repMM_table_se.txt', sep='\t', quote=F,row.names=F)



sortBP = function(vc, vc.se, ovc=NULL) {
    colnames(vc)=gsub('1_', '', colnames(vc))
    colnames(vc)=gsub('_1', '', colnames(vc))
    # vc = bpAA
   # vc.se = bpAA.se
    vc.cum  =apply(vc, 2, cumsum)
    if(is.null(ovc)) {    ovc = order(vc.cum[3,]) }
    vc=vc[,ovc]
    vc.se=vc.se[,ovc]
    vc.cum = vc.cum[,ovc]
    par(mar=c(12,4,4,2))
    bp=barplot(vc[1:3,], las=2, col=c('lightblue', 'lightgreen', 'pink' ), ylim=c(0,1), 
               ylab='fraction of phenotypic variance', legend=c('additive (A)', 'interaction (A x A)', 'residual strain repeatibility'), args.legend=list(x='topleft'))
    segments(bp-.225,  vc.cum[1,]-vc.se[1,], bp-.225, vc.cum[1,]+vc.se[1,], lwd=1.5, col='black')
    segments(bp,      vc.cum[2,]-vc.se[2,], bp,     vc.cum[2,]+vc.se[2,], lwd=1.5, col='black')
    segments(bp+.225,  vc.cum[3,]-vc.se[3,], bp+.225, vc.cum[3,]+vc.se[3,], lwd=1.5, col='black')
    return(ovc)
}
sortAVG = function(vc, vc.se, ovc=NULL) {
    colnames(vc)=gsub('1_', '', colnames(vc))
    colnames(vc)=gsub('_1', '', colnames(vc))
    # vc = bpAA
    # vc.se = bpAA.se
    vc.cum  =apply(vc, 2, cumsum)
    if(is.null(ovc)) {    ovc = order(vc.cum[2,]) }
    vc=vc[,ovc]
    vc.se=vc.se[,ovc]
    vc.cum = vc.cum[,ovc]
    par(mar=c(12,4,4,2))
    bp=barplot(vc[1:2,], las=2, col=c('lightblue', 'lightgreen' ), ylim=c(0,1), 
               ylab='fraction of phenotypic variance', legend=c('A', 'AA'), args.legend=list(x='topleft'))
    segments(bp-.05,  vc.cum[1,]-vc.se[1,], bp-.05, vc.cum[1,]+vc.se[1,], lwd=2, col='black')
    segments(bp,      vc.cum[2,]-vc.se[2,], bp,     vc.cum[2,]+vc.se[2,], lwd=2, col='black')
    return(ovc)
}

sortBP4 = function(vc, vc.se, ovc=NULL) {
   colnames(vc)=gsub('1_', '', colnames(vc))
   colnames(vc)=gsub('_1', '', colnames(vc))
   # vc = bpAA
   # vc.se = bpAA.se
    vc.cum  =apply(vc, 2, cumsum)
    if(is.null(ovc)) {    ovc = order(vc.cum[3,]) }
    vc=vc[,ovc]
    vc.se=vc.se[,ovc]
    vc.cum = vc.cum[,ovc]
    par(mar=c(12,4,4,2))
    bp=barplot(vc[1:4,], las=2, col=c('lightblue', 'lightgreen', 'orange','pink' ), ylim=c(0,1), 
               ylab='fraction of phenotypic variance', legend=c('A', 'AA' ,'AAA', 'R'), args.legend=list(x='topleft'))
    segments(bp-.05,  vc.cum[1,]-vc.se[1,], bp-.05, vc.cum[1,]+vc.se[1,], lwd=2, col='black')
    segments(bp,      vc.cum[2,]-vc.se[2,], bp,     vc.cum[2,]+vc.se[2,], lwd=2, col='black')
    segments(bp+.05,  vc.cum[3,]-vc.se[3,], bp+.05, vc.cum[3,]+vc.se[3,], lwd=2, col='black')
    segments(bp+.1,   vc.cum[4,]-vc.se[4,], bp+.1,  vc.cum[4,]+vc.se[4,], lwd=2, col='black')
    return(ovc)
}
bpA = (sapply(sigmas, function(x) x$mm.A))
bpA.se = (sapply(sigma.cov, function(x) x$mm.A))

bpAA=(sapply(sigmas, function(x) x$mm.AA))
bpAA.se=sapply(sigma.cov, function(x) x$mm.AA)

bpAAA=(sapply(sigmas, function(x) x$mm.AAA))
bpAAA.se=sapply(sigma.cov, function(x) x$mm.AAA)

bpAAp=(sapply(sigmas, function(x) x$mm.peaks.AA))
bpAAp.se=sapply(sigma.cov, function(x) x$mm.peaks.AA)

bpAApg=(sapply(sigmas, function(x) x$mm.peaks.AAg))
bpAApg.se=sapply(sigma.cov, function(x) x$mm.peaks.AAg)

pdf(file='~/Dropbox/4000BYxRM/Figures/VC_barplot_replicates.pdf', width=10, height=8)
ovc=sortBP(bpAA, bpAA.se)
dev.off()
x11()
sortBP4(bpAAA, bpAAA.se, ovc)
bpmean=t(t((apply(bpAA, 1, mean))))
bpp=barplot(t(t((apply(bpAA, 1, mean)))), las=2, col=c('lightblue', 'lightgreen', 'pink','white' ),
          ylab='fraction of phenotypic variance',
         legend=c('additive', 'two-way epistatic', 'Residual Strain Repeatibility'),args.legend=list(x='topleft')
         )
cvt=t(t(cumsum(apply(bpAA, 1, mean))))
bpp.se=t(t((apply(bpAA, 1, sd))))/sqrt(20)
segments(bpp,  cvt-bpp.se, bpp, cvt+bpp.se, lwd=2, col='black')

pdf(file='~/Dropbox/4000BYxRM/Figures/VC_pie_replicates.pdf', width=10, height=10)
pie(bpmean, labels=c('additive (A)', 'interaction (A x A)', 'residual strain repeatibility', ''),
    col=c('lightblue', 'lightgreen', 'pink','white' ))
dev.off()

pdf(file='~/Dropbox/4000BYxRM/VC_pie_replicates_scaled.pdf', width=10, height=10)
pie(bpmean[1:3,]/sum(bpmean[1:3,]), labels=c('additive (A)', 'two-way epistatic (A x A)', 'residual strain repeatibility'),
    col=c('lightblue', 'lightgreen', 'pink' ))
dev.off()

epistatic.table=rbind(bpAA[2,], bpAApg[2,],bpAAp[2,])
rownames(epistatic.table)=c('whole genome','additive QTL X genome', 'additive QTL')
epistatic.table.se=rbind(bpAA.se[2,], bpAApg.se[2,],bpAAp.se[2,])
colnames(epistatic.table.se)=colnames(epistatic.table)

et=epistatic.table
et.se=epistatic.table.se
colnames(et)=gsub('1_', '', colnames(et))
colnames(et)=gsub('_1', '', colnames(et))
et=et[,ovc]
et.se=et.se[,ovc]
par(mar=c(12,4,4,2))
bp=barplot(et, las=2, beside=T, legend=rownames(et), args.legend=list(x='topleft'),ylim=c(0, max(et)+.025),ylab='fraction of phenotypic variance')
segments(bp-.05,  et-et.se, bp-.05, et+et.se, lwd=2, col='black')
# Start 2 -locus section ##################################################################################################

peakfinder2D = function(i2d, threshold,pdist.m=500) {
    #   i2d   = iLODsh[trait,,]
        ipeaks  = which(i2d>threshold, arr.ind=T)
        ipeaks = cbind(ipeaks, i2d[ipeaks])
        ipeaks = cbind(ipeaks, rep(NA, nrow(ipeaks)))
        ipeaks = cbind(ipeaks, ipeaks[,3])
        colnames(ipeaks) = c('x', 'y', 'lod', 'group', 'glod')
        g=1
        while(sum(is.na(ipeaks[,'group']))>0) {
            peak = which.max(ipeaks[,'glod'])
            pdist = sqrt(abs(ipeaks[peak,'x']-ipeaks[,'x'])^2 + abs(ipeaks[peak,'y']-ipeaks[,'y'])^2)
            gpeaks = which(pdist< pdist.m & is.na(ipeaks[,'group']))
            if(length(gpeaks)>2) {
               ipeaks[gpeaks, 'group']=g
               ipeaks[gpeaks, 'glod']=0
               g=g+1
               print(g)
            } else{
               ipeaks[gpeaks, 'group']=0
               ipeaks[gpeaks, 'glod']=0
            }
        }
        ipeaks=ipeaks[ipeaks[,'group']>0,]
        ips = split(data.frame(ipeaks), ipeaks[,'group'])
        ipsp=data.frame(t(sapply(ips, function(x) {
                lmax= which.max(x$lod)
                x[lmax, c(1,2,3,4)]
                })))
        ipsp=        apply(ipsp, 2, as.numeric)
        return(ipsp)
}

nbr=matrix(NA,4390,20)
all.strain.names=names(pheno_raw[[1]])

# extract blup residuals and properly assign to strains##############################
nbr=matrix(NA,4390,20)
for(i in 1:length(pheno_raw) ) {
    ind.to.put=which(!is.na(cross$pheno[,i]))
    nbr[ind.to.put,i]=as.vector(newMM.blupResids[[i]])
}
colnames(nbr)=names(pheno_raw)

#testing with copper
copper.example=matrix(NA,4390,16)
for(i in 1:16 ) {
    ind.to.put=which(!is.na(cross$pheno[,2]))
    copper.example[ind.to.put,i]=as.vector(newMM.cblups[[2]][,i])
}
cop.lods.chr=get.LOD.by.COR(4276,copper.example, gdata)
#######################################################################################

# how scan was executed :::: 
#s2mat=matrix(0,28220, 28220)
#array(0, dim=c(10, 28220,28220))
for(i in 1:(max.marker-50)) {
    print(i)
    gint=gdata[,i]*gdata[,(i+50):max.marker]
    s2=get.LOD.by.COR(n.pheno, nbr, gint,doGPU=T)
   # save(s2, file=paste('/data/4000BYxRM/s2blup/', i, sep=''))
}

n.pheno   = countStrainsPerTrait(nbr)
#gdata.sub = gdata[,seq(1, ncol(gdata), 4)]
#rm(cross)
#rm(gdata)
gsub.name=do.call('rbind', strsplit(colnames(gdata) ,'_'))
gsub.pos=as.numeric(gsub.name[,1])
tick.spot=c(0, cumsum(rle(gsub.name[,2])$lengths)[-16])+1
tick.name=rle(gsub.name[,2])$values

gsub.ind.chr.split= split(1:nrow(gsub.name), match(gsub('chr', '', gsub.name[,2]), as.character(as.roman(c(1:16)))))


# Two locus scan with blup residuals ---------------------------------------------------------------
# calculate 2-locus model with residuals from additive polygenic model
ddata=list(gdata=gdata, nbr.plus.perms=nbr.plus.perms)
#save(ddata, file='/data/4000BYxRM/data.RData')
#load('/data/4000BYxRM/data.RData')
attach(ddata)
n.pheno.p=apply(nbr.plus.perms, 2, function(x) {sum(!is.na(x) )} )          
max.marker=ncol(gdata)
ddata=list(gdata=gdata, nbr=nbr)
#save(ddata, file='/data/4000BYxRM/dataBResid.RData')
#nbr.p=replicate(5, nbr[sample(1:4390),])
#nbr.plus.perms=cbind(nbr, nbr.p[,,1], nbr.p[,,2], nbr.p[,,3], nbr.p[,,4], nbr.p[,,5]) 
#save.image(file='/data/4000BYxRM/032314.RData')
#n.pheno.p = rep(n.pheno, 6)
#foreach(i=1:(max.marker-50)) %dopar% {
for(i in 1:(max.marker-50)) {
#for(i in 2883:6000) {
    #(max.marker-50)) {
    print(i)
    gint=gdata[,i]*gdata[,(i+50):max.marker]
    s2=get.LOD.by.COR(n.pheno.p, nbr.plus.perms, gint)
    #save(s2, file=paste('/data/4000BYxRM/s2blup/', i, sep=''))
}
##################33#############################################33################################


# 2D peak detection
blp.full.real.2D.peaks=list()
blp.marg.real.2D.peaks=list()
for(i in 1:20) {
    print(i)
    load(paste('/data/4000BYxRM/s2blupmat_real/', i, sep=''))
    blp.full.real.2D.peaks[[i]]=peakfinder2D(s2mat,3)
    s2.sp=Matrix(0, 28220, 28220, sparse=T)
    for( n in newMM.peaks[[i]] ) {
        print(n)
        s2.sp[n,] = s2mat[n,]
        s2.sp[,n] = s2mat[,n]
    }
    blp.marg.real.2D.peaks[[i]]=peakfinder2D(s2.sp,3)
    rm(s2mat)
}
#save(blp.full.real.2D.peaks, file = '/data/4000BYxRM/blp_full_real2Dpeaks.RData')
#save(blp.marg.real.2D.peaks, file = '/data/4000BYxRM/blp_marg_real2Dpeaks.RData')

blp.marg.null.2D.peaks=list()
blp.full.null.2D.peaks=list()
for(i in 21:120) {
    print(i)
    load(paste('/data/4000BYxRM/s2blupmat_perm/', i, sep=''))
    blp.full.null.2D.peaks[[i]]=peakfinder2D(s2mat,3)
    tval = i %%20
    if(tval==0) {tval=20 } 
    print(i)
    print(tval)
    s2.sp=Matrix(0, 28220, 28220, sparse=T)
    for( n in newMM.peaks[[tval]] ) {
        print(n)
        s2.sp[n,] = s2mat[n,]
        s2.sp[,n] = s2mat[,n]
    }
    blp.marg.null.2D.peaks[[i]]=peakfinder2D(s2.sp,3)
    rm(s2mat)
}
#save(blp.full.null.2D.peaks, file = '/data/4000BYxRM/blp_full_null2Dpeaks.RData')
#save(blp.marg.null.2D.peaks, file = '/data/4000BYxRM/blp_marg_null2Dpeaks.RData')


load('/data/4000BYxRM/blp_full_real2Dpeaks.RData') #blp.full.real.2D.peaks
load('/data/4000BYxRM/blp_full_null2Dpeaks.RData') #blp.full.null.2D.peaks
load('/data/4000BYxRM/blp_marg_real2Dpeaks.RData') #blp.marg.real.2D.peaks
load('/data/4000BYxRM/blp_marg_null2Dpeaks.RData') #blp.marg.null.2D.peaks

load('/data/4000BYxRM/raw_full_real2Dpeaks.RData') #blp.full.real.2D.peaks
load('/data/4000BYxRM/raw_full_null2Dpeaks.RData') #blp.full.null.2D.peaks
load('/data/4000BYxRM/raw_marg_real2Dpeaks.RData') #blp.marg.real.2D.peaks
load('/data/4000BYxRM/raw_marg_null2Dpeaks.RData') #blp.marg.null.2D.peaks


# FDR calculator ... assumes 20 traits and permutations labeled 21-120
doFDR.2locus = function(r, n) {
    perm.lods.s2=lapply(n,                      
                        function(x) { if(is.vector(x)){x[3]} else {x[,3]} } )
    iecnt=sapply(perm.lods.s2, function(x) { sapply(seq(2,7,.05), function(tt) {sum(x>tt)}) })
    rownames(iecnt)=seq(2,7,.05)

    p1=apply(iecnt[,21:40], 1, sum)
    p2=apply(iecnt[,41:60], 1, sum)
    p3=apply(iecnt[,61:80], 1, sum)
    p4=apply(iecnt[,81:100], 1, sum)
    p5=apply(iecnt[,101:120], 1, sum)
    iecnt=apply(cbind(p1,p2,p3,p4,p5), 1, mean)

    real.lods.s2=lapply(r, function(x) { x[,3] } )
    obscnt =sapply(real.lods.s2, function(x) { sapply(seq(2,7,.05), function(tt) {sum(x>tt)}) })
    rownames(obscnt)=seq(2,7,.05)
    obscnt=apply(obscnt, 1, sum)
    iecnt/obscnt
}

# 1 5 10 20 50
blp.full.fdr=doFDR.2locus(blp.full.real.2D.peaks,blp.full.null.2D.peaks)  #6.9 5.75 5.3 4.75 3.6
blp.marg.fdr=doFDR.2locus(blp.marg.real.2D.peaks,blp.marg.null.2D.peaks)   

raw.full.fdr=doFDR.2locus(raw.full.real.2D.peaks,raw.full.null.2D.peaks)
raw.marg.fdr=doFDR.2locus(raw.marg.real.2D.peaks,raw.marg.null.2D.peaks)

blp.full2D.fdr10=lapply(blp.full.real.2D.peaks, function(x) { x[x[,3]>5.3,] })
blp.marg2D.fdr10=lapply(blp.marg.real.2D.peaks, function(x) { x[x[,3]>3.75,] })

raw.full2D.fdr10=lapply(raw.full.real.2D.peaks, function(x) {(x[x[,3]>5.6,]) })
raw.marg2D.fdr10=lapply(raw.marg.real.2D.peaks, function(x) { x[x[,3]>4.2,] })


sum(sapply(blp.full2D.fdr10, nrow))
sum(sapply(blp.marg2D.fdr10, nrow))

sum(sapply(raw.full2D.fdr10, function(x) if(is.vector(x)){1} else { nrow(x) } ))
sum(sapply(raw.marg2D.fdr10, function(x) if(is.vector(x)){1} else { nrow(x) } ))


plot(blp.full2D.fdr10[[12]][,'x'], blp.full2D.fdr10[[12]][,'y'], cex=blp.full2D.fdr10[[12]][,'lod']/3, xaxs='i', yaxs='i',
     xlim=c(0,28220), ylim=c(0,28220) , main='Full 2D scan')
points(raw.full2D.fdr10[[12]][,'x'], raw.full2D.fdr10[[12]][,'y'], cex=raw.full2D.fdr10[[12]][,'lod']/3, xaxs='i', yaxs='i', col='red')
abline(0,1)

plot(blp.marg2D.fdr10[[12]][,'x'], blp.marg2D.fdr10[[12]][,'y'], cex=blp.marg2D.fdr10[[12]][,'lod']/3, xaxs='i', yaxs='i',
     xlim=c(0,28220), ylim=c(0,28220) , main='Marginal 2D scan')
points(raw.marg2D.fdr10[[12]][,'x'], raw.marg2D.fdr10[[12]][,'y'], cex=raw.marg2D.fdr10[[12]][,'lod']/3, xaxs='i', yaxs='i', col='red')

plot(blp.full2D.fdr10[[12]][,'x'], blp.full2D.fdr10[[12]][,'y'], cex=blp.full2D.fdr10[[12]][,'lod']/3, xaxs='i', yaxs='i',
     xlim=c(0,28220), ylim=c(0,28220) , main='Full vs Marginal 2D scan')
points(blp.marg2D.fdr10[[12]][,'x'], blp.marg2D.fdr10[[12]][,'y'], cex=blp.marg2D.fdr10[[12]][,'lod']/3, xaxs='i', yaxs='i', col='red')

plot(blp.full2D.fdr10[[1]][,'x'], blp.full2D.fdr10[[1]][,'y'], cex=blp.full2D.fdr10[[1]][,'lod']/3, xaxs='i', yaxs='i',
     xlim=c(0,28220), ylim=c(0,28220) , main='Full vs Marginal 2D scan')
points(blp.marg2D.fdr10[[1]][,'x'], blp.marg2D.fdr10[[1]][,'y'], cex=blp.marg2D.fdr10[[1]][,'lod']/3, xaxs='i', yaxs='i', col='red')

j=20
plot(blp.full2D.fdr10[[j]][,'x'], blp.full2D.fdr10[[j]][,'y'], cex=blp.full2D.fdr10[[j]][,'lod']/3, xaxs='i', yaxs='i',
     xlim=c(0,28220), ylim=c(0,28220) , main='Full vs Marginal 2D scan')
points(blp.marg2D.fdr10[[j]][,'x'], blp.marg2D.fdr10[[j]][,'y'], cex=blp.marg2D.fdr10[[j]][,'lod']/3, xaxs='i', yaxs='i', col='red')


c2= apply(cross$pheno, 2, scale)

#f2d=list()
f2dm=list()
for(i in 1:20) {
    # bit for additive QTL  --------------------------------------------
    apeaks = newMM.peaks[[i]]
    ax= paste('gdata[,', apeaks,']', sep='')
    aq=paste(ax, collapse= ' + ')
    am=lm(paste('c2[,' , i,']' , ' ~ ', (aq), '-1'))
    aov.a = anova(am)
    tssq=sum(aov.a[,2])
    a.effs=(aov.a[1:(nrow(aov.a)-1),2]/tssq)
    va.a = (tssq-aov.a[nrow(aov.a),2])/(tssq)
    ap.pos = do.call('rbind', strsplit(colnames(gdata)[apeaks], '_'))[,c(2,3)]
    colnames(ap.pos)=c('chr', 'pos')
    aq.table=data.frame(index=apeaks, ap.pos, lm_coefficient=coefficients(am), variance_explained=a.effs,row.names=NULL) 
  
    #f2d[[i]]=list(a=va.a, am=am, aq.table=aq.table, ai=0, i=0, n.int=0,iq.table=rep(NA,12) )
    f2dm[[i]]=list(a=va.a, am=am, aq.table=aq.table, ai=0, i=0, n.int=0,iq.table=rep(NA,12) )
    if(i==7) next ; 
    # -------------------------------------------------------------------
    #!!!!!!!@@@@@@@!! for full 2d (f2d)
    #ipeaks = blp.full2D.fdr10[[i]]
    #!!!!!!!!@@@@@@!! for marginal 2d (f2dm)
     ipeaks = blp.marg2D.fdr10[[i]]

    gint=gdata[,ipeaks[,'x']]*gdata[,ipeaks[,'y']]

    snapit.x = sapply(ipeaks[,'x'], function(x) { abs(x-apeaks) })
    snapit.y = sapply(ipeaks[,'y'], function(x) { abs(x-apeaks) })

    ca.x = apeaks[apply(snapit.x, 2, which.min)]
    ca.y = apeaks[apply(snapit.y, 2, which.min)]

    d.x=apply(snapit.x, 2, min)
    d.y= apply(snapit.y, 2, min)

    ix=paste('gdata[,', ipeaks[,'x'] ,']', sep='')
    iy=paste('gdata[,', ipeaks[,'y'] ,']', sep='')

    iq= paste(ix, '*', iy, sep='')
    iq= paste(iq, collapse = ' + ') 
    im=lm(paste('c2[,' , i,']' , ' ~ ', aq, '+', iq , '-1'))

    aov.i = anova(im)
    int.ind= grep(':', rownames(aov.i))
    i.effs=(aov.i[int.ind,2]/tssq)

    # calculate only variance from interactions
    i.ve=sum(aov.i[int.ind,2])/tssq
    ai.ve=sum(aov.i[1:(min(int.ind)-1),2])/tssq
    #va.i = (tssq-aov.i[nrow(aov.i),2])/(tssq)
    n.int =length(int.ind)

    ip.pos1 = do.call('rbind', strsplit(colnames(gdata)[ipeaks[,'x']], '_'))[,c(2,3)]
    colnames(ip.pos1)=c('Q1_chr', 'Q1_pos')
    ip.pos2 = do.call('rbind', strsplit(colnames(gdata)[ipeaks[,'y']], '_'))[,c(2,3)]
    colnames(ip.pos2)=c('Q2_chr', 'Q2_pos')
    iq.table=data.frame(Q1_index=ipeaks[,'x'], Q2_index=ipeaks[,'y'],  ip.pos1, ip.pos2,LOD=ipeaks[,'lod'],
                        lm_coefficient=coefficients(im)[int.ind], variance_explained=i.effs, closest.additive.Q1=ca.x, closest.additive.Q2=ca.y,
                        QQ.class=((d.x>200) + (d.y>200)),   row.names=NULL) 

    # gint is very important for the VC analysis
   # f2d[[i]]=list(a=va.a, am=am, aq.table=aq.table, ai=ai.ve, i = i.ve, n.int=n.int, gint=gint, im=im, iq.table=iq.table)
    f2dm[[i]]=list(a=va.a, am=am, aq.table=aq.table, ai=ai.ve, i = i.ve, n.int=n.int, gint=gint, im=im, iq.table=iq.table)
}

#names(f2d) =  names(pheno_raw)
names(f2dm) = names(pheno_raw)

save(f2d, file='~/Dropbox/4000BYxRM/f2d.RData')
save(f2dm, file='~/Dropbox/4000BYxRM/f2dm.RData')

# analysis of residuals #####################################################
residfm=matrix(NA, 4390, 19)
rownames(residfm)=seq(1,4390)
all.resids=lapply(f2dm[c(1:6, 8:20)], function(x) { residuals(x$im) })
for (i in 1:19) {
    residfm[as.numeric(names(all.resids[[i]])),i]=as.vector(all.resids[[i]])
}
#residfm[is.na(residfm)]=-10
bad.segs=unique(sort(which(is.na(residfm), arr.ind=T)[,1]))
residfm=residfm[-bad.segs, ]
dresid = dist(residfm)
hd= hclust(dresid)
plot(hd)
#############################################################################

# table of additive QTL
add.qtl.table=do.call('rbind', lapply(f2d, function(x) x$aq.table))
nn.add.qtl.table=gsub('_1.*', '', rownames(add.qtl.table))
nn.add.qtl.table=gsub('1_', '', (nn.add.qtl.table))
add.qtl.table=data.frame(trait=nn.add.qtl.table, add.qtl.table)
write.table(add.qtl.table, file='~/Dropbox/4000BYxRM/add_qtl_table.txt', sep='\t', quote=F,row.names=F)


# table of interacting QTL
miq.qtl.table=do.call('rbind', lapply(f2dm[c(1:6,8:20)], function(x) x$iq.table))
nn.miq.qtl.table=gsub('_1.*', '', rownames(miq.qtl.table))
nn.miq.qtl.table=gsub('1_', '', (nn.miq.qtl.table))
miq.qtl.table=data.frame(trait=nn.miq.qtl.table, miq.qtl.table)
write.table(miq.qtl.table, file='~/Dropbox/4000BYxRM/marginal_int_qtl_table.txt', sep='\t', quote=F,row.names=F)

iq.qtl.table=do.call('rbind', lapply(f2d[c(1:6,8:20)], function(x) x$iq.table))
nn.iq.qtl.table=gsub('_1.*', '', rownames(iq.qtl.table))
nn.iq.qtl.table=gsub('1_', '', (nn.iq.qtl.table))
iq.qtl.table=data.frame(trait=nn.iq.qtl.table, iq.qtl.table)
write.table(iq.qtl.table, file='~/Dropbox/4000BYxRM/full2d_int_qtl_table.txt', sep='\t', quote=F,row.names=F)

#calculate overlap
m2dtable.by.trait= split(miq.qtl.table, miq.qtl.table$trait)
f2dtable.by.trait= split(iq.qtl.table, iq.qtl.table$trait)
for(i in 1:19) {
    #m2dtable.by.trait[[i]]$f2d.Q1o=
    a = sapply(m2dtable.by.trait[[i]]$Q1_index, function(x){ (abs(x-f2dtable.by.trait[[i]]$Q1_index))  })
    b = sapply(m2dtable.by.trait[[i]]$Q2_index, function(x){ (abs(x-f2dtable.by.trait[[i]]$Q2_index))  })
    cc = apply(a+b,2,min)
    c.ind=apply(a+b,2,which.min)
    m2dtable.by.trait[[i]]$f2d.QQmin.ind=c.ind
    m2dtable.by.trait[[i]]$f2d.QQmin.distsum=cc
    m2dtable.by.trait[[i]]$f2d.o= m2dtable.by.trait[[i]]$f2d.QQmin.distsum<400
   #m2dtable.by.trait[[i]]$f2d.Q2o=sapply(m2dtable.by.trait[[i]]$Q2_index, function(x){ (abs(x-f2dtable.by.trait[[i]]$Q2_index))  })
   #m2dtable.by.trait[[i]]$f2d.QQdist=(m2dtable.by.trait[[i]]$f2d.Q1o<200 & m2dtable.by.trait[[i]]$f2d.Q2o<200)
}
for(i in 1:19) {
    a=sapply(f2dtable.by.trait[[i]]$Q1_index, function(x){(abs(x-m2dtable.by.trait[[i]]$Q1_index))})
    b=sapply(f2dtable.by.trait[[i]]$Q2_index, function(x){(abs(x-m2dtable.by.trait[[i]]$Q2_index))})
    cc = apply(a+b,2,min)
    c.ind=apply(a+b,2,which.min)
    f2dtable.by.trait[[i]]$m2d.QQmin.ind=c.ind
    f2dtable.by.trait[[i]]$m2d.QQmin.distsum=cc
    f2dtable.by.trait[[i]]$m2d.o=f2dtable.by.trait[[i]]$m2d.QQmin.distsum<400
}

sum(sapply(m2dtable.by.trait, function(x) sum(x$f2d.o))) #191 from marginal overlap the full scan    #out of 266 overlap from full 2D scan
sum(sapply(f2dtable.by.trait, function(x) sum(x$m2d.o))) 
write.table(do.call('rbind', m2dtable.by.trait), row.names=T, sep='\t', quote=F,  file='~/Desktop/m.txt')
write.table(do.call('rbind', f2dtable.by.trait), row.names=T, sep='\t', quote=F,  file='~/Desktop/f.txt')

sum(sapply(f2dtable.by.trait, function(x) sum(x$m2d.o)))
sum(sapply(m2dtable.by.trait, function(x) sum(x$f2d.o)))


hist(as.vector(unlist(sapply(m2dtable.by.trait, function(x) x$f2d.QQmin.distsum))) ,breaks=500, xlim=c(0,1000))

        
        #library(igraph)
#x=cbind(iq.qtl.table$closest.additive.Q1,iq.qtl.table$closest.additive.Q2)
#y=x[order(x[,1]),]
#plot(1:27819, seq(1,27819)/27819, type='n')
#segments(y[,1], 0, y[,2], 1 ) 
#segments(y[,2], 0, y[,1], 1 ) 


library(ggplot2)
pdf(file='~/Dropbox/4000BYxRM/Figures/QQvarexpbyclass.pdf', width=6, height=6)
# note use 'theme' ... need to figure this out
ggplot(iq.qtl.table, aes((factor(abs(QQ.class-2))), variance_explained))+geom_boxplot(notch=F,colour = I("#3366FF"), alpha=.5, outlier.shape = NA) + geom_point(position=position_jitter(width=0.1),shape=1, alpha=1)+
    ylab('fraction of phenotypic variance explained') + xlab('interactor(s) with a significant additive effect') +
    theme_bw() +     
    opts(axis.line = theme_segment(colour = "black"),
    panel.grid.major = theme_blank(),
    panel.grid.minor = theme_blank(),
    panel.border = theme_blank(),
    panel.background = theme_blank()) 
dev.off()

full.ve.int=sapply(f2d, function(x) x$i)
marg.ve.int=sapply(f2dm, function(x) x$i)

full.n.int=sapply(f2d, function(x) x$n.int)
marg.n.int=sapply(f2dm, function(x) x$n.int)


####VC analysis on trait averages 


avgMM.results=list()
### recalculate variance components on the scale of average trait values 
for(phenotype in names(pheno_raw) ) {
    print(phenotype)
    s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'

    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------

    y=ny[!is.na(ny)]

    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 

    strains.with.phenos=match(unique.sn, all.strain.names)

    n.strains=length(strains.with.phenos)
   
    # calculate relatedness matrix from genomewide marker data
    # also, in this case, this is equivalent to (gdata %*% t(gdata))
    A  = A.mat(gdata[strains.with.phenos,], shrink=FALSE)/2
    AA = A*A
    AAA = A*A*A
    y.avg=as.vector(by(y, names(y), mean, na.rm=T))
    names(y.avg)=strains.with.phenos
    peaks=newMM.peaks[[phenotype]]
    
    #peak.A = A.mat(gdata[strains.with.phenos, sort(sample(ncol(gdata), length(peaks)))], shrink=FALSE)/2
    peak.A = A.mat(gdata[strains.with.phenos,peaks], shrink=FALSE)/2
    pAA    = peak.A*peak.A
    pAG    = peak.A *A 

    mm.A =calcMM(y.avg, B=list(A=A), alg='fs',reps=FALSE)
    mm.A = extractVarCompResultsJB(mm.A)

    mm.peaks.A =calcMM(y.avg, B=list(A=peak.A), alg='fs',reps=FALSE)
    mm.peaks.A = extractVarCompResultsJB(mm.peaks.A)

    mm.AA =calcMM(y.avg, B=list(A=A,AA=AA), alg='fs',reps=FALSE)
    mm.AA = extractVarCompResultsJB(mm.AA)
    
   # mm.AAA=calcMM(y.avg, B=list(A=A, AA=AA, AAA=AAA), alg='fs', reps=FALSE)
   # mm.AAA = extractVarCompResultsJB(mm.AAA)

   mm.peaks.AA = calcMM(y.avg, B=list(A=peak.A,AA=pAA), alg='fs',reps=FALSE)
   mm.peaks.AA = extractVarCompResultsJB(mm.peaks.AA)
    
   mm.peaks.AAg = calcMM(y.avg, B=list(A=peak.A,AA=pAG), alg='fs',reps=FALSE)
   mm.peaks.AAg = extractVarCompResultsJB(mm.peaks.AAg)

    # new 3/10/15
    mm.G.QQ = calcMM(y.avg, B=list(A=A, AA=pAA), alg='ai',reps=FALSE)
    mm.G.QQ=extractVarCompResultsJB( mm.G.QQ )
    avgMM.results[[phenotype]]$mm.G.QQ=mm.G.QQ
   
    mm.G.GQ = calcMM(y.avg, B=list(A=A, AA=pAG), alg='ai',reps=FALSE)
    mm.G.GQ=extractVarCompResultsJB(  mm.G.GQ )
    avgMM.results[[phenotype]]$mm.G.GQ=mm.G.GQ
    #####------------------------------------


   avgMM.results[[phenotype]]=list(mm.AA = mm.AA, 
                                    mm.peaks.AA =  mm.peaks.AA, 
                                    mm.peaks.AAg = mm.peaks.AAg,
                                    mm.A=mm.A,
                                   mm.peaks.A=mm.peaks.A)
}
#save(file='~/Dropbox/new_mm/4000BYxRM/Results_avgVC.RData',avgMM.results )

# Use mixed model to calculate variance explained by significant interacting QTL 
for(phenotype in names(pheno_raw) ) {
    if(phenotype ==  "1_Hydroxyurea_1") next;
    print(phenotype)
    s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'
    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------
    y=ny[!is.na(ny)]

    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    strains.with.phenos=match(unique.sn, all.strain.names)
    n.strains=length(strains.with.phenos)
   
    # calculate relatedness matrix from genomewide marker data
    # also, in this case, this is equivalent to (gdata %*% t(gdata ))
    A  = A.mat(gdata[strains.with.phenos,], shrink=FALSE)/2
    IQm = A.mat(f2dm[[phenotype]]$gint[strains.with.phenos,], shrink=FALSE)/2
    IQf = A.mat(f2d[[phenotype]]$gint[strains.with.phenos,],  shrink=FALSE)/2
    #IQ2 = (f2dm[[phenotype]]$gint %*% t(f2dm[[phenotype]]$gint))/ncol(f2dm[[phenotype]]$gint) 
   
    y.avg=as.vector(by(y, names(y), mean, na.rm=T))
    names(y.avg)=strains.with.phenos
    
    peaks=newMM.peaks[[phenotype]]
    peak.A = A.mat(gdata[strains.with.phenos,peaks], shrink=FALSE)/2
    pAA    = peak.A*peak.A
    pAG    = peak.A *A 

    mm.sigQQm.AA = calcMM(y.avg, B=list(A=A, AA=IQm), alg='fs',reps=FALSE)
    mm.sigQQm.AA = extractVarCompResultsJB(mm.sigQQm.AA)

    mm.sigQQf.AA = calcMM(y.avg, B=list(A=A, AA=IQf), alg='fs',reps=FALSE)
    mm.sigQQf.AA = extractVarCompResultsJB(mm.sigQQf.AA)

    avgMM.results[[phenotype]]$mm.sigQQm.AA= mm.sigQQm.AA
    avgMM.results[[phenotype]]$mm.sigQQf.AA= mm.sigQQf.AA
}

avgMM.results[["1_Hydroxyurea_1"]]$mm.sigQQm.AA$sigma=matrix(c(0,0,0))
avgMM.results[["1_Hydroxyurea_1"]]$mm.sigQQm.AA$sigma.cov=rbind(c(0,0,0),c(0,0,0),c(0,0,0))
avgMM.results[["1_Hydroxyurea_1"]]$mm.sigQQf.AA$sigma=matrix(c(0,0,0))
avgMM.results[["1_Hydroxyurea_1"]]$mm.sigQQf.AA$sigma.cov=rbind(c(0,0,0),c(0,0,0),c(0,0,0))


#load('~/Dropbox/new_mm/4000BYxRM/Results_avgVC.RData')
#avgMM.results )
#Extract and rescale VCs 
avg.sigmas=lapply(avgMM.results, function(x) lapply(x, function(y) y$sigma ))
avg.sigma.cov=lapply(avgMM.results, function(x) lapply(x, function(y) sqrt(diag(y$sigma.cov)) ))
avg.nf=lapply(avg.sigmas, function(x) {
             nf= sapply(x, sum)
             nf[nf<.0001]=1
             return(nf)
                                   } )
for(i in 1:20){
    for(j in 1:9) {
    avg.sigmas[[i]][[j]]=avg.sigmas[[i]][[j]]/avg.nf[[i]][j]
    avg.sigma.cov[[i]][[j]]=avg.sigma.cov[[i]][[j]]/avg.nf[[i]][j]
    }
}

# Avg MM table
avgMM.table=do.call('cbind', lapply(avg.sigmas, function(x) {do.call('rbind', x) }))
colnames(avgMM.table) = names(avg.sigmas)
colnames(avgMM.table)=gsub('1_', '', colnames(avgMM.table))
colnames(avgMM.table)=gsub('_1', '', colnames(avgMM.table))
#"mm.broad"     "mm.A"         "mm.AA"        "mm.peaks.AA"  "mm.peaks.AAg" "mm.AAA"
write.table(cbind(rep(names(avg.sigmas[[1]]), as.vector(sapply(avg.sigmas[[1]], nrow))), avgMM.table), file='~/Dropbox/4000BYxRM/avgMM_table_v2.txt', sep='\t', quote=F,row.names=F)
avgMM.se.table=do.call('cbind', lapply(avg.sigma.cov, function(x) do.call('c', x)))
write.table(cbind(rep(names(avg.sigmas[[1]]), as.vector(sapply(avg.sigmas[[1]], nrow))), avgMM.se.table), file='~/Dropbox/4000BYxRM/avgMM_table_se_v2.txt', sep='\t', quote=F,row.names=T)


# VC on avg data
avgA=(sapply(avg.sigmas, function(x) x$mm.A))
avgA.se=sapply(avg.sigma.cov, function(x) x$mm.A)

avgAp=(sapply(avg.sigmas, function(x) x$mm.peaks.A))
avgAp.se=sapply(avg.sigma.cov, function(x) x$mm.peaks.A)


avgAA=(sapply(avg.sigmas, function(x) x$mm.AA))
avgAA.se=sapply(avg.sigma.cov, function(x) x$mm.AA)

#avgAAA=(sapply(avg.sigmas, function(x) x$mm.AAA))
#avgAAA.se=sapply(avg.sigma.cov, function(x) x$mm.AAA)

#these could be replaced
# V1 with Q as additive term
#avgAAp=(sapply(avg.sigmas, function(x) x$mm.peaks.AA))
#avgAAp.se=sapply(avg.sigma.cov, function(x) x$mm.peaks.AA)

#avgAApg=(sapply(avg.sigmas, function(x) x$mm.peaks.AAg))
#avgAApg.se=sapply(avg.sigma.cov, function(x) x$mm.peaks.AAg)

# V2 with G as additive term
avgAAp=(sapply(avg.sigmas, function(x) x$mm.G.QQ))
avgAAp.se=sapply(avg.sigma.cov, function(x) x$mm.G.QQ)

avgAApg=(sapply(avg.sigmas, function(x) x$mm.G.GQ))
avgAApg.se=sapply(avg.sigma.cov, function(x) x$mm.G.GQ)

#ovc.a=sortAVG(avgAA[1:2,], avgAA.se[1:2,],ovc)
#x11()
#sortBP4(avgAAA, avgAAA.se, ovc.a)

#avg.epistatic.table=rbind(avgAA[2,], avgAApg[2,],avgAAp[2,])
#rownames(avg.epistatic.table)=c('whole genome','additive QTL X genome', 'additive QTL')
#avg.epistatic.table.se=rbind(avgAA.se[2,], avgAApg.se[2,],avgAAp.se[2,])
#colnames(avg.epistatic.table.se)=colnames(avg.epistatic.table)

#avg.et=avg.epistatic.table
#avg.et.se=avg.epistatic.table.se
#colnames(avg.et)=gsub('1_', '', colnames(avg.et))
#colnames(avg.et)=gsub('_1', '', colnames(avg.et))
#avg.et=avg.et[,ovc]
#avg.et.se=avg.et.se[,ovc]
#par(mar=c(12,4,4,2))
#avg=barplot(avg.et, las=2, beside=T, legend=rownames(avg.et), args.legend=list(x='topleft'),ylim=c(0, max(avg.et)+.025),ylab='fraction of phenotypic variance', main='strain averages')
#segments(avg-.05,  avg.et-avg.et.se, avg-.05, avg.et+avg.et.se, lwd=2, col='black')

# apply(avg.epistatic.table, 1, min)                                                                                                                                                                                               
#                 Additive (G) x Additive (G)                  Additive (Q) x Additive (G)                  Additive (Q) x Additive (Q)       significant QTL-QTL from marginal scan significant QTL-QTL from full two-locus scan 
#                                    0.040347                                     0.020879                                     0.006473                                     0.000000                                     0.000000 
#R> apply(avg.epistatic.table, 1, max)                                                                                                                                                                                               
#                 Additive (G) x Additive (G)                  Additive (Q) x Additive (G)                  Additive (Q) x Additive (Q)       significant QTL-QTL from marginal scan significant QTL-QTL from full two-locus scan 
#                                      0.2268                                       0.2152                                       0.1614                                       0.1275                                       0.1313 
#R> apply(avg.epistatic.table, 1, median)                                                                                                                                                                                            
#                 Additive (G) x Additive (G)                  Additive (Q) x Additive (G)                  Additive (Q) x Additive (Q)       significant QTL-QTL from marginal scan significant QTL-QTL from full two-locus scan 
#                                     0.10247                                      0.10216                                      0.04476                                      0.03934                                      0.03875 
#apply(avg.epistatic.table, 1, mean)
#                 Additive (G) x Additive (G)                  Additive (Q) x Additive (G)                  Additive (Q) x Additive (Q)       significant QTL-QTL from marginal scan significant QTL-QTL from full two-locus scan 
#                                     0.10030                                      0.10266                                      0.06038                                      0.04535                                      0.04350 

# VC analysis for interactions
library(gplots)
pdf(file='~/Dropbox/4000BYxRM/Figures/add_var_exp.pdf', width=10, height=10)
plotCI( avgA[1,], avgAp[1,], 
        uiw=avgA.se[1,] ,err='x',
        xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', sfrac=0, gap=.85,
     xlab='additive variance (whole genome)', ylab='phenotypic variance captured by QTL', barcol='darkgrey')
plotCI( avgA[1,], avgAp[1,], 
        uiw=avgAp.se[1,] ,err='y',
        xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i', sfrac=0, gap=.85,
       barcol='grey', add=T)
abline(0,1)
dev.off()

load('~/Dropbox/4000BYxRM/f2d.RData')
load('~/Dropbox/4000BYxRM/f2dm.RData')
full.n.int=sapply(f2d, function(x) x$n.int)
marg.n.int=sapply(f2dm, function(x) x$n.int)

avgAAmq=(sapply(avg.sigmas, function(x) x$mm.sigQQm.AA))
avgAAmq.se=sapply(avg.sigma.cov, function(x) x$mm.sigQQm.AA)

avgAAfq=(sapply(avg.sigmas, function(x) x$mm.sigQQf.AA))
avgAAfq.se=sapply(avg.sigma.cov, function(x) x$mm.sigQQf.AA)


avg.epistatic.table=rbind(avgAA[2,], avgAApg[2,],avgAAp[2,], 
                          #marg.ve.int, full.ve.int, 
                          avgAAmq[2,],avgAAfq[2,] )
rownames(avg.epistatic.table)=c('Genome x Genome',
                                'QTL x Genome',
                                'QTL x QTL', 
                                                                'significant QTL-QTL from marginal scan',
                                'significant QTL-QTL from full two-locus scan')
##'least squares sig QTL-QTL marginal',
                                ##'least squares sig QTL-QTL full',

avg.epistatic.table.se=rbind(avgAA.se[2,], avgAApg.se[2,],avgAAp.se[2,],
                             #rep(0,20), rep(0,20),
                             avgAAmq.se[2,],avgAAfq.se[2,] )
colnames(avg.epistatic.table)=paste(colnames(avg.epistatic.table), ' (', marg.n.int, ',', full.n.int, ')')
colnames(avg.epistatic.table.se)=colnames(avg.epistatic.table)

avg.et=avg.epistatic.table
oe = order(avg.et[1,])
avg.et.se=avg.epistatic.table.se
colnames(avg.et)=gsub('1_', '', colnames(avg.et))
colnames(avg.et)=gsub('_1', '', colnames(avg.et))
#universal ordering
#avg.et=avg.et[,ovc]
#avg.et.se=avg.et.se[,ovc]
# relative ordering
avg.et=avg.et[,oe]
avg.et.se=avg.et.se[,oe]

pdf(file='~/Dropbox/4000BYxRM/Figures/epistatic_var_full_breakdown_barplot.pdf', width=15, height=10)
par(mar=c(13,4,4,2))

                                           #'pink', 'purple',
avg=barplot(avg.et, las=2, beside=T, col=c('grey', 'lightblue', 'lightgreen', 
                                           'orange','violet'),
          legend=rownames(avg.et), args.legend=list(x='topleft'),ylim=c(0, max(avg.et)+.025),ylab='fraction of phenotypic variance', main='')
segments(avg,  avg.et-avg.et.se, avg, avg.et+avg.et.se, lwd=1, col='black')
dev.off()

aem=apply(avg.epistatic.table[,-7], 1, mean)
ae.se=apply(avg.epistatic.table[,-7], 1, sd)/sqrt(20)

pdf(file='~/Dropbox/4000BYxRM/Figures/epistatic_var_avg_barplot.pdf', width=10/1.5, height=16/1.5)
par(mar=c(20,8,8,4))
ae.bp=barplot(aem, beside=T, col=c('grey', 'lightblue', 'lightgreen', 'orange','violet'), ylim=c(0, .13), las=2, 
              ylab='fraction of phenotypic variance', main='')
segments(ae.bp,  aem-ae.se, ae.bp, aem+ae.se, lwd=1, col='black')
dev.off()


ape=c(mean(avgA[1,]), mean(avgAp[1,]), aem)
ape.se=c(sd(avgA[,1])/20, sd(avgAp[,1])/20, ae.se)
names(ape)=c('whole genome A', 'QTL A', paste(names(aem), 'AA'))
par(mar=c(20,8,8,4))
ape.bp=barplot(ape, beside=T, #col=c('black', 'blue', 'grey', 'lightblue', 'lightgreen', 'pink', 'purple','orange','violet')
                col=c('black', 'black', 'grey','grey','grey','grey','grey','grey','grey') 
               , ylim=c(0, .5), las=2, 
              ylab='fraction of phenotypic variance', main='variance summary')
segments(ape.bp,  ape-ape.se, ape.bp, ape+ape.se, lwd=1, col='black')

hist(avg.epistatic.table[6,]/avg.epistatic.table[1,], xlab='fraction of epistatic variance explained by sig QTL-QTL', main='')

# histograms of effect size
#a.effs=do.call('c', sapply(f2dm, function(x) x$a.effs))
#i.effs=do.call('c', sapply(f2dm, function(x) x$i.effs))
a.effs=do.call('c', sapply(f2dm, function(x) x$aq.table$variance_explained))
i.effs=do.call('c', sapply(f2dm[c(1:6,8:20)], function(x) x$iq.table[,'variance_explained']))
nnn.additive.qtl=sapply(newMM.peaks, length)
hist(nnn.additive.qtl, xlab='number of additive QTL', ylab='traits', main='')
# marginal nint
nnn.int.qtl=sapply(f2dm, function(x) {tryCatch(x$n.int, error=function(e){0})} )
nnn.int.qtl[[7]]=0
nnn.int.qtl=unlist(nnn.int.qtl)
hist(nnn.int.qtl, xlab='number of QTL-QTL', ylab='traits', main='')

# fulll nint
fnn.int.qtl=sapply(f2d, function(x) {tryCatch(x$n.int, error=function(e){0})} )
fnn.int.qtl[[7]]=0
fnn.int.qtl=unlist(fnn.int.qtl)
hist(fnn.int.qtl, xlab='number of QTL-QTL', ylab='traits', main='')

library(qtlDesign)
siglev = pchisq(2.5*4.60517,1,lower.tail=FALSE)
siglev.int = pchisq(3.85*4.60517,1,lower.tail=FALSE)
dd     = seq(0.01,10,.01)
n16000 = power.t.test(n=8000,  delta=seq(0.01,10,.01),sig.level=siglev)$power
n16000.int = power.t.test(n=8000,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
n8000 = power.t.test(n=4000,  delta=seq(0.01,10,.01),sig.level=siglev)$power
n8000.int = power.t.test(n=4000,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
n4000 = power.t.test(n=2000,  delta=seq(0.01,10,.01),sig.level=siglev)$power
n4000.int = power.t.test(n=2000,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
n1000 = power.t.test(n=500,  delta=seq(0.01,10,.01),sig.level=siglev)$power
n1000.int = power.t.test(n=500,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
n500 = power.t.test(n=250,  delta=seq(0.01,10,.01),sig.level=siglev)$power
n500.int = power.t.test(n=250,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power
n100 = power.t.test(n=50,  delta=seq(0.01,10,.01),sig.level=siglev)$power
n100.int = power.t.test(n=50,  delta=seq(0.01,10,.01),sig.level=siglev.int)$power

pvar= prop.var('ri', dd/2,1)
pdf(file='~/Dropbox/4000BYxRM/Figures/power_n_all.pdf', width=10, height=5)
par(xaxs='i', yaxs='i')
plot(pvar,n4000, type='l',col='red',lwd=2,xlim=c(0,0.05), ylab='power',
     xlab='phenotypic variance explained', xaxt='n')
axis(1, at=c(0,0.01, 0.02,0.03,0.04,0.05))#, 0.1, 0.15))
points(pvar,n4000.int, type='l',col='red',lwd=2, lty=2)

points(pvar,n1000, type='l',col='blue',lwd=2)
points(pvar,n1000.int, type='l',col='blue',lwd=2,lty=2)

#points(pvar,n500, type='l',col='purple',lwd=2)
#points(pvar,n500.int, type='l',col='purple',lwd=2, lty=2)

points(pvar,n100, type='l',col='green',lwd=2)
points(pvar,n100.int, type='l',col='green',lwd=2, lty=2)
#'n=500','n=100'),
legend('right', legend=c('n=4000', 'n=1000', 'n=100'),        fill=c('red', 'blue', 'green')) 
#  'purple', 'green'))
dev.off()

#pink', 'blue', 'cyan' 
#       ,'purple', 'violet', 'green', 'lightgreen'
#       ))

power4000 = data.frame(pvar, n4000)
p4000fx   = approxfun(n4000, pvar, rule=2)

power4000.int = data.frame(pvar, n4000.int)
p4000.intfx   = approxfun(n4000.int, pvar, rule=2)


dat = data.frame(variance_explained=c(sqrt(a.effs), sqrt(i.effs)), lines=              
                 c(rep('QTL', length(a.effs)), rep('QTL-QTL', length(i.effs))))
ggplot(dat, aes(x = variance_explained, fill = lines)) + geom_density(alpha = 0.2)

plot.multi.dens(list(Q=sqrt(a.effs), QQ=sqrt(i.effs)))

#ha=hist(a.effs,freq=F)
#hi=hist(i.effs, freq=F)

# histogram phenotypic var
x11()
par(mar = c(5,5,2,5))
par(xaxs='i', yaxs='i')
ha=hist((a.effs),col=rgb(0, 0, 0,0.3), main='effect size', xlab='fraction of phenotypic variance', breaks=100, xlim=c(0,0.03))
hi=hist((i.effs), add=T,col= rgb(1, 0, 0,0.8), breaks=20)
sf= max(ha$counts)
points((pvar), n4000*sf, type='l',col='blue',lwd=2, lty=1, xlab=NA, ylab=NA)
points((pvar), n4000.int*sf, type='l',col='red',lwd=2, lty=2,  xlab=NA, ylab=NA)
axis(side=4, at=seq(0,1, .2)*sf, labels= seq(0,1, .2))
mtext(side = 4, line = 3, "power")

# histogram abs allele diff
x11()
par(mar = c(5,5,2,5))
par(xaxs='i', yaxs='i')
ha=hist(sqrt(a.effs),col=rgb(0, 0, 0,0.3), main='effect size', xlab='fraction of phenotypic variance', breaks=100, xaxt='n',xlim=c(0,0.3))
hi=hist(sqrt(i.effs), add=T,col= rgb(1, 0, 0,0.8), breaks=50)
stp=.025
axis(side = 1,  at=seq(0, .5,stp), labels=round(seq(0, .5,stp)^2,3))
sf= max(ha$counts)
points(sqrt(pvar), n4000*sf, type='l',col='blue',lwd=2, lty=1, xlab=NA, ylab=NA)
points(sqrt(pvar), n4000.int*sf, type='l',col='red',lwd=2, lty=2,  xlab=NA, ylab=NA)
axis(side=4, at=seq(0,1, .2)*sf, labels= seq(0,1, .2))
mtext(side = 4, line = 3, "power")

# grey white 

# density phenotypic var (not square 
pdf(file='~/Dropbox/4000BYxRM/Figures/eff_dens.pdf', width=12, height=6)
#x11()
da=density((a.effs))
di=density((i.effs))
par(mar = c(4,4,2,4))
par(xaxs='i', yaxs='i')
plot(di, xlab='fraction of phenotypic variance', ylab='density', xlim=c(0, 0.025), main='', type='n')
polygon(da, col=rgb(0, 0, 1,0.3), border=NA)
polygon(di, col=rgb(1, 0, 0,0.3), border=NA)
sf= max(c(da$y, di$y))
points((pvar), n4000*sf, type='l',col='blue',lwd=3, lty=1, xlab=NA, ylab=NA)
points((pvar), n4000.int*sf, type='l',col='red',lwd=3, lty=1,  xlab=NA, ylab=NA)
axis(side=4, at=seq(0,1, .2)*sf, labels= seq(0,1, .2))
mtext(side = 4, line = 3, "power")
dev.off()

library(gplots)
#pdf(file='~/Dropbox/4000BYxRMvarcomp.pdf', width=10, height=10)
#plotCI( avg.epistatic.table[1,], avg.epistatic.table[6,], 
#        uiw=avg.epistatic.table.se[1,] ,err='x',
#        xlim=c(0,.27), ylim=c(0,.27), xaxs='i', yaxs='i', sfrac=0, gap=.85,
#     xlab='AA variance (whole genome)', ylab='AA variance (significant QTL-QTL)', barcol='darkgrey')
#plotCI( avg.epistatic.table[1,], avg.epistatic.table[6,],
#        uiw=avg.epistatic.table.se[6,] ,err='y',
#        xaxs='i', yaxs='i', sfrac=0, gap=.85,
#     xlab='A genome', ylab='A QTL', barcol='grey', add=T)
#abline(0,1)


pdf(file='/home/jbloom/Dropbox/4000BYxRM/Figures/2B_int_var_exp.pdf', width=10, height=10)
plotCI( avg.epistatic.table[1,], avg.epistatic.table[4,], 
        uiw=avg.epistatic.table.se[1,] ,err='x',
        xlim=c(0,.27), ylim=c(0,.27), xaxs='i', yaxs='i', sfrac=0, gap=.85,
     xlab='Genome x Genome variance', ylab='phenotypic variance captured by QTL-QTL', barcol='darkgrey')
plotCI( avg.epistatic.table[1,], avg.epistatic.table[4,],
        uiw=avg.epistatic.table.se[4,] ,err='y',
        xaxs='i', yaxs='i', sfrac=0, gap=.85,
     xlab='A genome', ylab='A QTL', barcol='grey', add=T)
abline(0,1)
dev.off()





load('~/Dropbox/new_mm/4000BYxRM/Results_avgVC.RData')
simVC.A=sapply(sim.avgMM.results, function(x) {x$mm.peaks.A[1]})
simVC.A=lapply(sim.avgMM.results, function(x) { do.call('cbind', x$mm.peaks.A)[1,] } )

bp.names=gsub('1_', '', names(sim.avgMM.results))
bp.names=gsub('_1', '', bp.names)
nnn.additive.qtl=sapply(newMM.peaks, length)
hist(nnn.additive.qtl, xlab='number of additive QTL', ylab='traits', main='')
# marginal nint
nnn.int.qtl=sapply(f2dm, function(x) {tryCatch(x$n.int, error=function(e){0})} )
nnn.int.qtl[[7]]=0
nnn.int.qtl=unlist(nnn.int.qtl)
hist(nnn.int.qtl, xlab='number of QTL-QTL', ylab='traits', main='')

# fulll nint
fnn.int.qtl=sapply(f2d, function(x) {tryCatch(x$n.int, error=function(e){0})} )
fnn.int.qtl[[7]]=0
fnn.int.qtl=unlist(fnn.int.qtl)
hist(fnn.int.qtl, xlab='number of QTL-QTL', ylab='traits', main='')



# some simulations

load('~/Dropbox/new_mm/4000BYxRM/Results.RData')
newMM.peaks=results[[1]]
sim.avgMM.results=list()
### recalculate variance components on the scale of average trait values 
for(phenotype in names(pheno_raw)[c(1:6,7:20)] ) {
    print(phenotype)
    s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'

    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------

    y=ny[!is.na(ny)]

    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 
    strains.with.phenos=match(unique.sn, all.strain.names)
    n.strains=length(strains.with.phenos)
       A  = A.mat(gdata[strains.with.phenos,], shrink=FALSE)/2
       y.avg=as.vector(by(y, names(y), mean, na.rm=T))
    names(y.avg)=strains.with.phenos
    peaks=newMM.peaks[[phenotype]]

    # uncomment for all sims  (Commented to run QQ slice)
    lmm.peaks.A=list()
    lmm.peaks.AA=list()
    lmm.peaks.AAg=list()
    lmm.peaks.QQ =list()
    for(i in 1:50) { 
        print(i)
        peak.A = A.mat(gdata[strains.with.phenos, sort(sample(ncol(gdata), length(peaks)))], shrink=FALSE)/2
        pAA    = peak.A*peak.A
        pAG    = peak.A *A 

        mm.peaks.A = regress(y.avg~1, ~peak.A, pos=c(T,T), verbose=F,tol=1e-3)
        lmm.peaks.A[[i]] = mm.peaks.A$sigma/sum(mm.peaks.A$sigma)
        print(phenotype)
        print(lmm.peaks.A[[i]])
        mm.peaks.AA = regress(y.avg~1, ~peak.A+pAA, pos=c(T,T,T), verbose=F,tol=1e-3)
        lmm.peaks.AA[[i]] = mm.peaks.AA$sigma/sum(mm.peaks.AA$sigma)
        print(phenotype)
        print(lmm.peaks.AA[[i]])

        mm.peaks.AAg =regress(y.avg~1, ~peak.A+pAG, pos=c(T,T,T), verbose=F,tol=1e-3)
        lmm.peaks.AAg[[i]] =  mm.peaks.AAg$sigma/sum(mm.peaks.AAg$sigma)
        print(phenotype)
        print(lmm.peaks.AAg[[i]])
        
        nint =(f2dm[[phenotype]]$n.int)

        pQQ.1 = sample(1:ncol(gdata), ifelse(nint==0, 1, nint ))
        pQQ.2 = sample(1:ncol(gdata), ifelse(nint==0, 1, nint ))
        gint.pQQ = gdata[,pQQ.1]* gdata[,pQQ.2]
        if(nint==0) {  QQ = A.mat(gint.pQQ[strains.with.phenos], shrink=FALSE)/2 }
        else {     QQ = A.mat(gint.pQQ[strains.with.phenos,], shrink=FALSE)/2 }

        mm.peaks.QQ =regress(y.avg~1, ~A+QQ, pos=c(T,T,T), verbose=T,tol=1e-3)
        lmm.peaks.QQ[[i]] =  mm.peaks.QQ$sigma/sum(mm.peaks.QQ$sigma)
        print(phenotype)
        print(lmm.peaks.QQ[[i]])

   }
   sim.avgMM.results[[phenotype]]=list(mm.peaks.A   = lmm.peaks.A,
                                        mm.peaks.AA  = lmm.peaks.AA, 
                                        mm.peaks.AAg = lmm.peaks.AAg,
                                        mm.peaks.QQ = lmm.peaks.QQ)
}
#save(sim.avgMM.results, file = '~/Dropbox/4000BYxRM/simVC_avgMM_Results.RData')
#save(sim.avgMM.results, file = '~/Dropbox/4000BYxRM/simVC_avgMM_Results1.RData')
load('~/Dropbox/4000BYxRM/simVC_avgMM_Results.RData')
load('~/Dropbox/4000BYxRM/simVC_avgMM_Results1.RData')
load('~/Dropbox/new_mm/4000BYxRM/Results_avgVC.RData')

nnn.additive.qtl=sapply(newMM.peaks, length)

simVC.A=sapply(sim.avgMM.results, function(x) {x$mm.peaks.A[1]})
simVC.A=lapply(sim.avgMM.results, function(x) { do.call('cbind', x$mm.peaks.A)[1,] } )

bp.names=gsub('1_', '', names(sim.avgMM.results))
bp.names=gsub('_1', '', bp.names)
bp.names=paste(bp.names, nnn.additive.qtl, sep=' : ')
#par(oma=c(5,1,1,1))

pdf(file='~/Dropbox/4000BYxRM/Figures/add_rand_VC.pdf', width=15, height=8)
par(oma=c(5,1,1,1))
x= boxplot(simVC.A, las=3, ylim=c(0,1), ylab='fraction of phenotypic variance', names=bp.names)
x=t(sapply(avgMM.results, function(x) x$mm.A$sigma))
mm.A.genome=x[,1]/rowSums(x)
x=t(sapply(avgMM.results, function(x) x$mm.peaks.A$sigma))
mm.peaks.A=x[,1]/rowSums(x)
points(mm.A.genome, col='blue',pch=20)
points(mm.peaks.A, col='red',pch=20)
legend('topright', legend=c('additive variance (whole genome)','additive variance captured by QTL','additive variance captured by randomly selelcted markers'), col=c('blue', 'red', 'black'), pch=20)
dev.off()

pdf(file='~/Dropbox/4000BYxRM/Figures/add_rand_VC_scatter.pdf', width=8, height=8)
par(oma=c(1,1,1,1))
plot( sapply(simVC.A, median)/mm.A.genome, mm.peaks.A/mm.A.genome,  xlim=c(0,1.1), ylim=c(0,1.1), xlab='median additive variance captured by randomly selected markers',
     ylab='additive variance captured by QTL')
abline(0,1)
dev.off()





simVC.A=lapply(sim.avgMM.results, function(x) { do.call('cbind', x$mm.peaks.AA)[2,] } )
pdf(file='~/Dropbox/4000BYxRM/Figures/int_rand_QTL_VC.pdf', width=15, height=8)
par(oma=c(5,1,1,1))
x= boxplot(simVC.A, las=3, ylim=c(0,.25), ylab='fraction of phenotypic variance', names=bp.names)
x=t(sapply(avgMM.results, function(x) x$mm.AA$sigma))
mm.A.genome=x[,2]/rowSums(x)
x=t(sapply(avgMM.results, function(x) x$mm.peaks.AA$sigma))
mm.peaks.A=x[,2]/rowSums(x)
points(mm.A.genome, col='blue',pch=20)
points(mm.peaks.A, col='red',pch=20)
legend('topright', legend=c('interaction variance (whole genome)','interaction variance (additive QTL)','interaction variance captured by randomly selelcted markers'), col=c('blue', 'red', 'black'), pch=20)
dev.off()

#plot( sapply(simVC.A, median)/mm.A.genome, mm.peaks.A/mm.A.genome,  xlim=c(0,1.1), ylim=c(0,1.1), xlab='median interaction variance captured by randomly selected markers',
#     ylab='interaction variance captured by additive QTL')
#abline(0,1)

simVC.A=lapply(sim.avgMM.results, function(x) { do.call('cbind', x$mm.peaks.AAg)[2,] } )
pdf(file='~/Dropbox/4000BYxRM/Figures/int_rand_QTLbyGenome_VC.pdf', width=15, height=8)
par(oma=c(5,1,1,1))
x= boxplot(simVC.A, las=3, ylim=c(0,.25), ylab='fraction of phenotypic variance', names=bp.names)
x=t(sapply(avgMM.results, function(x) x$mm.AA$sigma))
mm.A.genome=x[,2]/rowSums(x)
x=t(sapply(avgMM.results, function(x) x$mm.peaks.AAg$sigma))
mm.peaks.A=x[,2]/rowSums(x)
points(mm.A.genome, col='blue',pch=20)
points(mm.peaks.A, col='red',pch=20)
legend('topright', legend=c('interaction variance (whole genome)','interaction variance (additive QTL * genome)','interaction variance captured by randomly selelcted markers * genome'), col=c('blue', 'red', 'black'), pch=20)
dev.off()


simVC.A=lapply(sim.avgMM.results, function(x) { do.call('cbind', x$mm.peaks.QQ)[2,] } )
pdf(file='~/Dropbox/4000BYxRM/Figures/int_rand_QTLbyQTL_VC.pdf', width=15, height=8)
par(oma=c(5,1,1,1))
x= boxplot(simVC.A, las=3, ylim=c(0,.25), ylab='fraction of phenotypic variance', names=bp.names)
x=t(sapply(avgMM.results, function(x) x$mm.AA$sigma))
mm.A.genome=x[,2]/rowSums(x)
x=t(sapply(avgMM.results, function(x) x$mm.peaks.AAg$sigma))
mm.peaks.A=x[,2]/rowSums(x)
x=t(sapply(avgMM.results, function(x) x$mm.sigQQm.AA$sigma))
mm.peaks.AA=x[,2]/rowSums(x)
points(mm.A.genome, col='blue',pch=20)
points(mm.peaks.AA, col='red',pch=20)
#points(mm.peaks.A, col='orange',pch=20)

legend('topright', legend=c('interaction variance (whole genome)','interaction variance captured by significant QTL-QTL',
                            'interaction variance captured by randomly selected QTL-QTL'),
       col=c('blue', 'red', 'black'), pch=20)
dev.off()

simVC.A=lapply(sim.avgMM.results, function(x) { do.call('cbind', x$mm.peaks.QQ)[2,] } )
x=t(sapply(avgMM.results, function(x) x$mm.AA$sigma))
mm.AA.genome=x[,2]/rowSums(x)
#x=t(sapply(avgMM.results, function(x) x$mm.peaks.AAg$sigma))
#mm.peaks.A=x[,2]/rowSums(x)
x=t(sapply(avgMM.results, function(x) x$mm.sigQQm.AA$sigma))
mm.peaks.AA=x[,2]/rowSums(x)
#points(mm.A.genome, col='blue',pch=20)
#points(mm.peaks.AA, col='red',pch=20)

pdf(file='~/Dropbox/4000BYxRM/Figures/int_rand_QTLbyQTL_VC_scatter.pdf', width=8, height=8)
par(oma=c(1,1,1,1))
plot( sapply(simVC.A, median)/mm.AA.genome, mm.peaks.AA/mm.AA.genome,  xlim=c(0,1.1), ylim=c(0,1.1), xlab='median interaction variance captured by randomly selected markers',
     ylab='interaction variance captured by QTL-QTL')
abline(0,1)
dev.off()


load('~/Dropbox/new_mm/4000BYxRM/Results.RData')
newMM.peaks=results[[1]]
sim.avgMM.results=list()
### recalculate variance components on the scale of average trait values 
for(phenotype in names(pheno_raw) ) {
    print(phenotype)
    s=(pheno_raw[[phenotype]])
    # fix this -----------------------------------
    srle=sapply(s, length)
    sr=list()
    sr$lengths=as.vector(srle)
    sr$values=names(srle)
    attr(sr, 'class')='rle'

    ny=as.vector(unlist(s))
    names(ny)=inverse.rle(sr)
    #-----------------------------------------------

    y=ny[!is.na(ny)]

    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)
    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
   
    #for constructing Strain Variance component
    Strain      = diag(strain.cnt)
    Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 
    strains.with.phenos=match(unique.sn, all.strain.names)
    n.strains=length(strains.with.phenos)
       A  = A.mat(gdata[strains.with.phenos,], shrink=FALSE)/2
       y.avg=as.vector(by(y, names(y), mean, na.rm=T))
    names(y.avg)=strains.with.phenos
    peaks=newMM.peaks[[phenotype]]
    
    lmm.peaks.A=list()
    lmm.peaks.AA=list()
    lmm.peaks.AAg=list()
    for(i in 1:50) { 
        print(i)
        peak.A = A.mat(gdata[strains.with.phenos, sort(sample(ncol(gdata), length(peaks)))], shrink=FALSE)/2
        pAA    = peak.A*peak.A
        pAG    = peak.A *A 

        mm.peaks.A = regress(y.avg~1, ~peak.A, pos=c(T,T), verbose=F,tol=1e-3)
        lmm.peaks.A[[i]] = mm.peaks.A$sigma/sum(mm.peaks.A$sigma)
        print(phenotype)
        print(lmm.peaks.A[[i]])
        mm.peaks.AA = regress(y.avg~1, ~peak.A+pAA, pos=c(T,T,T), verbose=F,tol=1e-3)
        lmm.peaks.AA[[i]] = mm.peaks.AA$sigma/sum(mm.peaks.AA$sigma)
        print(phenotype)
        print(lmm.peaks.AA[[i]])

        mm.peaks.AAg =regress(y.avg~1, ~peak.A+pAG, pos=c(T,T,T), verbose=F,tol=1e-3)
        lmm.peaks.AAg[[i]] =  mm.peaks.AAg$sigma/sum(mm.peaks.AAg$sigma)
        print(phenotype)
        print(lmm.peaks.AAg[[i]])

   }
    sim.avgMM.results[[phenotype]]=list(mm.peaks.A   = lmm.peaks.A,
                                        mm.peaks.AA  = lmm.peaks.AA, 
                                        mm.peaks.AAg = lmm.peaks.AAg)
}
save(sim.avgMM.results, file = '~/Dropbox/4000BYxRM/simVC_avgMM_Results.RData')

#pdf(file='~/Dropbox/4000BYxRMvarcomp.pdf', width=10, height=10)
#plotCI( avg.epistatic.table[1,], avg.epistatic.table[6,], 
#        uiw=avg.epistatic.table.se[1,] ,err='x',
#        xlim=c(0,.27), ylim=c(0,.27), xaxs='i', yaxs='i', sfrac=0, gap=.85,
#     xlab='AA variance (whole genome)', ylab='AA variance (significant QTL-QTL)', barcol='darkgrey')
#plotCI( avg.epistatic.table[1,], avg.epistatic.table[6,],
#        uiw=avg.epistatic.table.se[6,] ,err='y',
#        xaxs='i', yaxs='i', sfrac=0, gap=.85,
#     xlab='A genome', ylab='A QTL', barcol='grey', add=T)
#abline(0,1)

