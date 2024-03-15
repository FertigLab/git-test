library(GEOquery)
sticMethDat=getGEO("GSE155311")
sticMethDat=sticMethDat[[1]]
hgscMethDat=getGEO("GSE155760")

sticBeta=exprs(sticMethDat)
sticClin=pData(phenoData(sticMethDat))

save.image()
save("sticBeta","sticClin",file="Data/sticDat.Rda")


#TCGA methylation data
library(data.table)
methFile="~/PROJECTS/DATA/BTCmeth/Data/gdac.broadinstitute.org_OV.Merge_methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.Level_3.2016012800.0.0/OV.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt"
header=read.table(methFile,nrow=2,sep="\t",as.is=T)
betaRows=grep("Beta",header[2,])
dat=as.matrix(as.data.frame(fread(methFile,header=F,sep="\t",skip=2,select=c(betaRows))))
cpgs=as.data.frame(fread(methFile,header=F,sep="\t",skip=2,select=1))

betaRows=grep("Beta",dat[2,])
colnames(dat)=header[1,betaRows]
rownames(dat)=cpgs[,1]
type=substr(colnames(dat),14,15)
rownames(dat)=dat[-(1:2),1]

datM.N=datM[,which(type=="11")]
datM.T=datM[,which(type=="01")]
colnames(datM.N)=substr(colnames(datM.N),1,12)
colnames(datM.T)=substr(colnames(datM.T),1,12)



attach("~/PROJECTS/DATA/ovarianMeth/enhancer/epicData.rda")

##############################################
### june data

library(minfi)
targets=read.metharray.sheet("Data/juneData",pattern=".csv")
for(i in unique(targets$Slide)) targets$Slide[grep(i,targets$Name.Bar.Code.Position)]=i
targets$Slide=paste0(substr(targets$Slide,1,8),"0",substr(targets$Slide,9,11))
#and slide4=sapply(strsplit(targets4$Sample_Name,"_"),geti)
targets$Basename=file.path("Data/juneData",targets$Slide,paste(targets$Slide,targets$Array,sep="_"))
rownames(targets)=paste(targets$Slide,targets$Array,sep="_")

#### and field effect data
targets2=read.metharray.sheet("Data/Field Effect",pattern=".csv")
rownames(targets2)=paste(targets2$Slide,targets2$Array,sep="_")

colnames(targets2)[1]="Name.Bar.Code.Position"
comNames=intersect(colnames(targets),colnames(targets2))
targs=rbind(targets[,comNames],targets2[,comNames])

##### process
RGset=read.metharray.exp(targ=targs,force=T)

annotation(RGset)["array"]="IlluminaHumanMethylationEPIC"
annotation(RGset)["annotation"]="ilm10b2.hg19"
getManifest(RGset)
#prob=preprocessIllumina(RGset)
probF=preprocessFunnorm(RGset)
betaJ=getBeta(probF)
betaJA=betaJ[which(!(Locations[rownames(betaJ),"chr"]%in%c("chrX","chrY"))),]


detectionP=detectionP(RGset)
pvalScores=apply(detectionP>0.01,2,mean)
pvalScores001=apply(detectionP>0.001,2,mean)


pvalCols=rep("black", length(pvalScores))
pvalCols[which(pvalScores>0.01)]="yellow"
  pvalCols[which(pvalScores>0.1)]="red"  
    
    names(pvalCols)=names(pvalScores)=rownames(targs)
    plot(density(betaJA[,1]),type="l",lwd=3,xlab="",main="",ylim=c(0,3),col=pvalCols[1])
    for(i in 2:ncol(betaJA)) lines(density(betaJA[,i]),lwd=3,col=pvalCols[i])
      
    ################ 
    ## bisulfite QC
    
    library(sesame)
    
    
    #### get probe info for these probes
    extC <- sesameDataGet(paste0("EPIC", '.probeInfo'))$typeI.extC
    extT <- sesameDataGet(paste0("EPIC", '.probeInfo'))$typeI.extT
    
    #### find relevant portions of the RGset
    ### get green matrix
    greenMat=getGreen(RGset)
    
    ### get red probes
    red1Pb=getProbeInfo(IlluminaHumanMethylationEPICmanifest, type = "I-Red")
    rownames(red1Pb)=red1Pb[,1]
    
    ### reduce red probes to those in greenmat
    red1Pb=red1Pb[which(red1Pb[,"AddressA"]%in%rownames(greenMat) & red1Pb[,"AddressB"]%in%rownames(greenMat)),]
    
    extC=intersect(extC,rownames(red1Pb))
    extT=intersect(extT,rownames(red1Pb))
    Cmn=apply(greenMat[c(red1Pb[extC,"AddressA"],red1Pb[extC,"AddressB"]),],2,mean)
    Tmn=apply(greenMat[c(red1Pb[extT,"AddressA"],red1Pb[extT,"AddressB"]),],2,mean)
    #matC=rbind(greenMat[red1Pb[extC,"AddressA"],],greenMat[red1Pb[extC,"AddressB"],]
    #save("Cmn","Tmn",file="Data/bis.Rda")

    
   
    
    
    
    ##### clinical data
    library(readxl)
    #juneClin=data.frame(read_excel("Data/062123 BRCA EPICarray information.xlsx"))
    juneClin=read.csv("Data/070223 BRCA EPICarray information.csv")[1:115,]
    juneClinS=juneClin[grep("C01",juneClin[,"iDat.file.name"]),]
    rownames(juneClinS)=juneClinS[,"iDat.file.name"]
   
    ###### which samples
    ###### targets
    betaO=betaJ[,which(!(targs[colnames(betaJ),"Name.Bar.Code.Position"]%in%rownames(juneClinS)))]
    betaJ=betaJ[,which((targs[colnames(betaJ),"Name.Bar.Code.Position"]%in%rownames(juneClinS)))]
    colJune=colnames(betaJ)
    names(colJune)=targs[colnames(betaJ),"Name.Bar.Code.Position"]
    
    ###### GEO
    sticClinID=sapply(sticClin[,"description.2"],grep,x=rownames(juneClinS),value=T)
    sticBetaO=sticBeta[,which(nchar(sticClinID)==12)]
    sticBetaJ=sticBeta[,which(nchar(sticClinID)>12)]
    sticClinID=sticClinID[which(nchar(sticClinID)>12)]
    sticClinIDv=unlist(sticClinID)
   
    sticBetaJA=sticBetaJ[rownames(betaJA),]
    betaJA=betaJA[,colnames(betaJ)]
    
    

    
    ##### missing clinical data
    missingSamps=setdiff(rownames(juneClinS),c(sticClinIDv,names(colJune)))
    ### ok all good, nothing missing
    
    ### make sure which clinical row matches with beta column
    oldNamesStic=colnames(sticBetaJ)
    colnames(sticBetaJA)=colnames(sticBetaJ)=juneClinS[sticClinIDv,"New.Sample.Names"]
    oldNames=colnames(betaJ)
    colnames(betaJA)=colnames(betaJ)=juneClinS[targs[oldNames,"Name.Bar.Code.Position"],"New.Sample.Names"]
    #### replace all ids with new names
    juneClinF=juneClinS
    rownames(juneClinF)=juneClinS[,"New.Sample.Names"]
    
    ###. check clinical data and fill in age?
    samp=sapply(strsplit(rownames(juneClinF),split=" ",fixed=T),geti,i=1)
    smp=cbind(samp,juneClinF$Age)
    smp1=smp[which(smp[,2]!=""),]
    rownames(smp1)=smp1[,1]
    ageF=rep(NA,109)
    names(ageF)=rownames(juneClinF)
    for(i in 1:nrow(smp1)){
      
      ageF[which(samp==smp1[i,1])]=smp1[i,2]
    }
    ageF=as.numeric(ageF)
    ### subselect Cmn, Tmn, pvalCols,pvalScores
    CmnFull=Cmn
    TmnFull=Tmn
    #names(pvalCols)=names(pvalScores)
    pvalColsFull=pvalCols
    pvalScoresFull=pvalScores
    
  
    pvalScores=pvalScores[oldNames]
    pvalCols=pvalCols[oldNames]
    names(pvalScores)=names(pvalCols)=names(Cmn)=names(Tmn)=colnames(betaJA)
    rownames(juneClinS)[which(juneClinS$GEO.titles!="")]=juneClinS$GEO.titles[which(juneClinS$GEO.titles!="")]
    name=juneClinS[,"Name.on.Master.file"]
    name[which(name=="")]=juneClin[which(name==""),"Sample"]
    rownames(juneClinS)=juneClinS[,"iDat.file.name"]
    sampNamesJune=juneClinS[rbind(targets[,comNames],targets2[,comNames])$Name.Bar.Code.Position,"Name.on.Master.file"]
    
    
    
    ### batch correct
    #plot(apply(sticBetaJA[,which(juneClinF[colnames(sticBetaJA),"Type.of.lesion"]=="NFTE")],1,median,na.rm=T),
    #apply(betaJA[,which(juneClinF[colnames(betaJA),"Type.of.lesion"]=="NFTE")],1,median,na.rm=T),pch=".")
# abline(0,1,col="blue")



library(sva)
    lesion=data.frame("tol"=juneClinF[c(colnames(betaJA),colnames(sticBetaJA)),"Type.of.lesion"])
mod=model.matrix(~tol,data=lesion)
cb=ComBat(cbind(betaJA,sticBetaJA),batch=rep(0:1,c(ncol(betaJA),ncol(sticBetaJA))),mod=mod)
cb[which(cb<0)]=0
cb[which(cb>1)]=1

plot(apply(cb[,colnames(sticBetaJA)][,which(juneClinF[colnames(sticBetaJA),"Type.of.lesion"]=="NFTE")],1,median,na.rm=T),
apply(cb[,colnames(betaJA)][,which(juneClinF[colnames(betaJA),"Type.of.lesion"]=="NFTE")],1,median,na.rm=T),pch=".")
abline(0,1,col="blue")


plot(apply(cb[,colnames(sticBetaJA)][,which(juneClinF[colnames(sticBetaJA),"Type.of.lesion"]=="NFTE")],1,median,na.rm=T),
     apply(cb[,colnames(sticBetaJA)][,which(juneClinF[colnames(sticBetaJA),"Type.of.lesion"]=="STIC")],1,median,na.rm=T),pch=".")
abline(0,1,col="blue")

plot(apply(cb[,colnames(sticBetaJA)][,which(juneClinF[colnames(sticBetaJA),"Type.of.lesion"]=="STIC")],1,median,na.rm=T),
     apply(cb[,colnames(betaJA)][,which(juneClinF[colnames(betaJA),"Type.of.lesion"]=="STIC")],1,median,na.rm=T),pch=".")
abline(0,1,col="blue")
beta=cb
juneClinF=juneClinF[colnames(cb),]


plot(apply(beta[,colnames(sticBetaJA)][,which(juneClinF[colnames(sticBetaJA),"Type.of.lesion"]=="NFTE")],1,median,na.rm=T),
     apply(beta[,colnames(betaJA)][,which(juneClinF[colnames(betaJA),"Type.of.lesion"]=="NFTE")],1,median,na.rm=T),pch=".")
abline(0,1,col="blue")

    save(list=c("beta","betaJA","juneClinF","pvalScores","pvalCols","Cmn","Tmn","sticClin","colJune","oldNames"),file="Data/prelimDat.Rda")
    