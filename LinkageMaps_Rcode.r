setwd("E:/Data/Laupala/Rwork/")
library(qtl)
library(ggplot2)
library(gridExtra)

# read in raw genotype data
parkoh_raw<-read.delim("E:/Software/MapMaker/ParKoh_Shared_noSegDist_100k/ParKoh_Shared_noSegDist_100k_ext.raw", header=FALSE, skip=2) 
parkon_raw<-read.delim("E:/Software/MapMaker/ParKon_FB_noSegDist_100k/ParKon_FB_noSegDist100k.raw", header=FALSE, skip=2)
prukoh_raw<-read.delim("E:/Software/MapMaker/PruKoh_FB_noSegDist_100k/PruKoh_FB_noSegDist_100k.raw", header=FALSE, skip=2)

parkoh_raw$marker<-rownames(parkoh_raw)
parkon_raw$marker<-rownames(parkon_raw)
prukoh_raw$marker<-rownames(prukoh_raw)

 # read in MapMaker linkage maps at initial, comprehensive, and dense levels of density
combined_initial<-read.delim("combinedMaps_Initial.txt", header=FALSE)
combined_comprehensive<-read.delim("combinedMaps_Comprehensive.txt", header=FALSE)
combined_dense<-read.delim("combinedMaps_Dense.txt", header=FALSE)

colnames(combined_comprehensive)<-c("marker","contig","dist","cM","LG","order","cross")
colnames(combined_initial)<-c("marker","contig","dist","cM","LG","order","cross")
colnames(combined_dense)<-c("marker","contig","dist","cM","LG","order","cross")

# subset the combined maps to individual cross-specific and density-specific maps
parkoh_dense<-combined_dense[combined_dense[,"cross"]=="parkoh",] 
parkon_dense<-combined_dense[combined_dense[,"cross"]=="parkon",]
prukoh_dense<-combined_dense[combined_dense[,"cross"]=="prukoh",]
prukoh_initial<-combined_initial[combined_initial[,"cross"]=="prukoh",]
parkoh_initial<-combined_initial[combined_initial[,"cross"]=="parkoh",]
parkon_initial<-combined_initial[combined_initial[,"cross"]=="parkon",]
parkon_comprehensive<-combined_comprehensive[combined_comprehensive[,"cross"]=="parkon",]
parkoh_comprehensive<-combined_comprehensive[combined_comprehensive[,"cross"]=="parkoh",]
prukoh_comprehensive<-combined_comprehensive[combined_comprehensive[,"cross"]=="prukoh",]

# add genotype data
parkoh_dense<-merge(parkoh_dense,parkoh_raw,"marker",sort=FALSE) 
parkoh_initial<-merge(parkoh_initial,parkoh_raw,"marker",sort=FALSE)
parkoh_comprehensive<-merge(parkoh_comprehensive,parkoh_raw,"marker",sort=FALSE)

parkon_dense<-merge(parkon_dense,parkon_raw,"marker",sort=FALSE)
parkon_initial<-merge(parkon_initial,parkon_raw,"marker",sort=FALSE)
parkon_comprehensive<-merge(parkon_comprehensive,parkon_raw,"marker",sort=FALSE)

prukoh_dense<-merge(prukoh_dense,prukoh_raw,"marker",sort=FALSE)
prukoh_initial<-merge(prukoh_initial,prukoh_raw,"marker",sort=FALSE)
prukoh_comprehensive<-merge(prukoh_comprehensive,prukoh_raw,"marker",sort=FALSE)

## create segregation distortion plots using rQTL

# create and load rQTL crossfiles for comprehensive maps

parkoh_comprehensive_t<-t(parkoh_comprehensive)
write.table(parkoh_comprehensive_t,"ParKoh_Comprehensive_CrossFile.csv",sep=",", quote=FALSE, row.names=TRUE, col.names=FALSE)
# edit rows in Excel
parkoh_comprehensive_cross<-read.cross(format="csv",file="ParKoh_Comprehensive_CrossFile.csv",estimate.map=FALSE)

parkon_comprehensive_t<-t(parkon_comprehensive)
write.table(parkon_comprehensive_t,"ParKon_Comprehensive_CrossFile.csv",sep=",", quote=FALSE, row.names=TRUE, col.names=FALSE)
# edit rows in Excel
parkon_comprehensive_cross<-read.cross(format="csv",file="ParKon_Comprehensive_CrossFile.csv",estimate.map=FALSE)

prukoh_comprehensive_t<-t(prukoh_comprehensive)
write.table(prukoh_comprehensive_t,"PruKoh_Comprehensive_CrossFile.csv",sep=",", quote=FALSE, row.names=TRUE, col.names=FALSE)
# edit rows in Excel
prukoh_comprehensive_cross<-read.cross(format="csv",file="PruKoh_Comprehensive_CrossFile.csv",estimate.map=FALSE)

# create plots for segregation distortion

parkoh_comprehensive_genos <- geno.table(parkoh_comprehensive_cross, scanone.output=TRUE)
par(mfrow=c(2,1))
#plot(parkoh_comprehensive_genos, ylab=expression(paste(-log[10], " P-value")))
# sliding window: size =10, step = 2
neglog10P_slidewindow<-c()
for(j in 1:length(levels(parkoh_comprehensive_genos$chr))) {
	temp.df<-parkoh_comprehensive_genos[parkoh_comprehensive_genos[,"chr"]==as.character(levels(parkoh_comprehensive_genos$chr)[j]),]
	for(i in seq(2,nrow(temp.df),2)) { 
		if(i==2) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-1):i]),2))}
		if(i==4) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-3):i]),2))}
		if(i==6) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-5):i]),2))}
		if(i==8) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-7):i]),2))}
		if(i>=10) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-9):i]),2))}
		}
	neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i+1):nrow(temp.df)]),nrow(temp.df)-i))
	}
parkoh_comprehensive_genos$neglog10P<-neglog10P_slidewindow
plot(parkoh_comprehensive_genos, ylab=expression(paste(-log[10], " P-value")))
plot(parkoh_comprehensive_genos, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")

parkon_comprehensive_genos <- geno.table(parkon_comprehensive_cross, scanone.output=TRUE)
par(mfrow=c(2,1))
#plot(parkon_comprehensive_genos, ylab=expression(paste(-log[10], " P-value")))
# sliding window: size =10, step = 2
neglog10P_slidewindow<-c()
for(j in 1:length(levels(parkon_comprehensive_genos$chr))) {
	temp.df<-parkon_comprehensive_genos[parkon_comprehensive_genos[,"chr"]==as.character(levels(parkon_comprehensive_genos$chr)[j]),]
	for(i in seq(2,nrow(temp.df),2)) { 
		if(i==2) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-1):i]),2))}
		if(i==4) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-3):i]),2))}
		if(i==6) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-5):i]),2))}
		if(i==8) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-7):i]),2))}
		if(i>=10) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-9):i]),2))}
		}
	neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i+1):nrow(temp.df)]),nrow(temp.df)-i))
	}
parkon_comprehensive_genos$neglog10P<-neglog10P_slidewindow
plot(parkon_comprehensive_genos, ylab=expression(paste(-log[10], " P-value")))
plot(parkon_comprehensive_genos, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")


prukoh_comprehensive_genos <- geno.table(prukoh_comprehensive_cross, scanone.output=TRUE)
par(mfrow=c(2,1))
#plot(prukoh_comprehensive_genos, ylab=expression(paste(-log[10], " P-value")))
# sliding window: size =10, step = 2
neglog10P_slidewindow<-c()
for(j in 1:length(levels(prukoh_comprehensive_genos$chr))) {
	temp.df<-prukoh_comprehensive_genos[prukoh_comprehensive_genos[,"chr"]==as.character(levels(prukoh_comprehensive_genos$chr)[j]),]
	for(i in seq(2,nrow(temp.df),2)) { 
		if(i==2) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-1):i]),2))}
		if(i==4) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-3):i]),2))}
		if(i==6) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-5):i]),2))}
		if(i==8) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-7):i]),2))}
		if(i>=10) { neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i-9):i]),2))}
		}
	neglog10P_slidewindow<-c(neglog10P_slidewindow,rep(mean(temp.df$neglog10P[(i+1):nrow(temp.df)]),nrow(temp.df)-i))
	}
prukoh_comprehensive_genos$neglog10P<-neglog10P_slidewindow
plot(prukoh_comprehensive_genos, ylab=expression(paste(-log[10], " P-value")))
plot(prukoh_comprehensive_genos, lod=3:5, ylab="Genotype frequency")
abline(h=c(0.25, 0.5), lty=2, col="gray")


## Spearman's rank test for collinearity of maps 

rbind(c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="1","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="1","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="1","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="1","cMparkon"], method="spearman")$p.value),
c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="2","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="2","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="2","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="2","cMparkon"], method="spearman")$p.value),
c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="3","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="3","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="3","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="3","cMparkon"], method="spearman")$p.value),
c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="4","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="4","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="4","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="4","cMparkon"], method="spearman")$p.value),
c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="5","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="5","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="5","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="5","cMparkon"], method="spearman")$p.value),
c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="6","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="6","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="6","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="6","cMparkon"], method="spearman")$p.value),
c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="7","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="7","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="7","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="7","cMparkon"], method="spearman")$p.value),
c(cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="X","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="X","cMparkon"], method="spearman")$estimate,
cor.test(parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="X","cMparkoh"] , parkoh_parkon_cM_dense[parkoh_parkon_cM_dense[,"LGparkoh"]=="X","cMparkon"], method="spearman")$p.value))

rbind(c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="1","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="1","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="1","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="1","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="2","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="2","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="2","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="2","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="3","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="3","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="3","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="3","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="4","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="4","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="4","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="4","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="5","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="5","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="5","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="5","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="6","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="6","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="6","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="6","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="7","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="7","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="7","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="7","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="X","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="X","cMprukoh"], method="spearman")$estimate,
cor.test(parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="X","cMparkoh"] , parkoh_prukoh_cM_dense[parkoh_prukoh_cM_dense[,"LGparkoh"]=="X","cMprukoh"], method="spearman")$p.value))


rbind(c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="1","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="1","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="1","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="1","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="2","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="2","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="2","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="2","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="3","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="3","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="3","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="3","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="4","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="4","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="4","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="4","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="5","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="5","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="5","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="5","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="6","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="6","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="6","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="6","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="7","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="7","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="7","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="7","cMprukoh"], method="spearman")$p.value),
c(cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="X","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="X","cMprukoh"], method="spearman")$estimate,
cor.test(parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="X","cMparkon"] , parkon_prukoh_cM_dense[parkon_prukoh_cM_dense[,"LGparkon"]=="X","cMprukoh"], method="spearman")$p.value))

# export csv to ALLMAPS for map integration (pseudomolecule assembly)

parkoh_allmaps<-parkoh_dense[,c("V1","LG","cM")]
parkoh_allmaps<-cbind(data.frame(within(parkoh_allmaps, Scaffold<-data.frame(do.call('rbind', strsplit(as.character(V1), '_', fixed=TRUE))))[,4]),parkoh_allmaps[,c(2:3)])
colnames(parkoh_allmaps)<-c("contig","pos","LG","cM")
parkoh_allmaps$contig<-gsub("XS","S",parkoh_allmaps$contig)
write.table(parkoh_allmaps,"ParKoh_ALLMAPS.csv", sep=",", quote=FALSE, row.names=FALSE)
parkoh_allmaps$contig<-factor(as.character(parkoh_allmaps$contig))
parkoh_allmaps$pos<-as.numeric(as.character(parkoh_allmaps$pos))

parkon_allmaps<-parkon_dense[,c("V1","LG","cM")]
parkon_allmaps<-cbind(data.frame(within(parkon_allmaps, Scaffold<-data.frame(do.call('rbind', strsplit(as.character(V1), '_', fixed=TRUE))))[,4]),parkon_allmaps[,c(2:3)])
colnames(parkon_allmaps)<-c("contig","pos","LG","cM")
write.table(parkon_allmaps,"ParKon_ALLMAPS.csv", sep=",", quote=FALSE, row.names=FALSE)
parkon_allmaps$contig<-factor(as.character(parkon_allmaps$contig))
parkon_allmaps$pos<-as.numeric(as.character(parkon_allmaps$pos))


prukoh_allmaps<-prukoh_dense[,c("V1","LG","cM")]
prukoh_allmaps<-cbind(data.frame(within(prukoh_allmaps, Scaffold<-data.frame(do.call('rbind', strsplit(as.character(V1), '_', fixed=TRUE))))[,4]),prukoh_allmaps[,c(2:3)])
colnames(prukoh_allmaps)<-c("contig","pos","LG","cM")
write.table(prukoh_allmaps,"PruKoh_ALLMAPS.csv", sep=",", quote=FALSE, row.names=FALSE)
prukoh_allmaps$contig<-factor(as.character(prukoh_allmaps$contig))
prukoh_allmaps$pos<-as.numeric(as.character(prukoh_allmaps$pos))

## Calculate recombination rates and make Marey maps + recombination rate traces

genome<-read.delim("combinedMaps_newLG4.txt", header=FALSE) # this is a slightly modified chain output file from ALLMAPs (pseudomolucule assembly agp file)
colnames(genome)<-c("contig","length","direction","0","length2","LG","number","direction2","start","end","order")
genome$contig<-gsub("Lko057","",genome$contig)
genome$midpoint<-(genome$start+genome$end)/2

# restructure output depending on scaffold orientation
genome_new<-data.frame(matrix(nrow=1,ncol=ncol(genome)))
colnames(genome_new)<-colnames(genome)
for( i in 1:8) { 
	genome_sub<-genome[genome[,"LG"]==levels(genome$LG)[i],]; genome_sub$start=0; genome_sub$end=0
	for(j in 1:nrow(genome_sub)) {	
		if(j==1 && genome_sub$direction2[j]=="+") { genome_sub$start[j]=0; genome_sub$end[j]=genome_sub$length[j] }
		if(j==1 && genome_sub$direction2[j]=="-") { genome_sub$start[j]=genome_sub$length[j]; genome_sub$end[j]=0 }
		if(j > 1 && genome_sub$direction2[j]=="+") { genome_sub$start[j]=sum(genome_sub$length[1:(j-1)]); genome_sub$end[j]=genome_sub$start[j]+genome_sub$length[j]}
		if(j > 1 && genome_sub$direction2[j]=="-") { genome_sub$end[j]=sum(genome_sub$length[1:(j-1)]); genome_sub$start[j]=genome_sub$end[j]+genome_sub$length[j] }
		}
	genome_new<-rbind(genome_new,genome_sub)
	}
genome_new<-genome_new[-1,]
genome_new$midpoint<-(genome_new$start+genome_new$end)/2


# align pseudomolecule assembly with species-pair specific maps

parkoh_allmaps2<-merge(parkoh_allmaps,genome_new[,c("contig","midpoint")],"contig",sort=FALSE)
parkon_allmaps2<-merge(parkon_allmaps,genome_new[,c("contig","midpoint")],"contig",sort=FALSE)
prukoh_allmaps2<-merge(prukoh_allmaps,genome_new[,c("contig","midpoint")],"contig",sort=FALSE)

# calculate splines for physical position ~ genetic position  and first derivative (recombination rate)
# then plot genetic position against physical position, fitted splines, and first derivative of splines.

spline_parkoh_LG1<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="1","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="1","cM"], df=10)
deriv_spline_parkoh_LG1<-predict(spline_parkoh_LG1,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="1","midpoint"],deriv = 1)
deriv_spline_parkoh_LG1<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LG1$x),"cM"=as.numeric(deriv_spline_parkoh_LG1$y))
deriv_spline_parkoh_LG1<-deriv_spline_parkoh_LG1[order(deriv_spline_parkoh_LG1$pos),]

spline_parkon_LG1<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="1","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="1","cM"], df=10)
deriv_spline_parkon_LG1<-predict(spline_parkon_LG1,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="1","midpoint"],deriv = 1)
deriv_spline_parkon_LG1<-data.frame("pos"=as.numeric(deriv_spline_parkon_LG1$x),"cM"=as.numeric(deriv_spline_parkon_LG1$y))
deriv_spline_parkon_LG1<-deriv_spline_parkon_LG1[order(deriv_spline_parkon_LG1$pos),]

spline_prukoh_LG1<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="1","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="1","cM"], df=10)
deriv_spline_prukoh_LG1<-predict(spline_prukoh_LG1,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="1","midpoint"],deriv = 1)
deriv_spline_prukoh_LG1<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LG1$x),"cM"=as.numeric(deriv_spline_prukoh_LG1$y))
deriv_spline_prukoh_LG1<-deriv_spline_prukoh_LG1[order(deriv_spline_prukoh_LG1$pos),]

recombination_trace_LG1<-function() { par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="1","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="1","cM"], col="black", xlim=c(0,1.4e8), ylim=c(0,220),pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="1","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="1","cM"], col="grey70", xlim=c(0,1.4e8), ylim=c(0,220), pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="1","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="1","cM"], col="grey20", xlim=c(0,1.4e8), ylim=c(0,220), pch=5)
lines(spline_parkoh_LG1, col="black", xlim=c(0,1.4e8), ylim=c(0,220), lty=1)
lines(spline_parkon_LG1, col="grey70", xlim=c(0,1.4e8), ylim=c(0,220), lty=2)
lines(spline_prukoh_LG1, col="grey20", xlim=c(0,1.4e8), ylim=c(0,220), lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LG1,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.4e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LG1,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.4e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LG1,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.4e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

spline_parkoh_LG2<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="2","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="2","cM"], df=10)
deriv_spline_parkoh_LG2<-predict(spline_parkoh_LG2,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="2","midpoint"],deriv = 1)
deriv_spline_parkoh_LG2<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LG2$x),"cM"=as.numeric(deriv_spline_parkoh_LG2$y))
deriv_spline_parkoh_LG2<-deriv_spline_parkoh_LG2[order(deriv_spline_parkoh_LG2$pos),]

spline_parkon_LG2<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="2","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="2","cM"], df=10)
deriv_spline_parkon_LG2<-predict(spline_parkon_LG2,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="2","midpoint"],deriv = 1)
deriv_spline_parkon_LG2<-data.frame("pos"=as.numeric(deriv_spline_parkon_LG2$x),"cM"=as.numeric(deriv_spline_parkon_LG2$y))
deriv_spline_parkon_LG2<-deriv_spline_parkon_LG2[order(deriv_spline_parkon_LG2$pos),]

spline_prukoh_LG2<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="2","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="2","cM"], df=10)
deriv_spline_prukoh_LG2<-predict(spline_prukoh_LG2,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="2","midpoint"],deriv = 1)
deriv_spline_prukoh_LG2<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LG2$x),"cM"=as.numeric(deriv_spline_prukoh_LG2$y))
deriv_spline_prukoh_LG2<-deriv_spline_prukoh_LG2[order(deriv_spline_prukoh_LG2$pos),]

recombination_trace_LG2<-function() {
par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="2","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="2","cM"], col="black", xlim=c(0,1.2e8), ylim=c(0,220),  pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="2","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="2","cM"], col="grey70", xlim=c(0,1.2e8), ylim=c(0,220),  pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="2","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="2","cM"], col="grey20", xlim=c(0,1.2e8), ylim=c(0,220),  pch=5)
lines(spline_parkoh_LG2, col="black", lty=1)
lines(spline_parkon_LG2, col="grey70", lty=2)
lines(spline_prukoh_LG2, col="grey20", lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LG2,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.2e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LG2,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.2e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LG2,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.2e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

spline_parkoh_LG3<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="3","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="3","cM"], df=10)
deriv_spline_parkoh_LG3<-predict(spline_parkoh_LG3,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="3","midpoint"],deriv = 1)
deriv_spline_parkoh_LG3<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LG3$x),"cM"=as.numeric(deriv_spline_parkoh_LG3$y))
deriv_spline_parkoh_LG3<-deriv_spline_parkoh_LG3[order(deriv_spline_parkoh_LG3$pos),]

spline_parkon_LG3<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="3","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="3","cM"], df=10)
deriv_spline_parkon_LG3<-predict(spline_parkon_LG3,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="3","midpoint"],deriv = 1)
deriv_spline_parkon_LG3<-data.frame("pos"=as.numeric(deriv_spline_parkon_LG3$x),"cM"=as.numeric(deriv_spline_parkon_LG3$y))
deriv_spline_parkon_LG3<-deriv_spline_parkon_LG3[order(deriv_spline_parkon_LG3$pos),]

spline_prukoh_LG3<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="3","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="3","cM"], df=10)
deriv_spline_prukoh_LG3<-predict(spline_prukoh_LG3,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="3","midpoint"],deriv = 1)
deriv_spline_prukoh_LG3<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LG3$x),"cM"=as.numeric(deriv_spline_prukoh_LG3$y))
deriv_spline_prukoh_LG3<-deriv_spline_prukoh_LG3[order(deriv_spline_prukoh_LG3$pos),]

recombination_trace_LG3<- function() {
par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="3","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="3","cM"], col="black", xlim=c(0,1.1e8), ylim=c(0,220),  pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="3","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="3","cM"], col="grey70", xlim=c(0,1.1e8), ylim=c(0,220),  pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="3","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="3","cM"], col="grey20", xlim=c(0,1.1e8), ylim=c(0,220),  pch=5)
lines(spline_parkoh_LG3, col="black", lty=1)
lines(spline_parkon_LG3, col="grey70", lty=2)
lines(spline_prukoh_LG3, col="grey20", lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LG3,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.1e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LG3,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.1e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LG3,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.1e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

spline_parkoh_LG4<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="4","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="4","cM"], df=10)
deriv_spline_parkoh_LG4<-predict(spline_parkoh_LG4,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="4","midpoint"],deriv = 1)
deriv_spline_parkoh_LG4<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LG4$x),"cM"=as.numeric(deriv_spline_parkoh_LG4$y))
deriv_spline_parkoh_LG4<-deriv_spline_parkoh_LG4[order(deriv_spline_parkoh_LG4$pos),]

spline_parkon_LG4<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="4","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="4","cM"], df=10)
deriv_spline_parkon_LG4<-predict(spline_parkon_LG4,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="4","midpoint"],deriv = 1)
deriv_spline_parkon_LG4<-data.frame("pos"=as.numeric(deriv_spline_parkon_LG4$x),"cM"=as.numeric(deriv_spline_parkon_LG4$y))
deriv_spline_parkon_LG4<-deriv_spline_parkon_LG4[order(deriv_spline_parkon_LG4$pos),]

spline_prukoh_LG4<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="4","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="4","cM"], df=10)
deriv_spline_prukoh_LG4<-predict(spline_prukoh_LG4,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="4","midpoint"],deriv = 1)
deriv_spline_prukoh_LG4<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LG4$x),"cM"=as.numeric(deriv_spline_prukoh_LG4$y))
deriv_spline_prukoh_LG4<-deriv_spline_prukoh_LG4[order(deriv_spline_prukoh_LG4$pos),]

recombination_trace_LG4<- function () {
par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="4","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="4","cM"], col="black", xlim=c(0,.9e8), ylim=c(0,220),  pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="4","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="4","cM"], col="grey70", xlim=c(0,.9e8), ylim=c(0,220),  pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="4","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="4","cM"], col="grey20", xlim=c(0,.9e8), ylim=c(0,220),  pch=5)
lines(spline_parkoh_LG4, col="black", lty=1)
lines(spline_parkon_LG4, col="grey70", lty=2)
lines(spline_prukoh_LG4, col="grey20", lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LG4,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,.9e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LG4,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,.9e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LG4,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,.9e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

spline_parkoh_LG5<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="5","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="5","cM"], df=10)
deriv_spline_parkoh_LG5<-predict(spline_parkoh_LG5,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="5","midpoint"],deriv = 1)
deriv_spline_parkoh_LG5<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LG5$x),"cM"=as.numeric(deriv_spline_parkoh_LG5$y))
deriv_spline_parkoh_LG5<-deriv_spline_parkoh_LG5[order(deriv_spline_parkoh_LG5$pos),]

spline_parkon_LG5<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="5","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="5","cM"], df=10)
deriv_spline_parkon_LG5<-predict(spline_parkon_LG5,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="5","midpoint"],deriv = 1)
deriv_spline_parkon_LG5<-data.frame("pos"=as.numeric(deriv_spline_parkon_LG5$x),"cM"=as.numeric(deriv_spline_parkon_LG5$y))
deriv_spline_parkon_LG5<-deriv_spline_parkon_LG5[order(deriv_spline_parkon_LG5$pos),]

spline_prukoh_LG5<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="5","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="5","cM"], df=10)
deriv_spline_prukoh_LG5<-predict(spline_prukoh_LG5,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="5","midpoint"],deriv = 1)
deriv_spline_prukoh_LG5<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LG5$x),"cM"=as.numeric(deriv_spline_prukoh_LG5$y))
deriv_spline_prukoh_LG5<-deriv_spline_prukoh_LG5[order(deriv_spline_prukoh_LG5$pos),]

recombination_trace_LG5<- function() {
par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="5","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="5","cM"], col="black", xlim=c(0,0.7e8), ylim=c(0,220),  pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="5","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="5","cM"], col="grey70", xlim=c(0,0.7e8), ylim=c(0,220),  pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="5","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="5","cM"], col="grey20", xlim=c(0,0.7e8), ylim=c(0,220),  pch=5)
lines(spline_parkoh_LG5, col="black", lty=1)
lines(spline_parkon_LG5, col="grey70", lty=2)
lines(spline_prukoh_LG5, col="grey20", lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LG5,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.7e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LG5,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.7e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LG5,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.7e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

spline_parkoh_LG6<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="6","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="6","cM"], df=10)
deriv_spline_parkoh_LG6<-predict(spline_parkoh_LG6,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="6","midpoint"],deriv = 1)
deriv_spline_parkoh_LG6<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LG6$x),"cM"=as.numeric(deriv_spline_parkoh_LG6$y))
deriv_spline_parkoh_LG6<-deriv_spline_parkoh_LG6[order(deriv_spline_parkoh_LG6$pos),]

spline_parkon_LG6<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="6","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="6","cM"], df=10)
deriv_spline_parkon_LG6<-predict(spline_parkon_LG6,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="6","midpoint"],deriv = 1)
deriv_spline_parkon_LG6<-data.frame("pos"=as.numeric(deriv_spline_parkon_LG6$x),"cM"=as.numeric(deriv_spline_parkon_LG6$y))
deriv_spline_parkon_LG6<-deriv_spline_parkon_LG6[order(deriv_spline_parkon_LG6$pos),]

spline_prukoh_LG6<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="6","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="6","cM"], df=10)
deriv_spline_prukoh_LG6<-predict(spline_prukoh_LG6,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="6","midpoint"],deriv = 1)
deriv_spline_prukoh_LG6<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LG6$x),"cM"=as.numeric(deriv_spline_prukoh_LG6$y))
deriv_spline_prukoh_LG6<-deriv_spline_prukoh_LG6[order(deriv_spline_prukoh_LG6$pos),]

recombination_trace_LG6<-function() {
par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="6","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="6","cM"], col="black", xlim=c(0,0.6e8), ylim=c(0,220),  pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="6","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="6","cM"], col="grey70", xlim=c(0,0.6e8), ylim=c(0,220),  pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="6","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="6","cM"], col="grey20", xlim=c(0,0.6e8), ylim=c(0,220),  pch=5)
lines(spline_parkoh_LG6, col="black", lty=1)
lines(spline_parkon_LG6, col="grey70", lty=2)
lines(spline_prukoh_LG6, col="grey20", lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LG6,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.6e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LG6,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.6e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LG6,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.6e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

spline_parkoh_LG7<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="7","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="7","cM"], df=5)
deriv_spline_parkoh_LG7<-predict(spline_parkoh_LG7,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="7","midpoint"],deriv = 1)
deriv_spline_parkoh_LG7<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LG7$x),"cM"=as.numeric(deriv_spline_parkoh_LG7$y))
deriv_spline_parkoh_LG7<-deriv_spline_parkoh_LG7[order(deriv_spline_parkoh_LG7$pos),]

spline_parkon_LG7<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="7","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="7","cM"], df=5)
deriv_spline_parkon_LG7<-predict(spline_parkon_LG7,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="7","midpoint"],deriv = 1)
deriv_spline_parkon_LG7<-data.frame("pos"=as.numeric(deriv_spline_parkon_LG7$x),"cM"=as.numeric(deriv_spline_parkon_LG7$y))
deriv_spline_parkon_LG7<-deriv_spline_parkon_LG7[order(deriv_spline_parkon_LG7$pos),]

spline_prukoh_LG7<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="7","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="7","cM"], df=5)
deriv_spline_prukoh_LG7<-predict(spline_prukoh_LG7,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="7","midpoint"],deriv = 1)
deriv_spline_prukoh_LG7<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LG7$x),"cM"=as.numeric(deriv_spline_prukoh_LG7$y))
deriv_spline_prukoh_LG7<-deriv_spline_prukoh_LG7[order(deriv_spline_prukoh_LG7$pos),]

recombination_trace_LG7<- function() {
par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="7","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="7","cM"], col="black", xlim=c(0,0.3e8), ylim=c(0,220),  pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="7","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="7","cM"], col="grey70", xlim=c(0,0.3e8), ylim=c(0,220),  pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="7","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="7","cM"], col="grey20", xlim=c(0,0.3e8), ylim=c(0,220),  pch=5)
lines(spline_parkoh_LG7, col="black", lty=1)
lines(spline_parkon_LG7, col="grey70", lty=2)
lines(spline_prukoh_LG7, col="grey20", lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LG7,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.3e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LG7,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.3e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LG7,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,0.3e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

spline_parkoh_LGX<-smooth.spline(x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="X","midpoint"],y=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="X","cM"], df=10)
deriv_spline_parkoh_LGX<-predict(spline_parkoh_LGX,x=parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="X","midpoint"],deriv = 1)
deriv_spline_parkoh_LGX<-data.frame("pos"=as.numeric(deriv_spline_parkoh_LGX$x),"cM"=as.numeric(deriv_spline_parkoh_LGX$y))
deriv_spline_parkoh_LGX<-deriv_spline_parkoh_LGX[order(deriv_spline_parkoh_LGX$pos),]

spline_parkon_LGX<-smooth.spline(x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="X","midpoint"],y=parkon_allmaps2[parkon_allmaps2[,"LG"]=="X","cM"], df=10)
deriv_spline_parkon_LGX<-predict(spline_parkon_LGX,x=parkon_allmaps2[parkon_allmaps2[,"LG"]=="X","midpoint"],deriv = 1)
deriv_spline_parkon_LGX<-data.frame("pos"=as.numeric(deriv_spline_parkon_LGX$x),"cM"=as.numeric(deriv_spline_parkon_LGX$y))
deriv_spline_parkon_LGX<-deriv_spline_parkon_LGX[order(deriv_spline_parkon_LGX$pos),]

spline_prukoh_LGX<-smooth.spline(x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="X","midpoint"],y=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="X","cM"], df=10)
deriv_spline_prukoh_LGX<-predict(spline_prukoh_LGX,x=prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="X","midpoint"],deriv = 1)
deriv_spline_prukoh_LGX<-data.frame("pos"=as.numeric(deriv_spline_prukoh_LGX$x),"cM"=as.numeric(deriv_spline_prukoh_LGX$y))
deriv_spline_prukoh_LGX<-deriv_spline_prukoh_LGX[order(deriv_spline_prukoh_LGX$pos),]

recombination_trace_LGX<- function() { 
par(mar=c(5,4,4,5)+.1)
plot(parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="X","midpoint"],parkoh_allmaps2[parkoh_allmaps2[,"LG"]=="X","cM"], col="black", xlim=c(0,1.5e8), ylim=c(0,220),  pch=1)
points(parkon_allmaps2[parkon_allmaps2[,"LG"]=="X","midpoint"],parkon_allmaps2[parkon_allmaps2[,"LG"]=="X","cM"], col="grey70", xlim=c(0,1.5e8), ylim=c(0,220),  pch=2)
points(prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="X","midpoint"],prukoh_allmaps2[prukoh_allmaps2[,"LG"]=="X","cM"], col="grey20", xlim=c(0,1.5e8), ylim=c(0,220),  pch=5)
lines(spline_parkoh_LGX, col="black", lty=1)
lines(spline_parkon_LGX, col="grey70", lty=2)
lines(spline_prukoh_LGX, col="grey20", lty=3)
par(new=TRUE)
plot(deriv_spline_parkoh_LGX,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.5e8), ylim=c(0,8e-6), lty=1)
lines(deriv_spline_parkon_LGX,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.5e8), ylim=c(0,8e-6), lty=2)
lines(deriv_spline_prukoh_LGX,type="l",xaxt="n",yaxt="n",xlab="",ylab="", col="darkred", xlim=c(0,1.5e8), ylim=c(0,8e-6), lty=3)
axis(4)
}

cairo_pdf("recombination_trace.pdf",width=9, height=18)
par(mfrow=c(4,2))
recombination_trace_LG2()
recombination_trace_LG3()
recombination_trace_LG1()
recombination_trace_LG4()
recombination_trace_LG5()
recombination_trace_LG7()
recombination_trace_LG6()
recombination_trace_LGX()
dev.off()
