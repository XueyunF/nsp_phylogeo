# Sex chromosomes coverage

## Convert the vcf into a table of depths.
```
bcftools view -m2 -M2 -v snps -q 0.05:minor LG12.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' | bgzip -c > LG12_nonref_depth.txt.gz

bcftools view -S finpyo.sample.list -m2 -M2 -v snps -q 0.05:minor LG12.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP]\n' > LG12_finpyo_depth.txt
```

## Then in R, find the SNPs where the pattern in FIN-PYO samples matches the expectation and use those sites for all other individuals:
```
pyo=read.table("LG12_finpyo_depth.txt",header=T,colClasses=c("factor","integer","factor","factor",rep("integer",31)),na.strings=".")

rs = rowSums(pyo[,c(8,9,10,12,13)])
kp = rs>30 & rs<80
sum(kp,na.rm=T)
for(i in c(8,9,10,12,13)){
    rs = pyo[,i]
    kp = kp & rs>=5 & rs<=15
}

kpY=kp&pyo$POS<2e7&pyo$FIN.PYO.0<50&!is.na(pyo$FIN.PYO.0)
kpA=kp&pyo$POS>2.75e7&pyo$FIN.PYO.0<85&!is.na(pyo$FIN.PYO.0)


dat=read.table("LG12_nonref_depth.txt",header=T,colClasses=c("factor","integer","factor","factor",rep("integer",918)),na.strings=".")

y.a.ratio.sites <- function(dpt){
    sumY=sum(dpt[kpY],na.rm=T)
    sumA=sum(dpt[kpA],na.rm=T)
    rt=(sumY/sum(kpY,na.rm=T))/(sumA/sum(kpA,na.rm=T))  
    nsY=sum((dpt[kpY]>5),na.rm=T)
    nsA=sum((dpt[kpA]>5),na.rm=T)
    return(c(rt,nsY,nsA))
}

res=c()
for(i in 5:922){
 rt=y.a.ratio.sites(dat[,i])
 res=rbind(res,cbind(colnames(dat)[i],t(rt)))
}

res.df=data.frame(res)
res.df$sample=as.character(res.df$sample)

res.df$pop=(substr(res.df[,1],1,7))
res.df$ratio=as.numeric(as.character(res.df$ratio))

res.df$type = "U"
res.df$type[res.df$ratio>0.8] = "O"
res.df$type[res.df$ratio<0.6&res.df$ratio>0.35] = "H"
res.df$type[as.numeric(as.character(res.df$nsA))<250000] = "U"
```
