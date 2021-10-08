# Filtering data for population structure and phylogentic analysis 
## Continue from nsp888_minQ30_minDP8_minGQ20.vcf.gz, controlling read depth and missing data for ADMIXTURE and PCA
```
vcftools --gzvcf nsp888_minQ30_minDP8_minGQ20.vcf.gz --max-meanDP 25 --max-missing 0.9 --recode --recode-INFO-all -c | bcftools view -Oz -o nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90.vcf.gz
```
### LD pruning with plink
```
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90.vcf.gz -Oz -o nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated.vcf.gz

plink --vcf  nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated.vcf.gz --const-fid --allow-extra-chr --indep-pairwise 50 10 0.1 --out  nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune

plink --allow-extra-chr --extract nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_prune.prune.in --make-bed --out nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_prune_clean --recode vcf-iid --vcf nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_clean_annotated.vcf.gz

zcat nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune_clean.vcf.gz | sed 's/^LG//' | bgzip > nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune_clean_noLG.vcf.gz

mv nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune_clean_noLG.vcf.gz > SNPset2.vcf.gz
```
### Convert data for Plink and [ADMIXTURE](https://dalexander.github.io/admixture/download.html)
```
plink --vcf SNPset2.vcf.gz --recode --out SNPset2 --noweb

plink --vcf SNPset2.vcf.gz --make-bed --out SNPset2 --noweb 
```
## Run PCA and ADMIXTURE analysis
```
plink --bfile SNPset2 --distance 1-ibs square gz --pca --out SNPset2 --noweb

for K in `seq 2 10`; do admixture --cv SNPset2.bed $K | tee log${K}.out;  done 
```

# Phylogenomic analysis with RaxML and ASTRAL with *Pungitius.tymensis* as the outgroup (nuclear phylogeny)
## Select two individual from each popualtion and allow 50% missing data
```
bcftools view -S $downsampled.list nsp889_minQ30_minDP8_minGQ20.vcf.gz | vcftools --vcf - --max-missing 0.5 --recode --recode-INFO-all -c | bcftools view --min-ac=1 -Oz -o SNPset3.vcf.gz
```
Making even windows of 5000 SNPs each
``` 
bcftools query -f '%CHROM\t%POS\n' SNPset3.vcf.gz > positions.txt

for i in `seq 1 21`; do; cat position.txt |grep -w LG${i} >LG${i}.pos.txt ;done
```
## Extract the SNPs of each window 
### R script, subsample each window while discarding the windows in the end of chromosomes which has less than 5000 SNPs
```
for (f in c(1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20,21)){
infile<-read.table(paste0("LG",f,".pos.txt"),header = F)
for (i in 1:round(length(infile$V1)/5000)){
  a=(i-1)*5000+1
  b=i*5000
  df<-infile[a:b,]
  write.table(df, file=paste0("~/subsets/LG",f,"_region",i,".txt"), col.names = F, row.names = F, quote = F, sep = "\t")
}
}
```
### Extract the SNPs for each LG
```
vcftools --gzvcf SNPset3.vcf.gz --positions list_of_SNPs_LG#_WIN# --recode --recode-INFO-all -c | bcftools view -Oz -o ~/windowedvcf/LG#_WIN#.vcf.gz
```
Convert vcf to phylip format while defining Pungitius.tymensis as the outgroup using [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py)
```
python3 vcf2phylip.py -i ~/windowedvcf/${i}.vcf.gz -o PUN-TYM-50a
```
## Run RaxML for each of the windows
```
raxmlHPC-PTHREADS-AVX2 -T 12 -m GTRGAMMA --asc-corr=lewis -N 20 -s ${i} -n ${i}.out -p 12345
```
Merge all the best trees into one single file for ASTRAL
```
cat RAxML_bestTree*.tree > all_besttree.tree
```
Run ASTRAL with a predefined map file assgin the populations (e.g. POPA: SAMPLE1, SAPMLE2)
```
java -jar astral.5.7.4.jar -i all_besttree.tree -a astral.map -o ASTRAL_tree.speciestree 2>out1.log
```
Run ASTRAL without predefined map file and estimate the Local Posterior Probabilities for the main topology.
```
java -jar astral.5.7.4.jar -i all_besttree.tree -o all_besttree_astral.speciestree 2>out2.log
java -jar astral.5.7.4.jar -t 3 -q all_besttree_astral.speciestree -i all_besttree.tree -o all_besttree_astral_scored.speciestree 2>out3.log
```
# Mitochondria phylogeny
## Process the mitochondria data
### Fix the ploidy and apply filtering
```
bcftools plugin fixploidy -- in.vcf -s fake_sexes.txt -p mtDNA_ploidy.txt | bcftools view -Oz -o ~/mtDNA.vcf.gz
bcftools view -m2 -M2 -q 0.05 -Q 0.95 -Oz -o mtDNA_maf05.vcf.gz ~/mtDNA.vcf.gz 
```
#### mtDNA_ploidy.txt is a file assign the ploidy
```
Ppun_mitoSeq 1 16582 M 1
Ppun_mitoSeq 1 16582 F 1
```
#### fake_sexes.txt is a file specifies the gener, here we are doing this for mitochondria data so the gender does not matter
```
1-f M
10-f M
...
```
### Calculate with distance matrix with plink
```
plink \
  --allow-extra-chr
  --distance square gz 1-ibs
  --mind 0.1
  --out mtDNA
  --vcf mtDNA_maf05.vcf.gz
```
### Plot the tree in R:
```
library(ape)
library(phangorn)
nms=read.table("mtDNA.mdist.id",header=F,sep="\t")$V1
dst=read.table("mtDNA.mdist.gz",col.names=nms,header=F,fill=T,sep="\t")

njt = nj(as.dist(dst))
njt2 = root(njt, outgroup = "PUN.TYM", resolve.root = TRUE)
plot.phylo(njt2,cex=0.2,tip.color = cols)
```
### Generate the RAxML maximum likelihood tree for mitochondria data
## Convert vcf to phylip format
```
python3 vcf2phylip.py -i nsp91_mtDNA_maf01_SNPs.vcf -o PUN-TYM-50A -m 0
python3 vcf2phylip.py -i nsp889_mtDNA_maf01_SNPs.vcf -o PUN-TYM-50A -m 0
```
## Run RAxML using GTRGAMMA model, correct for ascertainment bias and run 1,000 bootstrap to infer the branch support
```
for i in $(ls *.phy)
do
echo $i
raxmlHPC-PTHREADS-AVX2 -T 12 -m GTRGAMMA --asc-corr=lewis -# 1000 -s ${i} -n ${i}.out -p 12345
done
```
## Infer the branch support
```
cat RAxML_result.nsp91_mtDNA_maf01_SNPs.min0.phy.out.RUN.* > RAxML_bestTree.nsp91_mtDNA_maf01_SNPs.min0.bootstrap.result
cat RAxML_result.nsp889_mtDNA_maf01_SNPs.min0.phy.out.RUN.* > RAxML_bestTree.nsp889_mtDNA_maf01_SNPs.min0.bootstrap.result

raxmlHPC-PTHREADS-AVX2 -T 12 -f b -m GTRGAMMA -n RAxML_bestTree.nsp91_mtDNA_maf01_SNPs.min0.bootstrap.tre -t RAxML_bestTree.nsp91_mtDNA_maf01_SNPs.min0.phy.out -z RAxML_bestTree.nsp91_mtDNA_maf01_SNPs.min0.bootstrap.result
raxmlHPC-PTHREADS-AVX2 -T 12 -f b -m GTRGAMMA -n RAxML_bestTree.nsp889_mtDNA_maf01_SNPs.min0.bootstrap.tre -t RAxML_bestTree.nsp889_mtDNA_maf01_SNPs.min0.phy.out -z RAxML_bestTree.nsp889_mtDNA_maf01_SNPs.min0.bootstrap.result
```
## R script to visulize the Cytonuclear incongruence
```
library(ape)
library(phangorn)
library(phytools)
t1 <- read.tree("ASTRAL.tre") 
t2 <- read.tree("RAxML_mtDNA.tre") 

#correct for the branch length of ASTRAL tree, if the branch lengths have 0 or NA values
t1$edge.length[is.na(t1$edge.length)] <- 0

#assign a long branch length to the outgroup (as the tree is unrooted and midpoint rooting will be used)
t1$edge.length[6] <- 20

#create an association file and specifies the sample should not be compared i.e. the outgroup, t1 and t2 must have identical sample names
assoc<-cbind(t1$tip.label,t1$tip.label) 
assoc<-assoc[-3,] #drop the outgroup sample

#compare the two trees while using midpoint rooting
obj<-cophylo(midpoint(t2),midpoint(t1),assoc=assoc)

#plot the results and annotate the nodes with support values
plot(obj,link.type="curved",link.lwd=2,
     link.col=make.transparent("blue",0.7),fsize=0.4)
nodelabels(t1$node.label,cex = .5,bg="transparent",frame="n",adj = c(1.5, -0.2))
nodelabels.cophylo(text=obj$trees[[1]]$node.label,
                   adj=c(1.5,-0.2),
                   frame="none",cex=0.5,which="left")
```
