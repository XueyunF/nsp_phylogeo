## Filtering data for population structure and phylogentic analysis 
```
# Continue from nsp888_minQ30_minDP8_minGQ20.vcf.gz, controlling read depth and missing data for ADMIXTURE and PCA

vcftools --gzvcf nsp888_minQ30_minDP8_minGQ20.vcf.gz --max-meanDP 25 --max-missing 0.9 --recode --recode-INFO-all -c | bcftools view -Oz -o nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90.vcf.gz

#LD pruning with plink

bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90.vcf.gz -Oz -o nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated.vcf.gz

plink --vcf  nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated.vcf.gz --const-fid --allow-extra-chr --indep-pairwise 50 10 0.1 --out  nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune

plink --allow-extra-chr --extract nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_prune.prune.in --make-bed --out nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_prune_clean --recode vcf-iid --vcf nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_clean_annotated.vcf.gz

zcat nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune_clean.vcf.gz | sed 's/^LG//' | bgzip > nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune_clean_noLG.vcf.gz

mv nsp888_minQ30_minDP8_minGQ20_maxmeanDP25_maxmissing90_annotated_prune_clean_noLG.vcf.gz > SNPset2.vcf.gz

# Convert data for Plink
plink --vcf SNPset2.vcf.gz --recode --out SNPset2 --noweb

plink --vcf SNPset2.vcf.gz --make-bed --out SNPset2 --noweb 

# Run PCA and ADMIXTURE analysis

plink --bfile SNPset2 --distance 1-ibs square gz --pca --out SNPset2 --noweb

for K in `seq 2 10`; do admixture --cv SNPset2.bed $K | tee log${K}.out;  done 
```

### Phylogenetic analysis with RaxML and ASTRAL with Pungitius.tymensis as the outgroup
```
# Select two individual from each popualtion and allow 50% missing data

bcftools view -S $downsampled.list nsp889_minQ30_minDP8_minGQ20.vcf.gz | vcftools --vcf - --max-missing 0.5 --recode --recode-INFO-all -c | bcftools view --min-ac=1 -Oz -o SNPset3.vcf.gz

# Making even windows of 5000 SNPs each
 
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
###Convert vcf to phylip format while defining Pungitius.tymensis as the outgroup using [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip/blob/master/vcf2phylip.py)
```
python3 vcf2phylip.py -i ~/windowedvcf/${i}.vcf.gz -o PUN-TYM-50a
```
# Run RaxML for each of the windows
```
raxmlHPC-PTHREADS-AVX2 -T 12 -m GTRGAMMA --asc-corr=lewis -N 20 -s ${i} -n ${i}.out -p 12345
```
# Put all the best trees into one sigle file for ASTRAL
```
cat RAxML_bestTree*.tree > all_besttree.tree
```
#Run ASTRAL with a predefined map file assgin the populations (e.g. POPA: SAMPLE1, SAPMLE2)
```
java -jar astral.5.7.4.jar -i all_besttree.tree -a astral.map -o ASTRAL_tree.speciestree 2>out2.log
```

# Mitochondria phylogeny


