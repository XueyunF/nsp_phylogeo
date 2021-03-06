# MSMC2 cross population analysis
## Prepare data for analysis
### Genetic map
```cat genetic_map_nosexchr.txt |awk -F "\t" 'BEGIN {OFS = FS}  {print $2,$1,$3}'  > map.reformat.txt
for i in `seq 1 21`; do ; cat map.reformat.txt |grep -w LG$i > LG$i.map; done
```
## Calculate the mean depth based on the largest linkeage group LG4
```
for i in $(ls ~/${pop}/*_realgn_sorted.bam); do ; echo $i ; samtools depth -r LG4 $i | awk '{sum += $3} END {print sum / NR}' ; done > ${pop}.DEPTH
```

## Variant calling
```
REF=~/NSP_V6b.fasta
cat ${pop}.DEPTH | while read -r a b 
do echo $a $b 
for CHR in {1..21}; do ; 
samtools mpileup -B -q 20 -Q 20 -C 50 -g -r LG$CHR -f $REF ~/${sample}_realgn_sorted.bam | bcftools call -c -V indels |python3.5 ~/bamCaller.py ${b} ~/${pop}/${a}.mask.LG$CHR.bed.gz |gzip -c > ~/$pop/${a}.LG$CHR.vcf.gz  ; done ; done
```

## Phase the vcf 
### Keep biallelic SNP only and merge the vcfs for each LG 
```
bcftools view -M 2 -O z $individual.LG${CHR}.vcf.gz > ~/unphasevcf/$individual.LG${CHR}.unphase.vcf.gz
ls ~/unphasevcf/*.LG${CHR}.unphase.vcf.gz > $pop.LG${CHR}.list
bcftools merge -l $pop.LG${CHR}.list -Oz -o $pop.LG${CHR}.unphase.vcf.gz 
```
### Phase the data with [SHAPEIT](https://odelaneau.github.io/shapeit4/)
```
tabix -p vcf $pop.LG${CHR}.unphase.vcf.gz
shapeit4.2 --input $pop.LG${CHR}.unphase.vcf.gz --map ~/LG${CHR}.map --region LG${CHR} --output $pop.LG${CHR}.phased.vcf.gz
```
## Genearte input data for MSMC2 (2 samples from each population)
```
bcftools view -s popA_sample1 popA.LG${CHR}.phased.vcf.gz -Oz -o popA_sample1.LG${CHR}.phased.vcf.gz
bcftools view -s popA_sample2 popA.LG${CHR}.phased.vcf.gz -Oz -o popA_sample2.LG${CHR}.phased.vcf.gz
bcftools view -s popB_sample1 popB.LG${CHR}.phased.vcf.gz -Oz -o popB_sample1.LG${CHR}.phased.vcf.gz
bcftools view -s popB_sample2 popB.LG${CHR}.phased.vcf.gz -Oz -o popB_sample2.LG${CHR}.phased.vcf.gz

python3.5 ~/generate_multihetsep.py  --chr LG${CHR} --negative_mask=~/mask/mask_LG${CHR}.bed.gz \
--mask=~/popA_sample1.mask.LG${CHR}.bed.gz --mask=~/popA_sample2.mask.LG${CHR}.bed.gz \
--mask=~/popB_sample1.mask.LG${CHR}.bed.gz --mask=~/popB_sample2.mask.LG${CHR}.bed.gz \
~/popA_sample1.LG${CHR}.phased.vcf.gz ~/popA_sample2.LG${CHR}.phased.vcf.gz \
~/popB_sample1.LG${CHR}.phased.vcf.gz  ~/popB_sample2.LG${CHR}.phased.vcf.gz  > ~/popA_popB.LG${CHR}.multihetsep.txt
```
## Run MSMC2 for each population and then for cross population analysis
```
msmc2 -s -t 15 -I 0,1 -o popA  ~/popA_popB.*.multihetsep.txt
msmc2 -s -t 15 -I 2,3 -o popB  ~/popA_popB.*.multihetsep.txt
msmc2 -s -t 15 -I 0-2,0-3,1-2,1-3 -o popA_popB  ~/popA_popB.*.multihetsep.txt
python3.5 ~/combineCrossCoal.py  popA_popB.final.txt popA.final.txt popB.final.txt > popA_popB.combined.msmc2.final.txt
```
## Plot the results in R 
```
mut <- 1.42e-8
gener <- 2
crossPopDat<-read.table("~/popA_popB.combined.msmc2.final.txt", header=TRUE)

plot(crossPopDat$left_time_boundary/mut*gener, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11),
     xlim=c(0,4000000),ylim=c(0,1), type="n", xlab="Years ago", ylab="Relative Cross-Coalescence Rate",main="Split time between populations")
lines(crossPopDat$left_time_boundary/mut*gener, 2 * crossPopDat$lambda_01 / (crossPopDat$lambda_00 + crossPopDat$lambda_11), type="s",lwd=2)
abline(h=0.5,lty=2,lwd=2)
legend(300000, 0.45, legend=c("popA v.s. popB"),
       col=c("orange"), lty=1, cex=0.5)
```

