# Reads alignment

## Process the reference and do alignments
### Index the reference
```
samtools dict ninespine_v6b.fa > ninespine_v6b.dict
samtools faidx ninespine_v6b.fa
bwa index ninespine_v6b.fa
```
### Align the raw reads
```
bwa mem -t4 -M -R "@RG\tID:SampleID tSM:SampleID tPL:illumina tLB:SampleID tPU:1" ninespine.fa SampleID_1.fq.gz SampleID_2.fq.gz | samtools view -F4 -h -Obam -o sample1_bwa.bam
```
### Sort reads by coordinate, duplicates remove and index
```
samtools fixmate -@ 14 -m SampleID_bwa.bam SampleID_bwa_fm.bam
samtools sort -@ 14 -O bam -o SampleID_bwa_fm_sorted.bam SampleID_bwa_fm.bam
samtools index -@ 14 SampleID_bwa_fm_sorted.bam SampleID_bwa_fm_sorted.bai
samtools markdup -@ 14 SampleID_bwa_fm_sorted.bam SampleID_bwa_fm_sorted_md.bam
samtools index -@ 14 SampleID_bwa_fm_sorted_md.bam SampleID_bwa_fm_sorted_md.bai
```
### Realignment around indels
```
gatk -T RealignerTargetCreator -R ninespine_v6b.fa -I SampleID_bwa_fm_sorted_md.bam -o SampleID_intervals.list

gatk -T IndelRealigner -R ninespine_v6b.fa -I SampleID_bwa_fm_sorted_md.bam -targetIntervals SampleID_intervals.list -o SampleID_realgn.bam
```
### Sort and index of cleaned bam
```
samtools sort -T /tmp/SampleID -O bam -o SampleID_realgn_sorted.bam SampleID_realgn.bam
samtools index SampleID_realgn_sorted.bam SampleID_realgn_sorted.bai
```

# Variant calling
## Individual call
```
gatk -T HaplotypeCaller -R ninespine_v6b.fa -I SampleID_realgn_sorted.bam -gt_mode DISCOVERY -ERC GVCF -stand_emit_conf 3 -stand_call_conf 10 -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQ 50 -variant_index_type LINEAR -variant_index_parameter 128000 -o /dev/stdout | bgzip -s - > data/SampleID.calls.gvcf.gz
```
## Joint calling of all samples
```
bcftools index -t SampleID.calls.gvcf.gz

gatk -T GenotypeGVCFs -R ninespine_v6b.fa -V Sample1.calls.gvcf.gz ..... -V Sample889.calls.gvcf.gz -o /dev/stdout | bcftools view -Oz -o Samples.calls.vcf.gz &
```

## Filter the raw output vcf file
### Remove duplicates and keep only biallelic SNPs with quality control
```
bcftools view -T ^outgroup_sites_to_remove.tab Samples.calls.vcf.gz -Oz -o Samples.ogrm.vcf.gz 

bcftools filter -g 10 Samples.ogrm.vcf.gz | bcftools view -m2 -M2 -v snps -Oz -o Samples.ogrm.clean1.vcf.gz
 
bcftools filter -g 20 Samples.ogrm.clean1.vcf.gz -Oz -o Samples.ogrm.clean2.vcf.gz

bcftools view -m2 -M2 -V mnps,indels -T ^repeats.tab Samples.ogrm.clean2.vcf.gz -Oz -o nsp889_clean.vcf.gz
```
### Exclude sex-chormosome data
```
bcftools view -t ^LG12 nsp889_clean.vcf.gz -Oz -o nsp889_clean_auto.vcf.gz
```
### Separate each population and do some filterings
```
bcftools view -S ~/pops/${pop} ~/nsp889_clean_auto.vcf.gz -Oz -o ${pop}.vcf.gz 
taibx -p vcf ${pop}.vcf.gz
vcftools --gzvcf ${pop}.vcf.gz --minGQ 20 --minQ 30 --minDP 8 --recode --recode-INFO-all -c | bcftools view -Oz -o ${pop}_minQ30_minDP8_minGQ20.vcf.gz
```
### Merge only the 888 ingroup samples
```
bcftools merge ~/*_minQ30_minDP8_minGQ20.vcf.gz -Oz -o nsp888_minQ30_minDP8_minGQ20.vcf.gz
```
### Find ingroup sites in the outgroup and merge individual files
```
bcftools view -R nsp888_minQ30_minDP8_minGQ20.vcf.gz PUN-TYM_minQ30_minDP8_minGQ20.vcf.gz  -Oz -o PUN-TYM_minQ30_minDP8_minGQ20_ingroup.vcf.gz
bcftools merge nsp888_minQ30_minDP8_minGQ20.vcf.gz PUN-TYM_minQ30_minDP8_minGQ20_ingroup.vcf.gz -Oz -o nsp889_minQ30_minDP8_minGQ20.vcf.gz
```
