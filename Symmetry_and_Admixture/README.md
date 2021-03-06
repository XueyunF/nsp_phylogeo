# Analysis with ADMIXTOOLS

## Data filtering for ADMIXTOOLS related analysis
```
vcftools --gzvcf nsp889_minQ30_minDP8_minGQ20.vcf.gz --max-meanDP 25 --max-missing 0.75 --thin 300 --stdout --recode --recode-INFO-all | bcftools view --min-ac=1 -Oz -o SNPset4.vcf.gz
```
### Remove the LG elements in the vcf file to avoid errors and then convert the format
```
zcat SNPset4.vcf.gz | sed 's/^LG//' | bgzip > SNPset4_noLG.vcf.gz
plink --vcf SNPset4_noLG.vcf.gz --recode --out SNPset4 --noweb
plink --vcf SNPset4_noLG.vcf.gz --make-bed --out SNPset4 --noweb 
```
### Then convert the files to eigenstrat format
```
convertf -p convert
```

#### Content of parameters in parameter file 'convert'
```
genotypename: SNPset4.bed
snpname:      SNPset4.bim
indivname:    SNPset4.ped
genotypeoutname: SNPset4.geno
snpoutname:      SNPset4.snp
indivoutname:    SNPset4.ind
familynames:     NO
```
## Analyses with [ADMIXTOOLS](https://github.com/DReichLab/AdmixTools), you can also use [AdmixR](https://github.com/martinsikora/admixr) for the same analyses
#### Then run analysis with default settings (default parameter files can be found [here](https://github.com/DReichLab/AdmixTools/tree/master/examples)
```
qp3Pop -p parF3 > parF3_outgroup.out
qpDstat -p parDstat > parDstat.out
qpF4ratio -p parF4ratio > parF4ratio.out
qpWave -p parWave > parWave.out
qpAdm -p parAdm > parAdm.out
```
### Modelling with qpGraph and plot the models
```
qpGraph -p parGraph -g input_Graph_model -o output_Graph_model.ggg -d output_Graph_model.dot > output_Graph_model.log
dot -Tpng output_Graph_model.dot -o output_Graph_model.png
```

#### Settings in parameter files parGraph
```
indivname:    SNPset4.ind  
snpname:      SNPset4.snp
genotypename: SNPset4.geno
outpop:  NULL
blgsize: 0.05
forcezmode: YES
lsqmode: YES
diag:  .0001
bigiter: 6
hires: YES
lambdascale: 1
initmix:      1000
precision:    .0001
zthresh:      3.0
terse:        NO
useallsnps: NO
```

## Analysis with [Dsuite](https://github.com/millanek/Dsuite)

### Data filtering for Dsuite related analysis
```
vcftools --gzvcf nsp889_minQ30_minDP8_minGQ20.vcf.gz --max-meanDP 25 --max-missing 0.9 --recode --recode-INFO-all -c | bcftools view --min-ac=1 -Oz -o SNPset5.vcf.gz
```
#### First calculate D and F4ratio tests 
```
Dsuite Dtrios SNPset5.vcf.gz sets.txt -o fullsetmm90_ASTRAL -t ASTRAL_root_Puntym.tree  #The file 'sets.txt' is a file which assign individuals to populations, and the outgroup for the analysis
```
#### Then calculate f-branch with the ASTRAL tree as the input phylogeny
```
Dsuite Fbranch ASTRAL_root_Puntym.tree fullsetmm90_ASTRAL__tree.txt > fbranch_ASTRAL.txt
```
#### Plot f-branch results with or without ordering, scripts can be found [here](https://github.com/millanek/Dsuite/tree/master/utils)
```
python3 dtools.py fbranch_ASTRAL.txt ASTRAL_root_Puntym.tree -n ASTRAL_root_Puntym_ladderized --ladderize
python3 dtools.py fbranch_ASTRAL.txt ASTRAL_root_Puntym.tree -n ASTRAL_root_Puntym_nonladderized
```

## Rare allele sharing statistics [(RASS)](https://github.com/TCLamnidis/RAStools)
### Convert data to required format with [Rarecoal-tools](https://github.com/stschiff/rarecoal-tools)
```
bcftools view -m2 -M2 -c1 -v snps SNPset4.vcf.gz | vcf2FreqSum | groupFreqSum -f groups.txt -m 0.25 > grouped_data.freqsum 
sed 's/^LG//' grouped_data.freqsum > grouped_data_noLG.freqsum
```
### Then run RASS analysis for each of the test populations 
```
python3 ~/RASCalculator.py -I grouped_data_noLG.freqsum \
-o PUN-TYM \
-a CAN-TEM,CAN-FLO,RUS-LEN,JAP-BIW,USA-HLA,RUS-LEV,GBR-GRO \
-L DEN-NOR,LAT-JAU,RUS-BOL,RUS-KRU,RUS-MAS,FIN-RYT,FIN-PYO,FIN-KRK,FIN-KAR,FIN-UKO,FIN-PUL,FIN-KEV,FIN-RII,SWE-KIR,SWE-NAV,SWE-HAN,SWE-ABB,SWE-BYN,FIN-KIV,FIN-HAM,SWE-BOL,EST-PUR,FIN-TVA,FIN-HEL,FIN-SEI,SWE-GOT,LAT-JAU,POL-GDY,GER-RUE,SWE-LUN,DEN-RES,SWE-FIS,NOR-KVN,NOR-UGE,DEN-NOR,FRA-VEY,BEL-MAL,SCO-HAR,NOR-ROY,NOR-ENG \
-C 21 \
--skipJackknife \
-c DEN-NOR,LAT-JAU,RUS-BOL,RUS-KRU,RUS-MAS,FIN-RYT,FIN-PYO,FIN-KRK,FIN-KAR,FIN-UKO,FIN-PUL,FIN-KEV,FIN-RII,SWE-KIR,SWE-NAV,SWE-HAN,SWE-ABB,SWE-BYN,FIN-KIV,FIN-HAM,SWE-BOL,EST-PUR,FIN-TVA,FIN-HEL,FIN-SEI,SWE-GOT,LAT-JAU,POL-GDY,GER-RUE,SWE-LUN,DEN-RES,SWE-FIS,NOR-KVN,NOR-UGE,DEN-NOR,FRA-VEY,BEL-MAL,SCO-HAR,NOR-ROY,NOR-ENG \
> RASS_TEST_allpopssingle.out
```
