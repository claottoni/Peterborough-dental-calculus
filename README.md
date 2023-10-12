# Metagenomic analyis of dental calculus samples from Peterborough, England
Main pipelines and command lines used for the analysis of Late Medieval ancient dental calculus samples from the leprosarium of St Leonard, Peterborough (England)
The bioinformatic analyses were done on the Galileo100 supercomputing cluster of Cineca, with the support of Elixir-Italy and the HPC@CINECA program.

## Table of contents
https://github.com/claottoni/Peterborough-dental-calculus/blob/main/README.md#pre-processing-of-raw-sequencing-data

## Pre-processing of raw sequencing data
Sequencing data from four individuals and two negative controls consisted of paired-end sequencing reads, which were quality-filtered, trimmed of the adapter sequences and merged with AdapterRemoval: 

```bash
AdapterRemoval --file1 filename_1.fastq.gz --file2 filename_2.fastq.gz --basename filename --minlength 30 --minquality 25 --trimns --trimqualities --gzip --threads 16 --collapse
```

## Taxonomic classification of the reads

### Construction of Kraken2 custom database
We constructed a custom database of complete genomes of Bacteria, Archaea, Fungi and Viruses from the NCBI RefSeq and GenBank with Kraken2. However, kraken-build downloads only genomes from RefSeq reported as "Complete Genome". We also included complete genomes from GenBank as follows: 

```bash
# download assembly summaries
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -O assembly_summary_gb_bacteria_Jan2022.txt
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt -O assembly_summary_rs_bacteria_Jan2022.txt
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt -O assembly_summary_gb_archaea_Jan2022.txt
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt -O assembly_summary_rs_archaea_Jan2022.txt

# print (and save in a file) the accession numbers of complete genomes in RefSeq and GenBank, removing the first four characters (GCA_ and GCF_ prefixes)
grep "Complete Genome" assembly_summary_gb_archaea_Jan2022.txt | awk -F "\t" '{print $1}' | cut -c 5- > accession_gb_archaea.txt
grep "Complete Genome" assembly_summary_rs_archaea_Jan2022.txt | awk -F "\t" '{print $1}' | cut -c 5- > accession_rs_archaea.txt 
grep "Complete Genome" assembly_summary_gb_bacteria_Jan2022.txt | awk -F "\t" '{print $1}' | cut -c 5- > accession_gb_bacteria.txt
grep "Complete Genome" assembly_summary_rs_bacteria_Jan2022.txt | awk -F "\t" '{print $1}' | cut -c 5- > accession_rs_bacteria.txt

# print (and save in a file) the accession numbers of complete genomes in GenBank that are not present in RefSeq (the genomes to add to the database)
grep -vFf accession_rs_archaea.txt accession_gb_archaea.txt > complete_to_add_from_gb_archaea.txt
grep -vFf accession_rs_bacteria.txt accession_gb_bacteria.txt > complete_to_add_from_gb_bacteria.txt

# print (and save in a file) the ftp directory paths where the complete genomes to add from GenBank are stored 
grep -Ff complete_to_add_from_gb_archaea.txt assembly_summary_gb_archaea_Jan2022.txt | awk -F "\t" '{print $20}' > ftpdirpaths_complete_to_add_archaea_gb_Jan2022.txt
sed 's/^......//' ftpdirpaths_complete_to_add_archaea_gb_Jan2022.txt > ftpdirpaths_complete_to_add_archaea_gb_Jan2022_final.txt
grep -Ff complete_to_add_from_gb_bacteria.txt assembly_summary_gb_bacteria_Jan2022.txt | awk -F "\t" '{print $20}' > ftpdirpaths_complete_to_add_bacteria_gb_Jan2022.txt
sed 's/^......//' ftpdirpaths_complete_to_add_bacteria_gb_Jan2022.txt > ftpdirpaths_complete_to_add_bacteria_gb_Jan2022_final.txt

# download the genomes in a folder named archaea 
while read i; do
  rsync --copy-links --recursive --times --verbose rsync:${i} archaea
done < ftpdirpaths_complete_to_add_archaea_gb_Jan2022_final.txt
# download the genomes in a folder named bacteria 
while read i; do
  rsync --copy-links --recursive --times --verbose rsync:${i} bacteria
done < ftpdirpaths_complete_to_add_bacteria_gb_Jan2022_final.txt
```

The genomes from RefSeq and GenBank (added sequences) were included in the database following the manual of Kraken2. Teh database was built setting the maximum database size to 300 Gb. The commands were run in the Galileo100 supercomputing cluster with a batch script submitted to Slurm setting 300 Gb of memory. 

```bash
DBNAME=customkraken2_300Gb_Jan2022
DBSIZE=300000000000	
THREADS=16
# Donwload taxonomy
kraken2-build --download-taxonomy --threads $THREADS --db $DBNAME
# Download libraries (complete genomes form the NCBI RefSeq).
kraken2-build --download-library bacteria --threads $THREADS --db $DBNAME --no-masking
kraken2-build --download-library viral --threads $THREADS --db $DBNAME
kraken2-build --download-library fungi --threads $THREADS --db $DBNAME
kraken2-build --download-library archaea --threads $THREADS --db $DBNAME
# Add the complete genomes from GenBank
for i in genomes_to_add/bacteria/fasta/*.fna
do 
  kraken2-build --add-to-library $i --threads $THREADS --db $DBNAME
done
for i in genomes_to_add/archaea/fasta/*.fna
do 
  kraken2-build --add-to-library $i --threads $THREADS --db $DBNAME
done
# All the genomes were masked for low complexity regions following the commands present in the mask_low_complexity.sh script of Kraken2
for i in `find $DBNAME/library \( -name '*.fna' -o -name '*.ffn' \)`
do
  dustmasker -in $i -infmt fasta -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > tempfile
  mv -f tempfile $i
done
# Finally, the database was built
kraken2-build --build --threads $THREADS --db $DBNAME --max-db-size $DBSIZE 
```
The final database size was 122 Gb. 

### Classification of the reads with Kraken2 and quantification with Bracken
Reads previously processed and merged with AdapterRemoval were submitted to taxonomic classification with Kraken2 against the custom database previously built. 

```bash
KRAKEN_DB=customkraken2_300Gb_Jan2022
OUTPUT=/path/to/output-folder
for i in $(find -name "*collapsed.gz" -type f)
do
  FILENAME=$(basename "$i")
  SAMPLE=${FILENAME%.gz}
  kraken2 --db $KRAKEN_DB $i --gzip-compressed --output ${OUTPUT}/${SAMPLE}.krk --report ${OUTPUT}/${SAMPLE}.krk.report
done
```
Then, we ran Bracken to obtain actual quantification data, as described in the manual of the program. We first built the Bracken database, and then ran the samples. 
We used a read length of 65 bp and a threshold of 50 reads. Modern metagenome data were run with read lentgh of 90 bp. 
```bash
KRAKEN_DB=customkraken2_50Gb_2020
THREADS=16
READ_LEN=65
THRESHOLD=50
# Build Bracken database
bracken-build -d ${KRAKEN_DB} -t ${THREADS} -l ${READ_LEN}
# Run Bracken from the kraken reports.
for i in $(find -name "*krk.report" -type f)
do
  FILENAME=$(basename "$i")
  SAMPLE=${FILENAME%.krk.report}
  bracken -d $KRAKEN_DB -i $i -o ${SAMPLE}.${READ_LEN}mer.bracken -r $READ_LEN -t $THRESHOLD
done
```

### Parse the Bracken outputs to abundance tables
To create an abundance table comprehensive of all the taxa abundances for every sample we used a custom R script, `brackenToAbundanceTable_v2.R` (available in toolbox repository), which merges the species abundances contained in the Bracken output of each sample in one table. The script needs as argument the path to the folder containing the Bracken results. When working in the folder containing the Bracken outputs: 
```bash
brackenToAbundanceTable_v2.R .
```
The script generates the abundance table file named taxa_abundance_bracken_names_IDs.txt, which contains the species names and the NCBI IDs in each row, and the samples in the columns.
All the fastq files from the literature used as comparative dataset as reported in the manuscript were submitted to the same workflow of taxonomic classification with Kraken2. All the abundance tables generated from each dataset were merged with the custom script `abundanceTablesMerger_v2.R`, as follows: 
```bash
abundanceTablesMerger_v2.R *.txt
```
The output is a file named abundance_table_IDs.merged, which is used for downstream analyses.

### Normalization of the abundance table for genome length and full taxonomy description
We normalized the abundance of each species for its genome lenth. To do that we used a Python script (`gL-normalizer-lite_v3.py`, available in the toolbox repository) that takes three arguments: 1) the abundance table with names of the species, 2) the table of genome lengths (parsed from the NCBI genome browser, available in the toolbox repository), 3) the name of the output file. The script searches the names of the species in the file with the list of genome lengths of species in the NCBI and divides the abundance of each species by the length of the assembly available in the NCBI. 
```bash
TABLE=gL-prokaryotes-viruses-Jun2022.table
python gL-normalizer-lite_v3.py abundance_table_IDs.merged $TABLE abundance_table_IDs.merged.norm
```
The output table is named `abundance_table_IDs.merged.norm`. The table reports the Species, the NCBI IDs, the Assembly stage of the genome in NCBI, the Length in Megabases, the kind of Match, exact, species (when a peculiar strain is not found in the list of lenghts), and genus (when a species is not found). Some species may not be contained in the table of the genome lengths (for example if the taxonomy names have been updated or have ambiguities in the species name). The species not normalized are removed from the table, you can check them and their frequency (they are normally very rare) in the `gLnotNormalized.log` file.

Then, we retrieved the full taxonomic ranks of the species in the abundance table. This is useful in the downstream steps to make analyses at taxonomic ranks higher than species (e.g. genus, phylum) in Phyloseq. To retrieve the full taxonomic ranks you must have the program `Taxaranks` installed (https://github.com/linzhi2013/taxonomy_ranks/blob/master/README.md). We incoroprated taxaranks in a custom script `getFullTaxaranks_v2.sh` that attaches the full taxonomy to every species ID in the normalized abudance table. Taxaranks takes the NCBI ID corresponding to each species in the table and searches the full taxonomy associated to that ID.
```bash
getFullTaxaranks_v2.sh -i abundance_table_IDs.merged.norm
```
The script generates three files: 1) `abundance_table_IDs.merged.norm.taxonomy` - the full taxonomic ranks (up to Kingdom) of the species entries; 2) abundance_table_IDs.merged.norm.taxonomy.err - the species for which no taxonomic ranks could be found; 3) abundance_table_IDs.merged.norm.taxonomy.final - the original abundance table with full taxonomy ranks attached. 
Finally, we focused the downstream analysis on Archaea and Bacteria by removing Viruses and Fungi taxa from the table: 
```bash
grep -v "Viruses\|Fungi" abundance_table_IDs.merged.norm.taxonomy.finall > abundance_table_IDs.merged.norm.taxonomy.final.Archaea_Bacteria
```

## Analysis of metagenomic data in R
We used the R package Phyloseq to analysed the normalized abundance table of Archaea and Bacteria. 

### Preparation of Phyloseq datasets
Phyloseq works with `biom` objects that store all the information required for the analyses: abundance of taxa, sample metadata, taxonomy. We imported the table with the species abundances and the full taxonomic ranks, which has the taxa in the rows and the samples in the columns. Note that we set the column #20, corresponding to the species names, as row names.

```R
library(phyloseq)
library(ggplot2)
library(vegan)
# import the taxonomy file generated by taxaranks
table = as.matrix(read.delim("abundance_table_IDs.merged.norm.taxonomy.final.Archaea_Bacteria", header=T, fill=T, row.names=20, sep="\t"))
# extract the taxonomy by selecting the corresponding columns.
taxonomy = (table[,c(3:19)])
# extract the otu by removing the columns corresponding to the taxonomy and other unwanted columns.
otu = (table[,-c(1:23)])
# make the matrix values as numeric
class(otu) <- "numeric"
# create OTU and TAX objects in Phyloseq
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxonomy)
my_biom = phyloseq(OTU, TAX)
# import mapping file with metadata
biom_metadata <- import_qiime_sample_data("metadata.txt")
# merge data
my_biom <- merge_phyloseq(my_biom, biom_metadata)
# filter low-abundance taxa below 0.02%
minTotRelAbun = 0.0002
x = taxa_sums(my_biom)
keepTaxa = taxa_names(my_biom)[which((x / sum(x)) > minTotRelAbun)]
my_biom_flt = prune_taxa(keepTaxa, my_biom)
# convert to relative freq
my_biom_rel = transform_sample_counts(my_biom_flt, function(x) x / sum(x))
```

### Non-metric Multidimensional Scaling
We set the following colour settings associated with the groups listed in the metadata:
```R
bg=c("royalblue",		
"blue",					
"mediumpurple 3",				
"#9C661F",				
"darkgreen",			
"blueviolet",			
"#9C661F",				
"#808080",				
"red",				
"palevioletred",				
"darkred")				

coul=c("skyblue 1",		
"blue",					
"mediumpurple",				
"#9C661F",				
"darkgreen",			
"blueviolet",			
"#9C661F",				
"#D3D3D3",				
"red",				
"palevioletred",				
"darkred")				
```

We made a  clr-transformation of the abundance data with the library `microbiome`, calculated the Aitchinson distances and the nMDS as follows: 
```R
library(microbiome)
my_biom_rel_clr = transform(my_biom_rel, "clr")
distAIT = distance(my_biom_rel_clr, method = "euclidean")
ordAIT = ordinate(my_biom_rel_clr, method = "NMDS", distance = distAIT)
```

Finally we plot the nMDS. We used the groups as reported in the `metadata.txt` file (column Group1).
```R
p=plot_ordination(my_biom_rel_clr, ordAIT, color = "Group1") 
# further refining of the chart
p=p + theme_bw() +
  geom_point(aes(shape=Group1, color=Group1, fill=Group1), size=2.2) +
  scale_shape_manual(values=c(21,21,21,10,3,4,0,21,8,2,13)) +
  scale_color_manual(values=bg) +
  scale_fill_manual(values=coul) +
  labs(color  = "Group1", shape = "Group1")   #merge legend
# remove overlapping layer
p$layers <- p$layers[-1]
```

### Differential taxonomic abundances with DESeq2
We used DESeq2 to identify differentially abundant species in the datasets examined. Pyloseq was used to select only the datasets of interest from the original `my_biom` object. 

```R
# subsample oral microbiomes 
subgroup = c("Peterborough","Dental_calculus_modern","Dental_calculus_human_18-19th_c.", "Dental_calculus_medieval","Oral_plaque")
subsample = subset_samples(my_biom, Group1 %in% subgroup)
# Filter low-abundance taxa (0.02% threshold).
minTotRelAbun = 0.0002
x = taxa_sums(subsample)
keepTaxa = taxa_names(subsample)[which((x / sum(x)) > minTotRelAbun)]
subsample_flt = prune_taxa(keepTaxa, subsample)
# make a +1 offset to run Deseq2 (which does not tolerate zero counts)
otu_table(subsample_flt) = otu_table(subsample_flt)+1  
# make deseq object from phyloseq object
ds = phyloseq_to_deseq2(subsample_flt, ~ Group1)
# Run DESeq2
dds.data = DESeq(ds)
```

We first compared Peterborough samples against the modern dental calculus samples: 
```R
# 1) Peterborough VS Dental_calculus_modern
# With the 'contrast' function you screen two different set of samples (based on your metadata) for differential taxa. 
res = results(dds.data, contrast=c("Group1","Peterborough","Dental_calculus_modern"))
# sort based on p-value adjusted:
resOrdered = res[order(res$padj, na.last=NA), ]
# set a threshold value for the False Discovery Rate (FDR):
alpha = 0.01
# get only significant taxa based on p-value adjusted (the FDR):
resSig <- subset(resOrdered, padj < alpha)
# sort significant values based on the log2-fold-change:
resfc = resSig[order(resSig$log2FoldChange),]
# sort significant values based on abundance (the base mean):
resbm = resSig[order(resSig$baseMean),]
# save the the tables 
write.table(as.data.frame(resbm), file="deseq2_Peterborough_modern.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=F)
  # convert the data in a dataframe
resSigMod = cbind(as(resSig, "data.frame"), as(tax_table(subsample_flt)[rownames(resSig), ], "matrix"))
# Plot log-fold-changes of the OTUs based on species
p1=ggplot(resSigMod, aes(x=species, y=log2FoldChange, color=phylum)) +
    geom_jitter(size=3, width = 0.2) +
    theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5))
```

Then we compared Peterborough samples against the 18th-19th century dental calculus samples from Radcliffe: 
```R
# 2) Peterborough VS Velsko
# With the 'contrast' function you screen two different set of samples (based on your metadata) for differential taxa. 
res = results(dds.data, contrast=c("Group1","Peterborough","Dental_calculus_human_18-19th_c."))
# sort based on p-value adjusted:
resOrdered = res[order(res$padj, na.last=NA), ]
# set a threshold value for the False Discovery Rate (FDR):
alpha = 0.01
# get only significant taxa based on p-value adjusted (the FDR):
resSig <- subset(resOrdered, padj < alpha)
# sort significant values based on the log2-fold-change:
resfc = resSig[order(resSig$log2FoldChange),]
# sort significant values based on abundance (the base mean):
resbm = resSig[order(resSig$baseMean),]
# save the the tables 
write.table(as.data.frame(resbm), file="deseq2_Peterborough_Velsko.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=F)
  # convert the data in a dataframe
resSigMod = cbind(as(resSig, "data.frame"), as(tax_table(subsample_flt)[rownames(resSig), ], "matrix"))
# Plot log-fold-changes of the OTUs based on species
p2=ggplot(resSigMod, aes(x=species, y=log2FoldChange, color=phylum)) +
    geom_jitter(size=3, width = 0.2) +
    theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5))
```

Finally we compared Peterborough with the Medeival dental calculus dataset form Ireland
```R
# 3) Peterborough VS Dental_calculus_medieval
# With the 'contrast' function you screen two different set of samples (based on your metadata) for differential taxa. 
res = results(dds.data, contrast=c("Group1","Peterborough","Dental_calculus_medieval"))
# sort based on p-value adjusted:
resOrdered = res[order(res$padj, na.last=NA), ]
# set a threshold value for the False Discovery Rate (FDR):
alpha = 0.01
# get only significant taxa based on p-value adjusted (the FDR):
resSig <- subset(resOrdered, padj < alpha)
# sort significant values based on the log2-fold-change:
resfc = resSig[order(resSig$log2FoldChange),]
# sort significant values based on abundance (the base mean):
resbm = resSig[order(resSig$baseMean),]
# save the the tables 
write.table(as.data.frame(resbm), file="deseq2_Peterborough_Velsko.tsv", sep="\t", row.names=TRUE, col.names=TRUE, quote=F)
  # convert the data in a dataframe
resSigMod = cbind(as(resSig, "data.frame"), as(tax_table(subsample_flt)[rownames(resSig), ], "matrix"))
# Plot log-fold-changes of the OTUs based on species
p3=ggplot(resSigMod, aes(x=species, y=log2FoldChange, color=phylum)) +
    geom_jitter(size=3, width = 0.2) +
    theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5))
```

We combined the plot with `ggarrange`
```R
ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 3, nrow = 1, widths = c(0.8,1,0.4))
```

### PERMANOVA and dispersion test
The PERMANOVA was done in R with `adonis` test. The dispersion test for variance homogeneity of the samples comared was done with `betadisper`. As already done in DEseq2, we contrasted the sample from Peterborough against those from modern dental calculus, 18th-19th century UK and Medieval Ireland. 
```R
#-1 Peterborough VS modern
subgroup1 = c("Peterborough","Dental_calculus_modern")
subsample1 = subset_samples(my_biom_rel_clr, Group1 %in% subgroup1)
# estimante Aitchinson dist (clr+euclidean)
distAIT1 = phyloseq::distance(subsample1, method = "euclidean")
# make permanova 
subsample1.df = data.frame(sample_data(subsample1))
adonis.results1 = adonis(distAIT1 ~ Group1, data = subsample1.df)
write.table(as.data.frame(adonis.results1$aov.tab),"adonis_Peterborough_modern.tsv", sep="\t",row.names=TRUE,col.names=TRUE)
#--- betadisper to test dispersion homogeneity
mod1 <- betadisper(distAIT1, subsample1.df$Group1)
dispersion1 = anova(mod1)  #Anova of distances to centroids to interpret the significance of F --> #NOT SIGNIFICANT: GROUP DISPERSIONs ARE HOMOGENEOUS
write.table(as.data.frame(dispersion1),"dispersion_Peterborough_modern.tsv", sep="\t",row.names=TRUE,col.names=TRUE)

#-2 Peterborough VS 18th-19th century UK
subgroup2 = c("Peterborough","Dental_calculus_human_18-19th_c.")
subsample2 = subset_samples(my_biom_rel_clr, Group1 %in% subgroup2)
# estimante Aitchinson dist (clr+euclidean)
distAIT2 = phyloseq::distance(subsample2, method = "euclidean")
# make permanova 
subsample2.df = data.frame(sample_data(subsample2))
adonis.results2 = adonis(distAIT2 ~ Group1, data = subsample2.df)
write.table(as.data.frame(adonis.results2$aov.tab),"adonis_Peterborough_Velsko.tsv", sep="\t",row.names=TRUE,col.names=TRUE)
#--- betadisper to test dispersion homogeneity
mod2 <- betadisper(distAIT2, subsample2.df$Group1)
dispersion2 = anova(mod2)  #Anova of distances to centroids to interpret the significance of F --> #NOT SIGNIFICANT: GROUP DISPERSIONs ARE HOMOGENEOUS
write.table(as.data.frame(dispersion2),"dispersion_Peterborough_Velsko.tsv", sep="\t",row.names=TRUE,col.names=TRUE)

#-3 Peterborough VS Medieval Ireland
subgroup3 = c("Peterborough","Dental_calculus_medieval")
subsample3 = subset_samples(my_biom_rel_clr, Group1 %in% subgroup3)
# estimante Aitchinson dist (clr+euclidean)
distAIT3 = phyloseq::distance(subsample3, method = "euclidean")
# make permanova 
subsample3.df = data.frame(sample_data(subsample3))
adonis.results3 = adonis(distAIT3 ~ Group1, data = subsample3.df)
write.table(as.data.frame(adonis.results3$aov.tab),"adonis_Peterborough_Mann.tsv", sep="\t",row.names=TRUE,col.names=TRUE)
#--- betadisper to test dispersion homogeneity
mod3 <- betadisper(distAIT3, subsample3.df$Group1)
dispersion3 = anova(mod3)  #Anova of distances to centroids to interpret the significance of F --> #NOT SIGNIFICANT: GROUP DISPERSIONs ARE HOMOGENEOUS
write.table(as.data.frame(dispersion3),"dispersion_Peterborough_Mann.tsv", sep="\t",row.names=TRUE,col.names=TRUE)
```

## Antibiotic microbial resistance analysis
We downloaded the `nucleotide_fasta_protein_homolog_model.fasta` from the CARD database and we replaced the spaces in the db with underscore, so as to keep the full info about the hits in the blastn analysis: 
```bash
sed 's/ /_/g' nucleotide_fasta_protein_homolog_model.fasta > nucleotide_fasta_protein_homolog_model.fasta_mod.fasta
```

We built the database with blast+
```bash
makeblastdb -in nucleotide_fasta_protein_homolog_model.fasta_mod.fasta -dbtype nucl -out card_nucl_Jan2023
```

We also created a database of three housekeeping genes to be used as control in the ancient samples (and normalized the results to account for possible taphonomic effects in chornologicaly different samples). We downloaded the sequences with `efetch` and `esearch` and built the database with blast+

```bash
esearch -db nucleotide -query '"recA"[Protein name]' | efetch -format fasta > recA_nu.fasta
esearch -db nucleotide -query '"gyrB"[Protein name]' | efetch -format fasta > gyrB_nu.fasta
esearch -db nucleotide -query '"rpoB"[Protein name]' | efetch -format fasta > rpoB_nu.fasta
cat *.fasta > recA_gyrB_rpoB_nu.fasta
# build the database
makeblastdb -in recA_gyrB_rpoB_nu.fasta -dbtype nucl -out recA_gyrB_rpoB_nu
```

We extracted the reads classified by Kraken2 as bacteria with `extract_kraken_reads.py`, from KrakenTools (https://github.com/jenniferlu717/KrakenTools). The bacterial reads were converted in `fasta` format.

```bash
KRK=filename.krk
FASTQ=filename.fastq.gz
OUTPUT=path/to/folder
FILENAME=$(basename "$KRK")
FNAME="${FILENAME%.krk}"
extract_kraken_reads.py -k $KRK -s $FASTQ -o ${OUTPUT}/${FNAME}.bacteria.fasta -t 2 --include-children -r ${KRK}.report
```

We ran the analysis in blast twice for each sample, once with the card database, and once with the housekeeping genes database. 

```bash
FASTA=fasta_file
OUTPUT=output_sample_name
# run against card db
DBNAME=/g100_scratch/userexternal/cottoni0/databases/card_aro_Jan2023/card_nucl_Jan2023
blastn -query $FASTA -db $DBNAME -out ${OUTPUT}.aro.blastn.out -outfmt 6
# run against housekeeping genes db
DBNAME=/g100_scratch/userexternal/cottoni0/databases/card_aro_housekeeping/recA_gyrB_rpoB_nu
blastn -query $FASTA -db $DBNAME -out ${OUTPUT}.hk.blastn.out -outfmt 6
```

A custom python script (available in toolbox repository) was used to normalize for the sequencing depth the hits to the housekeeping database for each sample of the comparative database used and generate the file `Peterborough_full_dataset_hk_genes.txt`. 
```bash
amr_blastn_parser_v2.py *.blastn.out > Peterborough_full_dataset_hk_genes.txt
```

We retrieved the aro categories for the hits in the card database and generated a `file.index` with a custom script `aro_blastn_parser.py`.
```bash
ARO_INDEX=card_download_Jan2023/aro_index.tsv
for i in $(find -name "*aro.blastn.out" -type f)
do 
  filename=$(basename "$i")
  fname="${filename%.out}"
  echo "parsing $i"
  aro_blastn_parser.py $i $ARO_INDEX > $(dirname "$i")/${fname}.out.index
done
```

Reads-hits in the blastn output are sorted by decreasing bitscore (last column). Sort and uniq the read-names (the first with the highest bitscore will be returned). 
```bash
for i in $(find -name "*.out.index" -type f); do echo $i; sort -u -k1,1 $i | awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$4}' > ${i}.uniq; done
```

The following analyses were done in R. We fits prepared a table of the samples by merging all the files from each dataset (from this study and the literature) to generate the file `Peterborough_full_dataset_amr_genes.txt` 

```R
files <- list.files(pattern="*.uniq$", full.names=T,recursive=FALSE)
for (i in files) {
  if (!exists("tabfinal")){
    tab=read.delim(i, header=F, fill=T, row.names=NULL, sep="\t")
    tabfinal = table(tab$V2)
    print(paste0("parsing ", files[c(1,2)]))
    } else {
  # merge all the others
  tab=read.delim(i, header=F, fill=T, row.names=NULL, sep="\t")
  tabfinal = merge(tabfinal, table(tab$V2), by=1, all=T)}
  #print(paste0("parsing ", i)
}
colnames(tabfinal) = c("gene",files)
write.table(tabfinal, file = 'Peterborough_full_dataset_amr_genes.txt', sep="\t", row.names=F, na="0", col.names=T, quote = FALSE)
```

Then we imported the two tables with the number of hits in the card database (AMR genes) of and in the housekeeping genes database and  in R and made some changes in the tables in order to make them fully consistent in terms of order of the samples. We did that in order to normalize the number of card hits in the following step. 
```R
amr_genes = read.delim("Peterborough_full_dataset_amr_genes.txt", header=T, fill=T, row.names=NULL, sep="\t")
hk3_genes = read.delim("Peterborough_full_dataset_hk_genes.txt", header=T, fill=T, row.names=NULL, sep="\t")
amr_genes[is.na(amr_genes)] <- 0
amr_genes.final = amr_genes		
row.names(amr_genes.final) = amr_genes.final[,1]
amr_genes.final = amr_genes.final[,-1]
# transpose the amr genes table
amr_genes.final = t(amr_genes.final)
#reorder hk3_genes samples to be consistent with amr 
hk3_genes.sorted = hk3_genes[c(2,1,4:6,8,9:53,57:87),]
write.table(hk3_genes.sorted, "Peterborough_full_dataset_hk_genes.sorted.txt", quote=F, sep="\t", row.names=F, col.names=T)
```

Finally we normalized the hits in the AMR genes with the housekeeping genes hits and created a dataframe
```R
amr_genes.final.norm = amr_genes.final/hk3_genes.sorted$Hits*1000000
amr_genes.final.norm = as.data.frame(amr_genes.final.norm)
write.table(amr_genes.final.norm, "Peterborough_AMR_hk3norm_cpm.tsv", quote=F, sep="\t", row.names=T, col.names=T)
```

To analyse the data in R we first created a vector with the group metadata of each sample and we attached it to the dataframe (coumn named "group")
```R
amr_genes.final.norm$group = c("Peterborough","NTC","Peterborough","Peterborough","Peterborough","NTC","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","UK-18-19th_c.","Modern","Modern","Modern","Modern","Modern","Modern","Modern","Modern","Modern","Modern","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Ireland-Medieval","Modern","Modern","Modern","Modern","Modern","Modern","Modern","Modern")
```

We made a Wilcoxon Rank Sum Tests with Benjamini & Hochberg adjusted P-values by selecting the columns in the dataframe corresponding to the genes. 
```R
# make wilcoxon test
wilcox = lapply(amr_genes.final.norm[,c(1:33)], function(x) pairwise.wilcox.test(x, amr_genes.final.norm$group, p.adjust.method = "BH"))
library(plyr)
out <- ldply(wilcox, function(x) x$p.value)
# save table
write.table(out, "Peterborough_wilcox_test_hk3_norm_full.tsv", quote=F, sep="\t", row.names=F, col.names=T)
```

To generate a multipanel plot of the results:
```R
# Plot all boxplots in multipanel
library(tidyverse)
#run pivot longer to prepare the data by including all the gene families to display in the plot.
df.long <- amr_genes.final.norm %>% 
  pivot_longer("ABC-F ATP-binding cassette ribosomal protection protein":"Rm3 family beta-lactamase", names_to = 'variable', values_to = 'value')
# we changed the order for displaying the charts
df.long$group <- factor(df.long$group, levels = c("NTC", "Peterborough", "Ireland-Medieval", "UK-18-19th_c.", "Modern"))
# generate the multipanel plots
g = ggplot(data = df.long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(lwd=0.1,
  outlier.size = 0.5,
  outlier.stroke = 0.1) +
  theme(axis.text.x = element_blank(),
  strip.text = element_text(size = 5),
  legend.title = element_text(colour = "black", size = 7, face = "bold"),
  legend.text = element_text(colour = "black", size = 6),
  axis.line = element_line(size=0.1),
  axis.ticks = element_line(size=0.1),
  axis.title.x = element_blank(),
  axis.text.y = element_text(size=4),  	#color="#993333", face="bold", angle=45
  axis.title.y = element_text( size = 12, face = "bold" )) +
  #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_wrap(facets = ~variable, scales = 'free', 
  #nrow=8, 
  ncol=5
  )



