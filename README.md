# Metagenomic analyis of dental calculus samples from Peterborough, England
Main pipelines and command lines used for the analysis of Late Medieval ancient dental calculus samples from the leprosarium of St Leonard, Peterborough (England)
The bioinformatic analyses were done on the Galileo100 supercomputing cluster of Cineca, with the support of Elixir-Italy and the HPC@CINECA program.

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
We used a read length of 65 bp and a threshold of 50 reads. 
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
