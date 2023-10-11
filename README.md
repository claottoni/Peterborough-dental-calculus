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

### Classification of the reads with Kraken2

