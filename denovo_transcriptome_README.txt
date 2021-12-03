## De novo transcriptome assembly  pipeline, version November 24, 2021
# Adapted by Michael Studivan (studivanms@gmail.com) and John Morris (john.morris@noaa.gov) based on repos by Misha Matz (https://github.com/z0on/annotatingTranscriptomes.git) and Brian Strehlow (https://github.com/bstrehlow/Transcriptome_assembly.git) for use on the FAU KoKo HPC


#------------------------------
## BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

# log onto cluster
ssh mstudiva@koko-login.hpc.fau.edu


#------------------------------
## Installing RNA-seq scripts and setting up the workspace

# switch to home directory
cd

# unless you have done it in the past, make directory called bin,
# all your scripts should go in there:
mkdir bin

# switch to bin:
cd bin

# clone github repositories
git clone https://github.com/mstudiva/tag-based_RNAseq.git
git clone https://github.com/jmorris181818/Cvarians-Gerakladium-denovo-transcriptome.git
git clone https://github.com/Eli-Meyer/transcriptome_utilities
git clone https://github.com/Eli-Meyer/sequence_utilities.git

# move files from subdirectories to bin/:
mv tag-based_RNAseq/* .
mv Cvarians-Gerakladium-denovo-transcriptome/* .
mv transcriptome_utilities/* .
mv sequence_utilities/* .

# remove the tag-based_RNAseq directory
rm -rf tag-based_RNAseq

chmod +x ~/bin/bs

## Installing Trinity
wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.13.2/trinityrnaseq-v2.13.2.FULL.tar.gz
tar -vxf trinityrnaseq-v2.13.2.FULL.tar.gz
cd trinityrnaseq-v2.13.2/
make install
make plugins
# test it
cd sample_data/test_Trinity_Assembly/
./runMe.sh

## Installing jellyfish
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-2.3.0.tar.gz
tar -vxf jellyfish-2.3.0.tar.gz
cd jellyfish-2.3.0
./configure
make
make install


#------------------------------
## Downloading files via ftp

sftp morris_6888@dnaseq2.igsp.duke.edu
get -r Morris2_6888_210518B6

## Get reference transcriptomes for de novo assembly (from Strehlow et al. 2021 Coral Reefs)
cd ~/annotate
mkdir cliona
wget -O Corientalis.fasta https://www.dropbox.com/sh/f74oan1eu32urd0/AADDRhtZApoJXZH7Nuv4GUYua/Cliona%20orientalis/cliona_transcriptome_with_names.fasta

cd ~/annotate
mkdir symG
wget -O Gendoclionum.fasta https://www.dropbox.com/sh/f74oan1eu32urd0/AACJrbfXnp2f11XgaIjGPb8Aa/Symbiodinium%20endoclionum/symbio_transcriptome_with_names.fasta

#------------------------------
## Unzipping reads with a launcher_creator script

# creating and launching a cluster job to unzip all files:
 ls *.gz | perl -pe 's/(\S+)/gunzip $1/' > gunz
 launcher_creator.py -j gunz -n gunz -q shortq7 -t 6:00:00 -e studivanms@gmail.com
 sbatch --mem=200GB gunz.slurm

# look at the reads:
# head -50 SampleName.fastq
# note that every read has four lines, the ID line starts with @A00

# this little one-liner will show sequence-only in file:
# head -100 SampleName.fastq | grep -E '^[NACGT]+$'

# to count the number of reads in all samples
echo "countreads.pl > countreads_raw.txt" > count
launcher_creator.py -j count -n count -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count.slurm

#------------------------------
## Trimming and quality filtering

# Create conda environment
# Uncomment and run below if you don't have a conda env. set up
# module load miniconda3-4.6.14-gcc-8.3.0-eenl5dj
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

conda create -n sctld cutadapt

## Removing adaptors and low quality reads
echo '#!/bin/bash' > trim.sh
echo 'conda activate sctld' >> trim.sh
for F in *.fastq; do
echo "cutadapt $F -b GGGGGGGG -b AGATCGG -q 15 -m 50 -o ${F/.fastq/}.trim" >>trim.sh;
done

# Does not work with launcher_creator, consider breaking up script and running multiple jobs
sbatch -o trim.o%j -e trim.e%j --mem=200GB trim.sh

# how the job is doing?
squeue -u mstudiva

# double check you have the same number of files as samples
ll *.trim | wc -l

# but did the trimming really work?
# Use the same one-liner as before on the trimmed file to see if it is different
# from the raw one that you looked at before:

# head -100 SampleName.fq | grep -E '^[NACGT]+$'

# head -100 SampleName.trim | grep -E '^[NACGT]+$'
# the long runs of base A should be gone

# to save time in case of issues, move the raw fastq files to backup directory
mv *.fastq ~/backup/

## to count the number of reads in trimmed samples
echo "countreads_trim.pl > countreads_trim.txt" > count_trim
launcher_creator.py -j count_trim -n count_trim -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_trim.slurm


#------------------------------
## Re-pairing, deduplicating, and assembling forward and reverse reads

# copy the re-pair.sh script from bin to current directory
cp ~/bin/re-pair.sh .
launcher_creator.py -j re-pair.sh -n re-pair -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch re-pair.slurm

cd ..
mkdir cleanReads
cd cleanReads
mv rawReads/Un_* .
mv rawReads/Un_* .

## De-duplicating paired reads (requires left [R1], right [R2], and unpaired read files)
cp ~/bin/dedup.sh .
launcher_creator.py -j dedup.sh -n dedup -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch dedup.slurm

## to count the number of reads in deduplicated samples
echo "countreads_dedup.pl > countreads_dedup.txt" > count_dedup
launcher_creator.py -j count_dedup -n count_dedup -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch count_dedup.slurm

## Assembling deduplicated reads (and renames)
cp ~/bin/assemble.sh .
launcher_creator.py -j assemble.sh -n assemble -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch assemble.slurm

# when it's done, move all the .dedup files to backup
mv *.dedup ../backup/


#------------------------------
## Transcriptome assembly using Trinity

# first, concatenate all forward and reverse reads into a single set of fastq files per species
echo "cat R1_Cliona* > R1_Cliona.fastq" > concat
echo "cat R2_Cliona* > R2_Cliona.fastq" >> concat
launcher_creator.py -j concat -n concat -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch concat.slurm

# create directories for each of the species in your annotate directory, and move the concatenated fastq files there
mkdir cliona

# Trinity
# the python module required by launcher_creator is messing up the steps involving numpy, so we'll do it the old school way
echo '#!/bin/bash' > trinity.sh
echo '#SBATCH --partition=mediumq7' >> trinity.sh
echo '#SBATCH -N 1' >> trinity.sh
echo '#SBATCH --exclusive' >> trinity.sh
echo '#SBATCH --mem-per-cpu=16000' >> trinity.sh
echo 'module load python3/3.7.7' >> trinity.sh
echo "~/bin/trinityrnaseq-v2.13.2/Trinity --seqType fq --left R1_Cliona.fastq --right R2_Cliona.fastq --CPU 20 --max_memory 100G --output cliona_trinity" >> trinity.sh
sbatch -o trinity.o%j -e trinity.e%j trinity.sh

# If any of the assemblies fail in the chrysalis step, find the output directory for each of the error files and delete them, or move them to your backup directory. They should look like this: "cliona_trinity/read_partitions/Fb_4/CBin_4670/c467359.trinity.reads.fa.out"
mv cliona_trinity/read_partitions/Fb_4/CBin_4670/c467359.trinity.reads.fa.out cliona_trinity/read_partitions/Fb_3/CBin_3088/c309109.trinity.reads.fa.out cliona_trinity/read_partitions/Fb_3/CBin_3690/c369317.trinity.reads.fa.out cliona_trinity/read_partitions/Fb_3/CBin_3701/c370414.trinity.reads.fa.out temp_backup/.

mv cliona_trinity.Trinity.fasta Cvarians.fasta
echo "seq_stats.pl Cvarians.fasta > seqstats_Cvarians.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Cvarians.fasta
-------------------------
812356 sequences.
801 average length.
26720 maximum length.
182 minimum length.
N50 = 1343
650.6 Mb altogether (650572316 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## Removing contigs <400 bp, per Kitchen et al. (2015) doi: 10.1534/g3.115.020164
# Assemblies include many small contigs that are unlikely to provide significant matches, so for analyses based on sequence homology we consider only contigs ≥400 bp.
# Use removesmalls.pl to get rid of contigs < specified length

cd ~/bin
git clone  https://github.com/drtamermansour/p_asteroides.git
mv p_asteroides/scripts/removesmalls.pl .
chmod +x removesmalls.pl

cd ~/annotate/cliona/
echo "perl ~/bin/removesmalls.pl  400 Cvarians.fasta > Cvarians_l400.fasta" > smalls
launcher_creator.py -j smalls -n smalls -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch smalls.slurm

Cvarians_l400.fasta
-------------------------
436654 sequences.
1244 average length.
26720 maximum length.
400 minimum length.
N50 = 1717
543.4 Mb altogether (543381748 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## Removing reads matching to rRNA and mitoRNA (contamination)

# Downloading reference data
# Ribosomal RNA
# Go to https://www.arb-silva.de/search/ and search for your species; if it's not available, pick a related species
# Add it to your cart using the checkbox on the left, then select Download in the top right
# Choose FASTA without gaps, and tar.gz
# Cliona varians is available in this GitHub repo as 'arb-silva.de_2021-12-01_id1089989'
# Originally from Redmond et al. (2013) doi: 10.1093/icb/ict078
cp ~/bin/arb-silva.de_2021-12-01_id1089989.tgz .
tar -vxf arb-silva.de_2021-12-01_id1089989.tgz

# Mitochondrial RNA
# Cliona varians is available in this GitHub repo as 'Cvarians_mitoRNA'
# Originally from Plese et al. (2021) doi: 10.1016/j.ympev.2020.107011
cp ~/bin/Cvarians_mitoRNA.fasta .

echo "perl ~/bin/RemoveContamSeq_blast+.pl type=blastn score=45 reads=Cvarians_l400.fasta contam=rRNA,arb-silva.de_2021-12-01_id1089989_tax_silva.fasta contam=Mt,Cvarians_mitoRNA.fasta table=Cvarians_contamination.txt passed=Cvarians_clean.fasta" > contam
launcher_creator.py -j contam -n contam -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch contam.slurm

Cvarians_clean.fasta
-------------------------
436654 sequences in input file
0 sequences look like contaminants
        rRNA    0
        Mt	0
436654 sequences passed all tests
-------------------------


#------------------------------
## Identify the most likely origin of each sequence by comparison to a protein DB from a single close relative and one or more databases of likely contaminants
# NOTE: before running blast, must shorten fasta headers to avoid error messages in blast output - try this fix: 'sed -e 's/>* .*$//' original.fasta > truncated.fasta'

sed -e 's/>* .*$//' Cvarians_clean.fasta > Cvarians_trunc.fasta

# Only one complete sponge genome is available through Uniprot (Amphimedon queenslandica)
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000007879/UP000007879_400682.fasta.gz

# Closest zoox relative is Symbiodinium microadriaticum
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000186817/UP000186817_2951.fasta.gz

gunzip *.gz &
mv UP000007879_400682.fasta Aqueenslandica.fasta
mv UP000186817_2951.fasta Smicroadriaticum.fasta

cat Aqueenslandica.fasta| grep '>' | wc -l
# 43437 -matches number in Uniprot
cat Smicroadriaticum.fasta| grep '>' | wc -l
# 43269 -matches number in Uniprot

# Truncate databases
sed -e 's/>* .*$//' Aqueenslandica.fasta > Aqueenslandica_trunc.fasta
sed -e 's/>* .*$//' Smicroadriaticum.fasta > Smicroadriaticum_trunc.fasta

# Making a blast database for each reference
echo "makeblastdb -in Aqueenslandica_trunc.fasta -dbtype prot" >mdb
echo "makeblastdb -in Smicroadriaticum_trunc.fasta -dbtype prot" >>mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e email@gmail.com
sbatch mdb.slurm

# Running the modified perl script to blast against reference proteomes
echo "perl ~/bin/CompareContamSeq_blast+.pl -q Cvarians_clean.fasta -s 45 -t Aqueenslandica_trunc.fasta -c Smicroadriaticum_trunc.fasta" > origin
launcher_creator.py -j origin -n origin -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch origin.slurm

436654 sequences input.
137973 of these matched Aqueenslandica_trunc.fasta more closely than any contaminants.
81578 matched contaminants more closely than Aqueenslandica_trunc.fasta.
217103 matched none of the supplied DB (nomatch.screened.fasta).


#------------------------------
## Determine most likely source for each adig sequence in transcriptome based on taxonomic ID of each sequence's best match in NCBI's nr database
# NOTE: nr database and taxdump files downloaded 2 December 2021

mkdir ~/annotate/ncbi/nr
cd ~/annotate/ncbi/nr

srun wget 'ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz'
ls *.gz | perl -pe 's/(\S+)/tar -zxvf  $1/' > unz
launcher_creator.py -j unz -n unz -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch --mem=200GB unz.slurm


cd ~/annotate/ncbi
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz
