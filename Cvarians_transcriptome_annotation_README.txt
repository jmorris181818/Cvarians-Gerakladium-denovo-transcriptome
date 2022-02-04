# Transcriptome annotation, version February 3, 2022
# Adapted by Michael Studivan (studivanms@gmail.com) based on a repo by Misha Matz (https://github.com/z0on/annotatingTranscriptomes.git) for use on FAU's HPC (KoKo)
# for use in generating transcriptome annotation files for the sponge Cliona varians
# also includes the separation of reads associated with C. varians and algal symbiont Gerakladium transcriptomes


#------------------------------
## BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
## Script and workplace setup

# To install Bioperl in your bin directory, please follow these instructions:
cd bin
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/mstudiva/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/jmorris181818/Cvarians-Gerakladium-denovo-transcriptome.git
mv Cvarians-Gerakladium-denovo-transcriptome/* .
rm -rf Cvarians-Gerakladium-denovo-transcriptome/

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate

# If working with host/symbiont transcriptomes, remember to create separate directories for each
# And do each of the species-specific steps in the respective directories
mkdir cliona
mkdir symG


#------------------------------
## Getting transcriptome

# Cliona varians and Gerakladium transcriptomes (January 2022)
wget https://www.dropbox.com/s/mjgen7g4qo6s3s8/Cvarians.fasta
wget https://www.dropbox.com/s/oqxmpttgd5x3u5m/Gerakladium.fasta

# use the stream editor to find and replace all instances of "TRINITY" with "Cvarians" in the host transcriptome, and with "Gerakladium" in the symbiont
sed -i 's/TRINITY_DN/Cvarians/g' Cvarians.fasta
sed -i 's/TRINITY_DN/Gerakladium/g' Gerakladium.fasta

# transcriptome statistics
echo "seq_stats.pl Cvarians.fasta > seqstats_Cvarians.txt" > seq_stats
echo "seq_stats.pl Gerakladium.fasta > seqstats_Gerakladium.txt" >> seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

Cvarians.fasta
-------------------------
119824 sequences.
1901 average length.
26720 maximum length.
60 minimum length.
N50 = 2655
227.8 Mb altogether (227777238 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

Gerakladium.fasta
-------------------------
28670 sequences.
1375 average length.
21103 maximum length.
400 minimum length.
N50 = 1672
39.4 Mb altogether (39418006 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------


#------------------------------
## Uniprot annotations with blast

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzipping
gunzip uniprot_sprot.fasta.gz &

# indexing the fasta database
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch mdb.slurm

# splitting the transcriptomes into 200 chunks
splitFasta.pl Cvarians.fasta 200
splitFasta.pl Gerakladium.fasta 200

# blasting all 200 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script (119824 sequences for Cvarians, 28670 for Gerakladium)

# combining all blast results
cat subset*br > myblast.br
mv subset* ~/annotate/backup/

# generating isogroup designations for each contig/component (there can be multiple contigs per isogroup)
grep ">" Cvarians.fasta | perl -pe 's/>Cvarians(\d+)(\S+)/Cvarians$1$2\tCvarians$1/'>Cvarians_seq2iso.tab
cat Cvarians.fasta | perl -pe 's/>Cvarians(\d+)(\S+).+/>Cvarians$1$2 gene=Cvarians$1/'>Cvarians_iso.fasta

grep ">" Gerakladium.fasta | perl -pe 's/>Gerakladium(\d+)(\S+)\s.+/Gerakladium$1$2\tGerakladium$1/'>Gerakladium_seq2iso.tab
cat Gerakladium.fasta | perl -pe 's/>Gerakladium(\d+)(\S+).+/>Gerakladium$1$2 gene=Gerakladium$1/'>Gerakladium_iso.fasta


#-------------------------
## Extracting coding sequences and corresponding protein translations:
echo "perl ~/bin/CDS_extractor_v2.pl Cvarians_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch cddd

echo "perl ~/bin/CDS_extractor_v2.pl Gerakladium_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch cddd


#------------------------------
## GO annotations

# selecting the longest contig per isogroup (also renames using isogroups based on annotations):
fasta2SBH.pl Cvarians_iso_PRO.fas >Cvarians_out_PRO.fas

fasta2SBH.pl Gerakladium_iso_PRO.fas >Gerakladium_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# Cvarians status: go on web to http://eggnog-mapper.embl.de/job_status?jobname=MM_5_x1i5b_

# Gerakladium status: go on web to http://eggnog-mapper.embl.de/job_status?jobname=MM_hcpw1mqd

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_5_x1i5b_/out.emapper.annotations

wget http://eggnog-mapper.embl.de/MM_hcpw1mqd/out.emapper.annotations

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Cvarians_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Cvarians_iso2geneName.tab

awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Gerakladium_iso2go.tab
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Gerakladium_iso2geneName.tab

# using the stream editor to replace 'isogroup' with 'Gerakladium'
sed -i 's/isogroup/Gerakladium/g' Gerakladium_iso2go.tab
sed -i 's/isogroup/Gerakladium/g' Gerakladium_iso2geneName.tab


#------------------------------
## KOG annotations

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Cvarians_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Cvarians_iso2kogClass1.tab > Cvarians_iso2kogClass.tab

awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Gerakladium_iso2kogClass1.tab
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Gerakladium_iso2kogClass1.tab > Gerakladium_iso2kogClass.tab

# using the stream editor to replace 'isogroup' with 'Gerakladium'
sed -i 's/isogroup/Gerakladium/g' Gerakladium_iso2kogClass1.tab
sed -i 's/isogroup/Gerakladium/g' Gerakladium_iso2kogClass.tab


#------------------------------
## KEGG annotations

# selecting the longest contig per isogroup:
srun fasta2SBH.pl Cvarians_iso.fasta >Cvarians_4kegg.fasta

srun fasta2SBH.pl Gerakladium_iso.fasta >Gerakladium_4kegg.fasta

# scp *4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*4kegg.fasta .
# use web browser to submit Cvarians_4kegg.fasta file to KEGG's KAAS server http://www.genome.jp/kegg/kaas/
# select SBH method, upload nucleotide query

# Cvarians status: https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1643933241&key=6JojAnPf

# Gerakladium status: https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1643928758&key=vEIZtnjJ

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1643933241/query.ko

wget https://www.genome.jp/tools/kaas/files/dl/1643928758/query.ko

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Cvarians_iso2kegg.tab

cat query.ko | awk '{if ($2!="") print }' > Gerakladium_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
## Done! Transfer files

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/* .
