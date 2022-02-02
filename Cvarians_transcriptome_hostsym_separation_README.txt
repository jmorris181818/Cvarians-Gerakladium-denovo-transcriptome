## De novo transcriptome assembly pipeline for the sponge Cliona varians, version December 27, 2021
# Adapted by Michael Studivan (studivanms@gmail.com) and John Morris (john.morris@noaa.gov) based on repos by Misha Matz (https://github.com/z0on/annotatingTranscriptomes.git), Eli Meyer (https://github.com/Eli-Meyer/sequence_utilities.git; https://github.com/Eli-Meyer/transcriptome_utilities.git), and  Brian Strehlow (https://github.com/bstrehlow/Transcriptome_assembly.git) for use on the FAU KoKo HPC


#------------------------------
## BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions – please make sure to read them before copy-pasting.

# log onto cluster
ssh mstudiva@koko-login.hpc.fau.edu


#------------------------------
## Identify the most likely origin of each sequence by comparison to a protein DB from a single close relative and one or more databases of likely contaminants
# NOTE: before running blast, must shorten fasta headers to avoid error messages in blast output - try this fix: 'sed -e 's/>* .*$//' original.fasta > truncated.fasta'

sed -e 's/>* .*$//' Cvarians.fasta > Cvarians_trunc.fasta

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

345923 sequences input.
119317 of these matched Aqueenslandica_trunc.fasta more closely than any contaminants.
74157 matched contaminants more closely than Aqueenslandica_trunc.fasta.
152449 matched none of the supplied DB (nomatch.screened.fasta).


#------------------------------
## Attempting to identify additional host/symbiont sequences in the no match assembly based on taxonomic ID of each sequence's best match in NCBI's nucleotide (nt) database

# NOTE: nt and taxdb databases downloaded 28 December 2021
mkdir ~/annotate/ncbi/nt
cd ~/annotate/ncbi/nt

echo 'update_blastdb.pl --decompress nt --passive' > get_nt
launcher_creator.py -j get_nt -n get_nt -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch get_nt.slurm

update_blastdb.pl taxdb
tar -xzf taxdb.tar.gz

# Now need to update your blastdb path to include the location of taxonomy database files
nano ~/.bashrc
# Add the following text:
export BLASTDB="$HOME/annotate/ncbi/nt:$BLASTDB"
# Save, then source to make the path active
source ~/.bashrc

# Split the no match assembly into 40 chunks to parallelize and decrease computing time per chunk
splitFasta.pl nomatch.screened.fasta 40

# Create list of commands for blasting each subset chunk
for i in `ls subset*nomatch*.fasta`; do echo blastn -query $i -db ~/annotate/ncbi/nt/nt -evalue 0.0001 -num_threads 4 -max_target_seqs 5 -outfmt "'6 qseqid sseqid evalue pident stitle staxids sscinames scomnames sblastnames sskingdoms salltitles stitle'" -out $i.br; done > bl_nomatch
launcher_creator.py -j bl_nomatch -n bl_nomatch -q mediumq7 -t 24:00:00 -e studivanms@gmail.com -N 10
sbatch bl_nomatch.slurm

# check blast progress
cat subset*.br | wc -l
# found matches in blast database for 29914 out of 152449 sequences

# generate combined blast report
cat subset*nomatch*.br > allblast.br

# scp the allblast.br file to your computer and run the taxonomizr.R script


#------------------------------
## Once back from R script, scp output .txt files to your HPC directory

# Use grep to generate a fasta file of the host/symbiont blast matches in the nomatch assembly
grep -w -A 1 -f nomatch_symbiont.txt nomatch.screened.fasta --no-group-separator > nomatch_symbiont.fasta

# The length should match the taxonomy matches in the taxonomizr.R script
cat nomatch_symbiont.fasta| grep '>' | wc -l
# 75

grep -w -A 1 -f nomatch_host.txt nomatch.screened.fasta --no-group-separator > nomatch_host.fasta
cat nomatch_host.fasta| grep '>' | wc -l
# 507

# Combine the host/symbiont nomatch assemblies with the original target/contam assemblies
cat Smicroadriaticum_trunc.screened.fasta nomatch_symbiont.fasta > Gerakladium.fasta

cat target.screened.fasta nomatch_host.fasta > Cvarians.fasta

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
## GC content with BBMap package

scp the transcriptomes to your computer and follow the BBMap installation instructions here: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/

sh bbmap/stats.sh in=Cvarians.fasta
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2706	0.2300	0.2288	0.2706	0.0000	0.0000	0.0000	0.4588	0.0569

sh bbmap/stats.sh in=Gerakladium.fasta
A	C	G	T	N	IUPAC	Other	GC	GC_stdev
0.2041	0.2978	0.2926	0.2056	0.0000	0.0000	0.0000	0.5904	0.0403


#------------------------------
## Transcriptome completeness with gVolante/BUSCO

# Upload the transcriptomes to gVolante here: https://gvolante.riken.jp/analysis.html
# Use a cutoff length of '1', sequence type of 'Coding/transcribed (nucleotide)', ortholog search pipeline of 'BUSCO v5', and ortholog set of 'Metazoa' for host and 'Eukaryota' for symbiont

Cvarians.fasta
-------------------------
Completeness Assessment Results:
	Total # of core genes queried:    954
	# of core genes detected
		Complete:    887 (92.98%)
		Complete + Partial:    904 (94.76%)
	# of missing core genes:    50 (5.24%)
	Average # of orthologs per core genes:    4.16
	% of detected core genes that have more than 1 ortholog:    83.54
	Scores in BUSCO format:    C:93.0%[S:15.3%,D:77.7%],F:1.8%,M:5.2%,n:954
-------------------------

Gerakladium.fasta
-------------------------
Completeness Assessment Results:
	Total # of core genes queried:    255
	# of core genes detected
		Complete:    47 (18.43%)
		Complete + Partial:    63 (24.71%)
	# of missing core genes:    192 (75.29%)
	Average # of orthologs per core genes:    1.17
	% of detected core genes that have more than 1 ortholog:    8.51
	Scores in BUSCO format:    C:18.5%[S:16.9%,D:1.6%],F:6.3%,M:75.2%,n:255
-------------------------


# Now follow the script 'Cvarians_transcriptome_annotation_README' for generating annotation files for differential expression analysis
