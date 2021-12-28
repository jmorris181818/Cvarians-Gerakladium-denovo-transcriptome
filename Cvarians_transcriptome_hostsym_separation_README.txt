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
## Determine most likely source for each target-matched sequence in assembly based on taxonomic ID of each sequence's best match in NCBI's non-redundant (nr) database

cd ~/bin
git clone https://github.com/ckenkel/Pseudo-nitzschia2bRAD.git
mv Pseudo-nitzschia2bRAD/taxfiles.pl .
chmod +x taxfiles.pl

# NOTE: nr database and taxdump files downloaded 28 December 2021
mkdir ~/annotate/ncbi
cd ~/annotate/ncbi

echo 'update_blastdb.pl --decompress nr --passive' > get_nr
launcher_creator.py -j get_nr -n get_nr -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch get_nr.slurm

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

cp ~/bin/taxfiles.pl .
mkdir DBdir
srun perl taxfiles.pl

# Split the target-aligned and no match transcriptomes into 40 chunks to parallelize and decrease computing time per chunk
splitFasta.pl target.screened.fasta 40
splitFasta.pl nomatch.screened.fasta 40

# Create list of commands for blasting each subset chunk
ls subset*target* | perl -pe 's{^(\S+)$}{blastx -query $1 -db ~/annotate/ncbi/nr -evalue 0\.00001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br}' > bl_target
launcher_creator.py -j bl_target -n bl_target -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch bl_target.slurm

ls subset*nomatch* | perl -pe 's{^(\S+)$}{blastx -query $1 -db ~/annotate/ncbi/nr -evalue 0\.00001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br}' > bl_nomatch
launcher_creator.py -j bl_nomatch -n bl_nomatch -q mediumq7 -t 24:00:00 -e studivanms@gmail.com
sbatch bl_nomatch.slurm




## Get reference transcriptomes for de novo assembly (from Strehlow et al. 2021 Coral Reefs)
cd ~/annotate
mkdir cliona
wget -O Corientalis.fasta https://www.dropbox.com/sh/f74oan1eu32urd0/AADDRhtZApoJXZH7Nuv4GUYua/Cliona%20orientalis/cliona_transcriptome_with_names.fasta

cd ~/annotate
mkdir symG
wget -O Gendoclionum.fasta https://www.dropbox.com/sh/f74oan1eu32urd0/AACJrbfXnp2f11XgaIjGPb8Aa/Symbiodinium%20endoclionum/symbio_transcriptome_with_names.fasta

# Now follow the script 'Cvarians_transcriptome_annotation_README' for generating annotation files for differential expression analysis
