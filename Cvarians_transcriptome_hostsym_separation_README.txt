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

# Split the transcriptome into 40 chunks to parallelize and decrease computing time per chunk
splitFasta.pl target.screened.fasta 40

# Create list of commands for blasting each subset chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db ~/annotate/ncbi/nr/nr\.fasta -evalue 0\.00001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/' > blast
launcher_creator.py -j blast -n blast -t 6:00:00 -q shortq7 -e email@gmail.com
sbatch blast.slurm

perl TaxaOriginByBlast.pl target.screened.fasta ~/annotate/ncbi/nr /export/scratch/bstrehlo/Sample_STR667_Cliona/bin 0.00001 kingdom false no  Cli_aqu2.chunk.100.blast

#generate combined blast report from files with prefix Cli and suffix .blast

cat Cli*.blast > Allblast.br
