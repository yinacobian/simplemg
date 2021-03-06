#!/bin/bash

#bash mg_runall_V2.sh IDS.txt /home/acobian/cobian2018_CFRR/mg 40

#To run do: thisscrit.sh [IDS.txt] [path to main folder] [number of threads to use in the system]

# $1 is a list of IDS
# $2 is the path to the main folder
# $3 is the number of threads to use : drogon 40

#This script requires:
#A folder to run the analysis : main folder
#A folder inside main folder, called P00_raw
#The raw reas need to be inside the folder P00_raw
#A list of IDS, which is the name of the raw reads files, everything thta is before .fastq 
#make a list of IDS:
#ls | grep '.fastq' | cut -f 1 -d '.' | cut -f 1,2 -d '_' | uniq > IDS.txt


#metagenome: CFRR metagenomes, only R1 is processed. 
#main path: 
#reads: 

#mkdir $2/P01_prinseq_output
#mkdir $2/P02_map_univec
#mkdir $2/P03_map_HG
#mkdir $2/P05_viral
#mkdir $2/P06_blastn_nt
#mkdir $2/P05_frap_viral

#mkdir $2/P04_vf
#mkdir $2/P04_vf/temp

#mkdir $2/P04_abx
#mkdir $2/P04_abx/temp

#mkdir $2/P04_exo/
#mkdir $2/P04_exo/temp

#mkdir $2/P07_denovo

#mkdir $2/P08_coverage_plots

#1.- Quality control

#quality filtering single end : prinseq++
#cat $1 | xargs -I{fileID} sh -c "prinseq++ -fastq $2/P00_raw/{fileID}_R1.fastq -ns_max 0 -derep -lc_entropy=0.5 -trim_qual_right=15 -trim_qual_left=15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -min_len 30 -min_qual_mean 20  -rm_header -out_name $2/P01_prinseq_output/{fileID} -threads $3 -out_format 1"

#quality filtering pair end : prinseq++
#cat $1 | xargs -I{fileID} sh -c 'prinseq++ -fastq /home/acobian/DH08102018/mg/P00_raw/{fileID}_L007_R1_001.fastq.gz -fastq2 /home/acobian/DH08102018/mg/P00_raw/{fileID}_L007_R2_001.fastq.gz -lc_entropy=0.5 -trim_qual_right=15 -trim_qual_left=15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -min_len 30 -min_qual_mean 20  -rm_header -out_name /home/acobian/DH08102018/mg/P01_prinseq_output/{fileID} -threads 40 -out_format 1'
#quality filtering pair end, only one end: prinseq++
#cat $1 | xargs -I{fileID} sh -c 'prinseq++ -fastq /home/acobian/DH08102018/mg/P00_raw/{fileID}_R1.fastq -lc_entropy=0.5 -trim_qual_right=15 -trim_qual_left=15 -trim_qual_type mean -trim_qual_rule lt -trim_qual_window 2 -min_len 30 -min_qual_mean 20  -rm_header -out_name /home/acobian/DH08102018/mg/P01_prinseq_output/{fileID} -threads 40 -out_format 1'

#remove UniVec : smalp map UniVec
#For one of the ends only:
#cat $1 | xargs -I{fileID} sh -c "smalt map -y 0.8 -n $3 -f samsoft -o $2/P02_map_univec/map_univec_qc_{fileID}.samsoft /home/DATABASES/UniVec/univec $2/P01_prinseq_output/{fileID}_good_out.fasta"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/samsoft2fasta_nohits.pl $2/P02_map_univec/map_univec_qc_{fileID}.samsoft > $2/P02_map_univec/nounivec_{fileID}.fasta"
##cat $1 | xargs -I{fileID} sh -c 'rm map_univec_qc_{fileID}.samsoft'

#2.- Host removal 

#remove human genome abundant sequences
#cat $1 | xargs -I{fileID} sh -c "smalt map -y 0.5 -n $3 -f samsoft -o $2/P03_map_HG/map_hgabundant_univec_qc_{fileID}.samsoft /home/DATABASES/HumanGenome/AbundantSequences/hgabundant $2/P02_map_univec/nounivec_{fileID}.fasta"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/samsoft2fasta_nohits.pl $2/P03_map_HG/map_hgabundant_univec_qc_{fileID}.samsoft > $2/P03_map_HG/nohgabundant_{fileID}.fasta"
##cat $1 | xargs -I{fileID} sh -c 'rm /home/acobian/DH08102018/mg/P03_map_HG/map_hgabundant_univec_qc_{fileID}.samsoft'

#remove human genome unplaced sequences 
#cat $1 | xargs -I{fileID} sh -c "smalt map -y 0.5 -n $3 -f samsoft -o $2/P03_map_HG/map_hgunplacedrandom_{fileID}.samsoft /home/DATABASES/HumanGenome/Chromosomes/hgunplacedrandom $2/P03_map_HG/nohgabundant_{fileID}.fasta"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/samsoft2fasta_nohits.pl $2/P03_map_HG/map_hgunplacedrandom_{fileID}.samsoft > $2/P03_map_HG/nohgunplacedrandom_{fileID}.fasta"
##cat $1 | xargs -I{fileID} sh -c 'rm /home/acobian/DH08102018/mg/P03_map_HG/map_hgunplacedrandom_{fileID}.samsoft'

#remove human genome chromosomal sequences 
#cat $1 | xargs -I{fileID} sh -c "smalt map -y 0.5 -n $3 -f samsoft -o $2/P03_map_HG/map_hgchromosomes_{fileID}.samsoft /home/DATABASES/HumanGenome/Chromosomes/hgchromosomes $2/P03_map_HG/nohgunplacedrandom_{fileID}.fasta"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/samsoft2fasta_nohits.pl $2/P03_map_HG/map_hgchromosomes_{fileID}.samsoft > $2/P03_map_HG/nohgchromosomes_{fileID}.fasta"
##cat $1 | xargs -I{fileID} sh -c 'rm /home/acobian/DH08102018/mg/P03_map_HG/map_hgunplacedrandom_{fileID}.samsoft'

#cat $1 | xargs -I{fileID} sh -c "cat $2/P03_map_HG/nohgchromosomes_{fileID}.fasta > $2/P03_map_HG/polished_{fileID}.fasta"

#3.- Custom databases 

#diamond makedb --in PATRIC_VF.faa -d PATRIC_VF
#diamond:virulence factors
#cat $1 | xargs -I{fileID} sh -c "diamond blastx -d /home/DATABASES/PATRIC/PATRIC_VF -q $2/P03_map_HG/polished_{fileID}.fasta -a $2/P04_vf/{fileID}_vs_patric -p $3 -t $2/P04_vf/temp"
#cat $1 | xargs -I{fileID} sh -c "diamond view -a $2/P04_vf/{fileID}_vs_patric.daa -o $2/P04_vf/{fileID}_vs_patric.m8"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P04_vf/{fileID}_vs_patric.m8 > $2/P04_vf/besthit_{fileID}_vs_patric.m8"
#cat $1 | xargs -I{fileID} sh -c "cut -f 1,2 $2/P04_vf/besthit_{fileID}_vs_patric.m8 | sort | uniq | cut -f 2 | sort | uniq -c | sort -nr  | sed -e 's/^ *//' | tr '"' '"' '"'\t'"' > $2/P04_vf/hits_virulencefactors_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "cut -f 2 $2/P04_vf/hits_virulencefactors_{fileID}.tab | cut -d '"'|'"' -f 2 | xargs -I{ID2} grep {ID2} /home/DATABASES/PATRIC/PATRIC_VF.faa > $2/P04_vf/names_hits_virulencefactors_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "paste $2/P04_vf/hits_virulencefactors_{fileID}.tab $2/P04_vf/names_hits_virulencefactors_{fileID}.tab > $2/P04_vf/OUT_virulencefactors_hits_and_names_{fileID}.tab"

#diamond makedb --in CARD.faa -d CARD
#diamond:antibiotic resistance genes
#rm -r /home/acobian/DH08102018/mg/P04_custom_dbs/temp
#mkdir /home/acobian/DH08102018/mg/P04_custom_dbs/temp

#cat $1 | xargs -I{fileID} sh -c "diamond blastx -d /home/DATABASES/PATRIC/CARD -q $2/P03_map_HG/polished_{fileID}.fasta -a $2/P04_abx/{fileID}_vs_CARD -p $3 -t $2/P04_abx/temp"
#cat $1 | xargs -I{fileID} sh -c "diamond view -a $2/P04_abx/{fileID}_vs_CARD.daa -o $2/P04_abx/{fileID}_vs_CARD.m8"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P04_abx/{fileID}_vs_CARD.m8 > $2/P04_abx/besthit_{fileID}_vs_CARD.m8"
#cat $1 | xargs -I{fileID} sh -c "cut -f 1,2 $2/P04_abx/besthit_{fileID}_vs_CARD.m8 | sort | uniq | cut -f 2 | sort | uniq -c | sort -nr  | sed -e 's/^ *//' | tr '"' '"' '"'\t'"' > $2/P04_abx/hits_abx_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "cut -f 2 $2/P04_abx/hits_abx_{fileID}.tab | cut -d '"'|'"' -f 2 | xargs -I{ID2} grep {ID2} /home/DATABASES/PATRIC/CARD.faa > $2/P04_abx/names_hits_abx_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "paste $2/P04_abx/hits_abx_{fileID}.tab $2/P04_vf/names_hits_abx_{fileID}.tab > $2/P04_abx/OUT_abx_hits_and_names_{fileID}.tab"

#diamond:Human_pathogenic_bacterial_exotoxin
#rm -r /home/acobian/DH08102018/mg/P04_custom_dbs/temp
#mkdir /home/acobian/DH08102018/mg/P04_custom_dbs/temp

#cat $1 | xargs -I{fileID} sh -c "diamond blastx -d /home/DATABASES/CustomDBS/Human_pathogenic_bacterial_exotoxin -q $2/P03_map_HG/polished_{fileID}.fasta -a $2/P04_exo/{fileID}_vs_exo -p $3 -t $2/P04_exo/temp"
#cat $1 | xargs -I{fileID} sh -c "diamond view -a $2/P04_exo/{fileID}_vs_exo.daa -o $2/P04_exo/{fileID}_vs_exo.m8"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P04_exo/{fileID}_vs_exo.m8 > $2/P04_exo/besthit_{fileID}_vs_exo.m8"
#cat $1 | xargs -I{fileID} sh -c "cut -f 1,2 $2/P04_exo/besthit_{fileID}_vs_exo.m8 | sort | uniq | cut -f 2 | sort | uniq -c | sort -nr  | sed -e 's/^ *//' | tr '"' '"' '"'\t'"' > $2/P04_exo/hits_exo_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "cut -f 2 $2/P04_exo/hits_exo_{fileID}.tab | cut -d '"'|'"' -f 2 | xargs -I{ID2} grep {ID2} /home/DATABASES/CustomDBS/Human_pathogenic_bacterial_exotoxin.faa > $2/P04_exo/names_hits_exo_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "paste $2/P04_exo/hits_exo_{fileID}.tab $2/P04_exo/names_hits_exo_{fileID}.tab > $2/P04_exo/OUT_exo_hits_and_names_{fileID}.tab"

#4.- BLASTn vs NT 

#add taxnomy file to blast folder : update_blastdb.pl taxdb
#add the database containing folder to the .bashrc file in your home: export BLASTDB=/home3/acobian/tax_db
#blastN vs NT
#cat $1 | xargs -I{fileID} sh -c "blastn -query $2/P03_map_HG/polished_{fileID}.fasta -db /home/DATABASES/blast/nt/nt -out $2/P06_blastn_nt/vs_NT_{fileID}.blastn -evalue 0.1 -num_threads $3 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatchgapopen qstart qend sstart send evalue bitscore sskingdoms sscinames'"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P06_blastn_nt/vs_NT_{fileID}.blastn > $2/P06_blastn_nt/besthit_vs_NT_{fileID}.blastn"
#cat $1 | xargs -I{fileID} sh -c "grep 'Viruses'  $2/P06_blastn_nt/besthit_vs_NT_{fileID}.blastn | cut -f 12 | sort -nr | uniq -c | sort -nr | sed -e 's/^ *//' | tr '"' '"' '"'\t'"' > $2/P06_blastn_nt/count_viral_hits_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "grep 'Bacteria' $2/P06_blastn_nt/besthit_vs_NT_{fileID}.blastn | cut -f 12 | sort -nr | uniq -c | sort -nr |  sed -e 's/^ *//' | tr '"' '"' '"'\t'"' > $2/P06_blastn_nt/count_species_bacteria_hits_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "grep 'Eukaryota' $2/P06_blastn_nt/besthit_vs_NT_{fileID}.blastn | cut -f 12 | sort -nr | uniq -c | sort -nr |  sed -e 's/^ *//' | tr '"' '"' '"'\t'"' > $2/P06_blastn_nt/count_eukaryota_hits_{fileID}.tab"
#cat $1 | xargs -I{fileID} sh -c "grep 'Bacteria' $2/P06_blastn_nt/besthit_vs_NT_{fileID}.blastn | cut -f 12 | cut -f 1 -d ' ' | sort -nr | uniq -c | sort -nr | sed -e 's/^ *//' | tr '"' '"' '"'\t'"' > $2/P06_blastn_nt/count_genus_bacteria_hits_{fileID}.tab"

#5.- Coverage plots for most abundant genus 

#get the most abundant genome from each genus from NCBI: 
#cat $1 | xargs -I{fileID} sh -c "mkdir $2/P08_coverage_plots/bacteria_plots_{fileID}"
#cat $1 | xargs -I{fileID} sh -c "head -n 5 $2/P06_blastn_nt/count_genus_bacteria_hits_{fileID}.tab | cut -f 2 > $2/P08_coverage_plots/bacteria_plots_{fileID}/topgenus_{fileID}.txt"
#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/topgenus_{fileID}.txt | xargs -I {genus} sh -c "'"'"grep {genus} $2/P06_blastn_nt/besthit_vs_NT_{fileID}.blastn | cut -f 2 | sort -n | uniq -c | sort -nr | head -n 1 | cut -f 4 -d '|' "'"'" > $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt"
#cat $1 | xargs -I{fileID} sh -c "mkdir $2/P08_coverage_plots/bacteria_plots_{fileID}/DB"
#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I{genomeID} sh -c "'"'"curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={genomeID}&retnode=text&rettype=fasta' > $2/P08_coverage_plots/bacteria_plots_{fileID}/DB/{genomeID}.fasta"'"'""

###need to fix the names of the reference genome files: 
#curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KC470543&retnode=text&rettype=fasta' > $(curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KC470543&retnode=xml&rettype=native' | xmllint --xpath "string(//Org-ref_taxname)" - | tr "\ " "_" ).fasta

cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {id} sh -c "'"'"grep {id} /home/acobian/cobian2018_CFRR/mg/P08_coverage_plots/bacteria_plots_{fileID}/DB/{id}.fasta"'"'" | cut -c 2- |  cut -f 1 -d ' ' --complement > $2/P08_coverage_plots/bacteria_plots_{fileID}/names_topspecies_{fileID}.txt"
cat $1 | xargs -I{fileID} sh -c "paste $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt $2/P08_coverage_plots/bacteria_plots_{fileID}/names_topspecies_{fileID}.txt > $2/P08_coverage_plots/bacteria_plots_{fileID}/ok_id_names_topspecies_{fileID}.txt"

#cat ids_topspecies_CF01mgD8.txt | xargs -I {id} sh -c 'grep -T {id} /home/acobian/cobian2018_CFRR/mg/P08_coverage_plots/bacteria_plots_CF01mgD8//DB/{id}.fasta | cut -c 2- ' |  cut -f 1 -d ' ' --complement  > names
#paste ids_topspecies_CF01mgD8.txt names > ok_names 


#map and get mpileups:

#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {genomeID} sh -c "'"'" smalt index -k 10 -s 5 $2/P08_coverage_plots/bacteria_plots_{fileID}/DB/{genomeID} $2/P08_coverage_plots/bacteria_plots_{fileID}/DB/{genomeID}.fasta"'"'" "
#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {genomeID} sh -c "'"'" smalt map -y 0.9 -n $3 -f sam -o $2/P08_coverage_plots/bacteria_plots_{fileID}/map_{genomeID}_{fileID}.sam $2/P08_coverage_plots/bacteria_plots_{fileID}/DB/{genomeID} $2/P03_map_HG/polished_{fileID}.fasta"'"'" "
#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {genomeID} sh -c "'"'" samtools view -b -T $2/P08_coverage_plots/bacteria_plots_{fileID}/DB/{genomeID}.fasta $2/P08_coverage_plots/bacteria_plots_{fileID}/map_{genomeID}_{fileID}.sam -o $2/P08_coverage_plots/bacteria_plots_{fileID}/map_{genomeID}_{fileID}.bam "'"'" "
#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {genomeID} sh -c "'"'" samtools sort -O bam -T toto $2/P08_coverage_plots/bacteria_plots_{fileID}/map_{genomeID}_{fileID}.bam  -o $2/P08_coverage_plots/bacteria_plots_{fileID}/sorted_map_{genomeID}_{fileID}.bam "'"'" "
#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {genomeID} sh -c "'"'" samtools mpileup -d 10000000 -a --reference $2/P08_coverage_plots/bacteria_plots_{fileID}/DB/{genomeID}.fasta $2/P08_coverage_plots/bacteria_plots_{fileID}/sorted_map_{genomeID}_{fileID}.bam -o $2/P08_coverage_plots/bacteria_plots_{fileID}/mpileup_map_{genomeID}_{fileID}.tab "'"'" "
#cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {genomeID} sh -c "'"'" cut -f 1,4 $2/P08_coverage_plots/bacteria_plots_{fileID}/mpileup_map_{genomeID}_{fileID}.tab > $2/P08_coverage_plots/bacteria_plots_{fileID}/cov_mpileup_map_{genomeID}_{fileID}.tab "'"'" "



#Coverage plots, needs the script PlotCoverage.R
cat $1 | xargs -I{fileID} sh -c "cat $2/P08_coverage_plots/bacteria_plots_{fileID}/ids_topspecies_{fileID}.txt | xargs -I {genomeID} sh -c "'"'" Rscript /home/acobian/bin/PlotCoverage.R {genomeID} {fileID} $2/P08_coverage_plots/bacteria_plots_{fileID} "'"'" "


#6.- Denovo assembly and comparison to NT 

#cat $1 | xargs -I{fileID} sh -c "spades.py -s $2/P03_map_HG/polished_{fileID}.fasta -t $3 --only-assembler -o $2/P07_denovo/spades_{fileID}"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/removesmalls.pl 900 $2/P07_denovo/spades_{fileID}/contigs.fasta > $2/P07_denovo/more900_contigs_{fileID}.fasta"
#cat $1 | xargs -I{fileID} sh -c "blastn -query $2/P07_denovo/more900_contigs_{fileID}.fasta -db /home/DATABASES/blast/nt/nt -out $2/P07_denovo/vs_NT_more900_contigs_{fileID}.blastn -evalue 0.1 -num_threads $3 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident length mismatchgapopen qstart qend sstart send evalue bitscore sskingdoms sscinames'"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl $2/P07_denovo/vs_NT_more900_contigs_{fileID}.blastn > $2/P07_denovo/besthit_vs_NT_more900_contigs_{fileID}.blastn"

#cat $1 | xargs -I{fileID} sh -c "cut -f 1 $2/P07_denovo/besthit_vs_NT_more900_contigs_{fileID}.blastn > $2/P07_denovo/toremove_{fileID}.txt"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/remove_fasta_id.pl $2/P07_denovo/toremove_{fileID}.txt $2/P07_denovo/more900_contigs_{fileID}.fasta > $2/P07_denovo/nohits_more900_contigs_{fileID}.fasta"
#cat $1 | xargs -I{fileID} sh -c "perl /home/acobian/bin/removesmalls.pl 2000 $2/P07_denovo/nohits_more900_contigs_{fileID}.fasta > $2/P07_denovo/nohits_more2000_contigs_{fileID}.fasta"

### rename the contigs to include the sample ID. 
### submitt to RAST and PATRIC for annotation 
### get the CDS 


#7.- FRAP viral refseq


#FRAP viral refseq
#sh -c "perl /home/acobian/bin/MYSCRIPTS/jmf3.pl /home/DATABASES/RefSeq/viral/all_viral_genomic.fna $2/P03_polished $2/P05_frap_viral smalt"
# this one worked: 
#perl /home/acobian/bin/MYSCRIPTS/jmf3.pl /home/DATABASES/RefSeq/viral/all_viral_genomic.fna /home/acobian/cobian2018_CFRR/mg/P03_polished results  smalt


#diamond:Viral proteins

#mkdir /home/acobian/DH08102018/mg/P05_viral/temp

#cat $1 | xargs -I{fileID} sh -c 'diamond blastx -d /home/DATABASES/RefSeq/viral/all_viral_proteins -q /home/acobian/DH08102018/mg/P03_map_HG/polished_{fileID}.fasta -a /home/acobian/DH08102018/mg/P05_viral/{fileID}_vs_viral_proteins -p 20 -t /home/acobian/DH08102018/mg/P05_viral/temp'
#cat $1 | xargs -I{fileID} sh -c 'diamond view -a /home/acobian/DH08102018/mg/P05_viral/{fileID}_vs_viral_proteins.daa -o /home/acobian/DH08102018/mg/P05_viral/{fileID}_vs_viral_proteins.m8'
#cat $1 | xargs -I{fileID} sh -c 'perl /home/acobian/bin/MYSCRIPTS/besthitblast.pl /home/acobian/DH08102018/mg/P05_viral/{fileID}_vs_viral_proteins.m8 > /home/acobian/DH08102018/mg/P05_viral/besthit_{fileID}_vs_viral_proteins.m8'
#cat $1 | xargs -I{fileID} sh -c 'cut -f 1,2 /home/acobian/DH08102018/mg/P05_viral/besthit_{fileID}_vs_viral_proteins.m8  | sort | uniq | cut -f2 | sort | uniq -c | sort -nr  | sed -e "s/^ *//" | tr " " "\t"  > /home/acobian/DH08102018/mg/P05_viral/hits_viral_proteins_{fileID}.tab'
#cat $1 | xargs -I{fileID} sh -c 'cut -f 2 /home/acobian/DH08102018/mg/P05_viral/hits_viral_proteins_{fileID}.tab | cut -d '"'"'|'"'"' -f 2 | xargs -I{ID2} grep {ID2} /home/acobian/DB/DB_CF_basic/all_viral_proteins.faa > /home/acobian/DH08102018/mg/P05_viral/names_hits_viral_proteins_{fileID}.tab'
#cat $1 | xargs -I{fileID} sh -c 'paste /home/acobian/DH08102018/mg/P05_viral/hits_viral_proteins_{fileID}.tab /home/acobian/DH08102018/mg/P05_viral/names_hits_viral_proteins_{fileID}.tab > /home/acobian/DH08102018/mg/P05_viral/OUT_viral_proteins_hits_and_names_{fileID}.tab'

##manually count hits
#ls | grep 'nounivec' | xargs -I {} sh -c ' echo {} ; grep "'">"'" {} | wc -l '




