			###############################################################################################
			      ############## A de novo Pipeline for the identification of lncRNA ##############
			###############################################################################################

source master_pipeline.source

########################################################################

### Pair .fq files (CH01)

########################################################################

$PAIRFQ_EXEC makepairs -f $acc1\_1 -r $acc1\_2 -fp $acc1\_paired_forward -rp $acc1\_paired_reverse -fs $acc1\_unpaired_forward -rs $acc1\_unpaired_reverse -s > pairfq_$acc2.log
$PAIRFQ_EXEC makepairs -f $acc2\_1 -r $acc2\_2 -fp $acc2\_paired_forward -rp $acc2\_paired_reverse -fs $acc2\_unpaired_forward -rs $acc2\_unpaired_reverse -s > pairfq_$acc1.log

acc1_PAIRED_FORWARD=$acc1\_paired_forward
acc1_PAIRED_REVERSE=$acc1\_paired_reverse
acc1_UNPAIRED_FORWARD=$acc1\_unpaired_forward
acc1_UNPAIRED_REVERSE=$acc1\_unpaired_reverse

acc2_PAIRED_FORWARD=$acc2\_paired_forward
acc2_PAIRED_REVERSE=$acc2\_paired_reverse
acc2_UNPAIRED_FORWARD=$acc2\_unpaired_forward
acc2_UNPAIRED_REVERSE=$acc2\_unpaired_reverse

########################################################################

### Assemble via Trinity (CH02)

########################################################################

$TRINITY_EXECUTABLE --seqType fq --max_memory 25G --SS_lib_type FR --CPU 1 --full_cleanup --output $acc1\_trinity --left $acc1_PAIRED_FORWARD,$acc1_UNPAIRED_FORWARD --right $acc1_PAIRED_REVERSE,$acc1_UNPAIRED_REVERSE > $acc1.Trinity.log
$TRINITY_EXECUTABLE --seqType fq --max_memory 25G --SS_lib_type FR --CPU 1 --full_cleanup --output $acc2\_trinity --left $acc2_PAIRED_FORWARD,$acc2_UNPAIRED_FORWARD --right $acc2_PAIRED_REVERSE,$acc2_UNPAIRED_REVERSE > $acc2.Trinity.log

sed -i "s/>TR/>$acc1.TR/g' $acc1\_trinity.Trinity.fasta
sed -i "s/>TR/>$acc2.TR/g' $acc2\_trinity.Trinity.fasta

########################################################################

### Cluster with cd-hit (CH03)

########################################################################

$CD_HIT_EXEC -i $acc1\_trinity.Trinity.fasta -o $acc1\_clust.Trinity.fasta -M 8000 -T 4
$CD_HIT_EXEC -i $acc2\_trinity.Trinity.fasta -o $acc2\_clust.Trinity.fasta -M 8000 -T 4

########################################################################

### Quantify with RSEM (CH04)

########################################################################

# Do quantification (CH04.1)

$TRINITY_DIR/util/align_and_estimate_abundance.pl --seqType fq --transcripts $acc1\_clust.Trinity.fasta --left $acc1_PAIRED_FORWARD,$acc1_UNPAIRED_FORWARD --right $acc1_PAIRED_REVERSE,$acc1_PAIRED_REVERSE --SS_lib_type FR --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir $acc1\_RSEM
$TRINITY_DIR/util/align_and_estimate_abundance.pl --seqType fq --transcripts $acc2\_clust.Trinity.fasta --left $acc2_PAIRED_FORWARD,$acc2_UNPAIRED_FORWARD --right $acc2_PAIRED_REVERSE,$acc2_PAIRED_REVERSE --SS_lib_type FR --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir $acc2\_RSEM

# Process quantifies transcripts (CH04.2)

cp $acc1\_RSEM/RSEM.isoforms.results ./$acc1\_RSEM.isoforms.results
cp $acc2\_RSEM/RSEM.isoforms.results ./$acc2\_RSEM.isoforms.results

awk ' $7 >= 1.50 ' $acc1\_RSEM.isoforms.results | cut -f1 > $acc1\_RSEM_filt
awk ' $7 >= 1.50 ' $acc2\_RSEM.isoforms.results | cut -f1 > $acc2\_RSEM_filt

# Generate FASTA files from filt list (CH04.3)

grabFasta $acc1\_RSEM_filt $acc1\_clust.Trinity.fasta $acc1\_filt_clust.Trinity.fasta
grabFasta $acc2\_RSEM_filt $acc2\_clust.Trinity.fasta $acc2\_filt_clust.Trinity.fasta

########################################################################

### Remove known proteins with Trinotate (CH05)

########################################################################

# Prepare ORF for BLAST, HMMER (CH05.1)

$TRANSDECODER_EXEC -t $acc1\_filt_clust.Trinity.fasta
$TRANSDECODER_EXEC -t $acc2\_filt_clust.Trinity.fasta

cp $acc1\_filt_clust.Trinity.fasta.transdecoder_dir/longest_orfs.pep $acc1\_filt_clust.Trinity.pep
cp $acc2\_filt_clust.Trinity.fasta.transdecoder_dir/longest_orfs.pep $acc2\_filt_clust.Trinity.pep

# Blast for dayz (CH05.2)

function sortBlast {
	sort -k1,1 -k12,12gr -k11,11g -k3,3gr $1 | sort -u -k1,1 --merge > SORTED_$1
}


queryx1=$acc1\_filt_clust.Trinity.fasta
queryx2=$acc2\_filt_clust.Trinity.fasta

queryp1=$acc1\_filt_clust.Trinity.pep
queryp2=$acc2\_filt_clust.Trinity.pep

makeblastdb -in $db1 -dbtype prot
makeblastdb -in $db2 -dbtype prot

blastx -query $queryx1 -db $db1 -outfmt 6 -out $acc1\_$db_short_name1.outfmt6x -num_threads 4  -evalue 1e-20
blastx -query $queryx1 -db $db2 -outfmt 6 -out $acc1\_$db_short_name2.outfmt6x -num_threads 4  -evalue 1e-20

blastp -query $queryp1 -db $db1 -outfmt 6 -out $acc1\_$db_short_name1.outfmt6p -num_threads 4  -evalue 1e-20
blastp -query $queryp1 -db $db2 -outfmt 6 -out $acc1\_$db_short_name2.outfmt6p -num_threads 4  -evalue 1e-20

blastx -query $queryx2 -db $db1 -outfmt 6 -out $acc2\_$db_short_name1.outfmt6x -num_threads 4  -evalue 1e-20
blastx -query $queryx2 -db $db2 -outfmt 6 -out $acc2\_$db_short_name2.outfmt6x -num_threads 4  -evalue 1e-20

blastp -query $queryp2 -db $db1 -outfmt 6 -out $acc2\_$db_short_name1.outfmt6p -num_threads 4  -evalue 1e-20
blastp -query $queryp2 -db $db2 -outfmt 6 -out $acc2\_$db_short_name2.outfmt6p -num_threads 4  -evalue 1e-20

sortBlast $acc1\_$db_short_name1.outfmt6x
sortBlast $acc1\_$db_short_name2.outfmt6x
sortBlast $acc1\_$db_short_name1.outfmt6p
sortBlast $acc1\_$db_short_name2.outfmt6p
sortBlast $acc2\_$db_short_name1.outfmt6x
sortBlast $acc2\_$db_short_name2.outfmt6x
sortBlast $acc2\_$db_short_name1.outfmt6p
sortBlast $acc2\_$db_short_name2.outfmt6p

# HMMER (CH05.3)

$HMMSCAN_EXEC --cpu $threads_per_node --domtblout $acc1\_PFAM.out Pfam-A.hmm $queryp1 > $acc1\_pfam.log
$HMMSCAN_EXEC --cpu $threads_per_node --domtblout $acc2\_PFAM.out Pfam-A.hmm $queryp1 > $acc1\_pfam.log


# Generate new gene_trans_map (froom filt_clust) (CH05.4)

$TRINITY_DIR/util/support_scripts/get_Trinity_gene_to_trans_map.pl $acc1\_filt_clust.Trinity.fasta > $acc1\_filt_clust.Trinity.fasta.gene_trans_map
$TRINITY_DIR/util/support_scripts/get_Trinity_gene_to_trans_map.pl $acc2\_filt_clust.Trinity.fasta > $acc2\_filt_clust.Trinity.fasta.gene_trans_map

# Compile into Trinotate DB (CH05.5)

$TRINOTATE_EXEC $trinotate_db1 init --gene_trans_map $acc1\_filt_clust.Trinity.fasta.gene_trans_map --transcript_fasta $acc1\_filt_clust.Trinity.fasta --transdecoder_pep $acc1\_filt_clust.Trinity.pep
$TRINOTATE_EXEC $trinotate_db1 LOAD_swissprot_blastp SORTED_$acc1\_Sprot.outfmt6p
$TRINOTATE_EXEC $trinotate_db1 LOAD_TrEMBL_blastp SORTED_$acc1\_Sprot.outfmt6p
$TRINOTATE_EXEC $trinotate_db1 LOAD_pfam $acc1\_PFAM.out
$TRINOTATE_EXEC $trinotate_db1 LOAD_swissprot_blastx SORTED_$acc1\_Sprot.outfmt6x
$TRINOTATE_EXEC $trinotate_db1 LOAD_TrEMBL_blastx SORTED_$acc1\_Sprot.outfmt6x

$TRINOTATE_EXEC $trinotate_db1 report > $acc1\_trinotate_annotation_report

$TRINOTATE_EXEC $trinotate_db2 init --gene_trans_map $acc2\_filt_clust.Trinity.fasta.gene_trans_map --transcript_fasta $acc2\_filt_clust.Trinity.fasta --transdecoder_pep $acc2\_filt_clust.Trinity.pep
$TRINOTATE_EXEC $trinotate_db2 LOAD_swissprot_blastp SORTED_$acc2\_Sprot.outfmt6p
$TRINOTATE_EXEC $trinotate_db2 LOAD_TrEMBL_blastp SORTED_$acc2\_Sprot.outfmt6p
$TRINOTATE_EXEC $trinotate_db2 LOAD_pfam $acc2\_PFAM.out
$TRINOTATE_EXEC $trinotate_db2 LOAD_swissprot_blastx SORTED_$acc2\_Sprot.outfmt6x
$TRINOTATE_EXEC $trinotate_db2 LOAD_TrEMBL_blastx SORTED_$acc2\_Sprot.outfmt6x

$TRINOTATE_EXEC $trinotate_db2 report > $acc2\_trinotate_annotation_report

# Filter Trinotate via SQLite (CH05.6)

#. Filter $acc1 (CH05.6.1)

sed 's/|/:/g' $acc1_trinotate_annotation_report | sed 's/\t/|/g' | sed 's/ //g' > scratch

sqlite3 $acc1_trinotate_filter.sqlite "CREATE TABLE report(gene_id, transcript_id, sprot_Top_BLASTX_hit, TrEMBL_Top_BLASTX_hit, RNAMMER, prot_id, prot_coords, sprot_Top_BLASTP_hit, TrEMBL_Top_BLASTP_hit, Pfam, SignalP, TmHMM, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript, peptide);"
sqlite3 $acc1_trinotate_filter.sqlite ".import scratch report"

sqlite3 $acc1_trinotate_filter.sqlite "SELECT * from report
        WHERE sprot_Top_BLASTX_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTX_hit like '%Viridiplantae%'
        OR sprot_Top_BLASTP_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTP_hit like '%Viridiplantae%';" > $acc1.annotated_genes

sqlite3 $acc1_trinotate_filter.sqlite "SELECT transcript_id from report
        WHERE sprot_Top_BLASTX_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTX_hit like '%Viridiplantae%'
        OR sprot_Top_BLASTP_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTP_hit like '%Viridiplantae%';" > $acc1.annotated_genes.headers

sqlite3 $acc1_trinotate_filter.sqlite "SELECT * from report
        WHERE sprot_Top_BLASTX_hit not like '%Viridiplantae%'
        AND TrEMBL_Top_BLASTX_hit not like '%Viridiplantae%'
        AND sprot_Top_BLASTP_hit not like '%Viridiplantae%'
        AND TrEMBL_Top_BLASTP_hit not like '%Viridiplantae%';" > $acc1.not_viridiplantae_report

cp $acc1.not_viridiplantae_report  not_viridiplantae.scratch

sqlite3 $acc1_trinotate_filter.sqlite "CREATE TABLE not_viridiplantae(gene_id, transcript_id, sprot_Top_BLASTX_hit, TrEMBL_Top_BLASTX_hit, RNAMMER, prot_id, prot_coords, sprot_Top_BLASTP_hit, TrEMBL_Top_BLASTP_hit, Pfam, SignalP, TmHMM, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript,  peptide);"
sqlite3 $acc1_trinotate_filter.sqlite ".import not_viridiplantae.scratch not_viridiplantae"

sqlite3 $acc1_trinotate_filter.sqlite "SELECT * from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit not like '.'
        OR TrEMBL_Top_BLASTX_hit not like '.'
        OR sprot_Top_BLASTP_hit not like '.'
        OR TrEMBL_Top_BLASTP_hit not like '.';" > $acc1.contaminant_genes

sqlite3 $acc1_trinotate_filter.sqlite "SELECT transcript_id from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit not like '.'
        OR TrEMBL_Top_BLASTX_hit not like '.'
        OR sprot_Top_BLASTP_hit not like '.'
        OR TrEMBL_Top_BLASTP_hit not like '.';" > $acc1.contaminant_genes.headers      # didn't make into fasta

sqlite3 $acc1_trinotate_filter.sqlite "SELECT * from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit like '.'
        AND TrEMBL_Top_BLASTX_hit like '.'
        AND sprot_Top_BLASTP_hit like '.'
        AND TrEMBL_Top_BLASTP_hit like '.';" > $acc1.no_annotation

sqlite3 $acc1_trinotate_filter.sqlite "SELECT transcript_id from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit like '.'
        AND TrEMBL_Top_BLASTX_hit like '.'
        AND sprot_Top_BLASTP_hit like '.'
        AND TrEMBL_Top_BLASTP_hit like '.';" > $acc1.no_annotation.headers

rm scratch
rm not_viridiplantae.scratch

for i in annotated_genes contaminant_genes no_annotation; do
        sed -i 's/|/\t/g' $acc1.$i
        sed -i 's/:/|/g' $acc1.$i
done

for i in annotated_genes.headers contaminant_genes.headers no_annotation.headers; do
        sed -i 's/:/|/g' $acc1.$i
done

grabFasta $acc1.no_annotation.headers $acc1\_filt_clust.Trinity.fasta $acc1.no_annotation.headers.fasta
grabFasta $acc1.annotated_genes.headers $acc1\_filt_clust.Trinity.fasta $acc1.annotated_genes.headers.fasta

#. Filter $acc2 (CH05.6.2)

sed 's/|/:/g' $acc2_trinotate_annotation_report | sed 's/\t/|/g' | sed 's/ //g' > scratch

sqlite3 $acc2_trinotate_filter.sqlite "CREATE TABLE report(gene_id, transcript_id, sprot_Top_BLASTX_hit, TrEMBL_Top_BLASTX_hit, RNAMMER, prot_id, prot_coords, sprot_Top_BLASTP_hit, TrEMBL_Top_BLASTP_hit, Pfam, SignalP, TmHMM, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript, peptide);"
sqlite3 $acc2_trinotate_filter.sqlite ".import scratch report"

sqlite3 $acc2_trinotate_filter.sqlite "SELECT * from report
        WHERE sprot_Top_BLASTX_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTX_hit like '%Viridiplantae%'
        OR sprot_Top_BLASTP_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTP_hit like '%Viridiplantae%';" > $acc2.annotated_genes

sqlite3 $acc2_trinotate_filter.sqlite "SELECT transcript_id from report
        WHERE sprot_Top_BLASTX_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTX_hit like '%Viridiplantae%'
        OR sprot_Top_BLASTP_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTP_hit like '%Viridiplantae%';" > $acc2.annotated_genes.headers

sqlite3 $acc2_trinotate_filter.sqlite "SELECT * from report
        WHERE sprot_Top_BLASTX_hit not like '%Viridiplantae%'
        AND TrEMBL_Top_BLASTX_hit not like '%Viridiplantae%'
        AND sprot_Top_BLASTP_hit not like '%Viridiplantae%'
        AND TrEMBL_Top_BLASTP_hit not like '%Viridiplantae%';" > $acc2.not_viridiplantae_report

cp $acc2.not_viridiplantae_report  not_viridiplantae.scratch

sqlite3 $acc2_trinotate_filter.sqlite "CREATE TABLE not_viridiplantae(gene_id, transcript_id, sprot_Top_BLASTX_hit, TrEMBL_Top_BLASTX_hit, RNAMMER, prot_id, prot_coords, sprot_Top_BLASTP_hit, TrEMBL_Top_BLASTP_hit, Pfam, SignalP, TmHMM, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript, peptide);"
sqlite3 $acc2_trinotate_filter.sqlite ".import not_viridiplantae.scratch not_viridiplantae"

sqlite3 $acc2_trinotate_filter.sqlite "SELECT * from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit not like '.'
        OR TrEMBL_Top_BLASTX_hit not like '.'
        OR sprot_Top_BLASTP_hit not like '.'
        OR TrEMBL_Top_BLASTP_hit not like '.';" > $acc2.contaminant_genes

sqlite3 $acc2_trinotate_filter.sqlite "SELECT transcript_id from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit not like '.'
        OR TrEMBL_Top_BLASTX_hit not like '.'
        OR sprot_Top_BLASTP_hit not like '.'
        OR TrEMBL_Top_BLASTP_hit not like '.';" > $acc2.contaminant_genes.headers

sqlite3 $acc2_trinotate_filter.sqlite "SELECT * from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit like '.'
        AND TrEMBL_Top_BLASTX_hit like '.'
        AND sprot_Top_BLASTP_hit like '.'
        AND TrEMBL_Top_BLASTP_hit like '.';" > $acc2.no_annotation

sqlite3 $acc2_trinotate_filter.sqlite "SELECT transcript_id from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit like '.'
        AND TrEMBL_Top_BLASTX_hit like '.'
        AND sprot_Top_BLASTP_hit like '.'
        AND TrEMBL_Top_BLASTP_hit like '.';" > $acc2.no_annotation.headers

rm scratch
rm not_viridiplantae.scratch

for i in annotated_genes contaminant_genes no_annotation; do
        sed -i 's/|/\t/g' $acc2.$i
        sed -i 's/:/|/g' $acc2.$i
done

for i in annotated_genes.headers contaminant_genes.headers no_annotation.headers; do
        sed -i 's/:/|/g' $acc2.$i
done

grabFasta $acc2.no_annotation.headers $acc2\_filt_clust.Trinity.fasta $acc2.no_annotation.headers.fasta
grabFasta $acc2.annotated_genes.headers $acc2\_filt_clust.Trinity.fasta $acc2.annotated_genes.headers.fasta

########################################################################

### Compare no_annotation fasta files (CH06)

########################################################################

# Determine long and short for BLAST (CH06.1)

acc1_len=`grep -c ">" $acc1.no_annotation.headers.fasta`
acc2_len=`grep -c ">" $acc2.no_annotation.headers.fasta`

if (( $acc1_len > $acc2_len ))

	then long_acc=$acc1 && short_acc=$acc2
	else long_acc=$acc2 && short_acc=$acc1

fi

# BLAST for commonalities (CH06.2)

makeblastdb -in $long_acc.no_annotation.headers.fasta -dbtype nucl
blastn -query $short_acc.no_annotation.headers.fasta -db $long_acc.no_annotation.headers.fasta -outfmt 6 -out $short_acc\_to_\$long_acc.outfmt6 -num_threads $threads_per_node
sortBlast $short_acc\_to_\$long_acc.outfmt6

# Filter Important hits (CHO6.3)

#. SUMMARIZE BLAST OUT (CHO06.3.1)

awk ' $4 >= 200 ' > $short_acc\_to_\$long_acc.outfmt6

#. GET HEADERS, OVERLAP_EVALUE (CHO06.3.2)

cut -f1 SORTED_$short_acc\_to_\$long_acc.outfmt6 > $short_acc.headers
cut -f2 SORTED_$short_acc\_to_\$long_acc.outfmt6 > $long_acc.headers
cut -f4,11 SORTED_$short_acc\_to_\$long_acc.outfmt6 > overlap_evalue

#. GET FASTA SEQS (CHO06.3.3)

grabFasta $short_acc.headers $short_acc.no_annotation.headers.fasta $short_acc.headers.fasta
grabFasta $long_acc.headers $long_acc.no_annotation.headers.fasta $long_acc.headers.fasta

#. GET LENGTHS (CHO06.3.4)

bioawk -c fastx '{ print $name, length($seq) }' $short_acc.headers.fasta > $short_acc.headers.len
bioawk -c fastx '{ print $name, length($seq) }' $long_acc.headers.fasta > $long_acc.headers.len

#. ASSEMBLE MASTER_LIST (CHO06.3.4)

paste $short_acc.headers $short_acc.headers.len $long_acc.headers $long_acc.headers.len overlap_evalue > master_almost
awk '{ print $5/$2 }' master_almost > percent_$short_acc
awk '{ print $5/$4 }' master_almost > percent_$long_acc
paste master_almost percent_$short_acc percent_$long_acc > master_$spec\_comparison

#. GENERATE master_longest_* (CH06.3.5)

awk ' $2 > $4 ' master_$spec\_comparison | cut -f1 > master_longest_$short_acc
awk ' $2 < $4 ' master_$spec\_comparison | cut -f3 > master_longest_$long_acc
awk ' $2 == $4 ' master_$spec\_comparison | cut -f1 > master_longest_same

grabFasta master_longest_$short_acc $short_acc.headers.fasta master_longest_$short_acc.tmp
grabFasta master_longest_$long_acc $long_acc.headers.fasta master_longest_$long_acc.tmp
grabFasta master_longest_same $long_acc.headers.fasta master_longest_same.tmp

cat *.tmp > master_longest_$spec.fasta
rm *.tmp

########################################################################

### Run Infernal to remove known ncRNAs (CH07)

########################################################################

# Run Infernal (CH07.1)

$CMSCAN_EXEC --tblout cmscan.tabout $RFAM_DIR/Rfam.cm master_longest_$spec.fasta > cmscan.out 2> cmscan.err

# Clean Infernal Output, filter using sqlite, get seqs (CH07.2)

awk ' $11 == “!” ' cmscan.tabout | awk '{ print $1 }' > filt_cmscan.headers
grep ">" master_longest_$spec.fasta | sed 's/>//g' > master_longest_$spec.headers

sort master_longest_$spec.headers | uniq | sed 's/|/:/g' > long.scratch
sort filt_cmscan.headers | uniq | sed 's/|/:/g' > short.scratch

sqlite3 filter_rfam.sqlite "CREATE TABLE long_list(transcript_id)"
sqlite3 filter_rfam.sqlite ".import long.scratch long_list"

sqlite3 filter_rfam.sqlite "CREATE TABLE short_list(transcript_id)"
sqlite3 filter_rfam.sqlite ".import short.scratch short_list"

sqlite3 filter_rfam.sqlite "SELECT transcript_id FROM long_list
	WHERE transcript_id NOT IN
		(SELECT transcript_id from short_list)" > $spec\_rfam_filt.headers

sed -i 's/:/|/g' $spec\_rfam_filt.headers

grabFasta $spec\_rfam_filt.headers master_longest_$spec.fasta $spec\_rfam_filt.headers.fasta

########################################################################

### Validate lncRNAs using CPC (CH08)

########################################################################

# Do CPC (CH08.1)

$CPC_EXEC $spec\_rfam_filt.headers.fasta $spec.out $CPC_HOME $spec

# Filter CPC output (CH08.2)

grep "nonc" $spec.out | cut -f1 > $spec\_master_lncRNA.headers

grabFasta $spec\_master_lncRNA.headers $spec\_rfam_filt.headers.fasta $spec\_master_lncRNA.headers.fasta
