#!/bin/bash

acc=
db_name=
file=Trinotate_annotation_report

sed 's/|/:/g' $file | sed 's/\t/|/g' | sed 's/ //g' > scratch

sqlite3 $db_name "CREATE TABLE report(gene_id, transcript_id, sprot_Top_BLASTX_hit, TrEMBL_Top_BLASTX_hit, RNAMMER, prot_id, prot_coords, sprot_Top_BLASTP_hit, TrEMBL_Top_BLASTP_hit, Pfam, SignalP, TmHMM, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript, peptide);"
sqlite3 $db_name ".import scratch report"

sqlite3 $db_name "SELECT * from report
	WHERE sprot_Top_BLASTX_hit like '%Viridiplantae%'
	OR TrEMBL_Top_BLASTX_hit like '%Viridiplantae%'
	OR sprot_Top_BLASTP_hit like '%Viridiplantae%'
	OR TrEMBL_Top_BLASTP_hit like '%Viridiplantae%';" > $acc.annotated_genes

sqlite3 $db_name "SELECT transcript_id from report
        WHERE sprot_Top_BLASTX_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTX_hit like '%Viridiplantae%'
        OR sprot_Top_BLASTP_hit like '%Viridiplantae%'
        OR TrEMBL_Top_BLASTP_hit like '%Viridiplantae%';" > $acc.annotated_genes.headers

sqlite3 $db_name "SELECT * from report
        WHERE sprot_Top_BLASTX_hit not like '%Viridiplantae%'
        AND TrEMBL_Top_BLASTX_hit not like '%Viridiplantae%'
        AND sprot_Top_BLASTP_hit not like '%Viridiplantae%'
        AND TrEMBL_Top_BLASTP_hit not like '%Viridiplantae%';" > $acc.not_viridiplantae_report

cp $acc.not_viridiplantae_report  not_viridiplantae.scratch

sqlite3 $db_name "CREATE TABLE not_viridiplantae(gene_id, transcript_id, sprot_Top_BLASTX_hit, TrEMBL_Top_BLASTX_hit, RNAMMER, prot_id, prot_coords, sprot_Top_BLASTP_hit, TrEMBL_Top_BLASTP_hit, Pfam, SignalP, TmHMM, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript, peptide);"
sqlite3 $db_name ".import not_viridiplantae.scratch not_viridiplantae"

sqlite3 $db_name "SELECT * from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit not like '.'
        OR TrEMBL_Top_BLASTX_hit not like '.'
        OR sprot_Top_BLASTP_hit not like '.'
        OR TrEMBL_Top_BLASTP_hit not like '.';" > $acc.contaminant_genes

sqlite3 $db_name "SELECT transcript_id from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit not like '.'
        OR TrEMBL_Top_BLASTX_hit not like '.'
        OR sprot_Top_BLASTP_hit not like '.'
        OR TrEMBL_Top_BLASTP_hit not like '.';" > $acc.contaminant_genes.headers

sqlite3 $db_name "SELECT * from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit like '.'
        AND TrEMBL_Top_BLASTX_hit like '.'
        AND sprot_Top_BLASTP_hit like '.'
        AND TrEMBL_Top_BLASTP_hit like '.';" > $acc.no_annotation

sqlite3 $db_name "SELECT transcript_id from not_viridiplantae
        WHERE sprot_Top_BLASTX_hit like '.'
        AND TrEMBL_Top_BLASTX_hit like '.'
        AND sprot_Top_BLASTP_hit like '.'
        AND TrEMBL_Top_BLASTP_hit like '.';" > $acc.no_annotation.headers

rm scratch
rm not_viridiplantae.scratch

for i in annotated_genes contaminant_genes no_annotation; do
	sed -i 's/|/\t/g' $acc.$i
	sed -i 's/:/|/g' $acc.$i
done

for i in annotated_genes.headers contaminant_genes.headers no_annotation.headers; do
	sed -i 's/:/|/g' $acc.$i
done
