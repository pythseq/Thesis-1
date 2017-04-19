#!/bin/bash

################################################################################################

# redefine any points to a specific file
# very inelegant. works.

#################################################################################################

db=iterative_rbh.sqlite

#################################################################################################

echo "#!/bin/bash"

echo "function sortBlast {
         sort -k1,1 -k12,12gr -k11,11g -k3,3gr \$1 | sort -u -k1,1 --merge > SORTED_\$1
}"

for i in {0..25}; do
echo "

## IMPORT ALL BLAST

sed -i 's/\t/|/g' iteration$i\_annotation2reference_query.outfmt6
sed -i 's/\t/|/g' iteration$i\_reference2annotation_query.outfmt6

sqlite3 $db 'CREATE TABLE iteration$i\_TOTAL_annotation2reference_query(query, target, percent_id, length, mismatches, gap_opens, query_start, query_end, target_start, target_end, evalue NUM, bit_score);'
sqlite3 $db 'CREATE TABLE iteration$i\_TOTAL_reference2annotation_query(query, target, percent_id, length, mismatches, gap_opens, query_start, query_end, target_start, target_end, evalue NUM, bit_score);'

sqlite3 $db '.import iteration$i\_annotation2reference_query.outfmt6 iteration$i\_TOTAL_annotation2reference_query'
sqlite3 $db '.import iteration$i\_reference2annotation_query.outfmt6 iteration$i\_TOTAL_reference2annotation_query'

sed -i 's/|/\t/g' iteration$i\_annotation2reference_query.outfmt6
sed -i 's/|/\t/g' iteration$i\_reference2annotation_query.outfmt6

## SORT BLAST

sortBlast iteration$i\_annotation2reference_query.outfmt6
sortBlast iteration$i\_reference2annotation_query.outfmt6

sed -i 's/\t/|/g' SORTED_iteration$i\_annotation2reference_query.outfmt6
sed -i 's/\t/|/g' SORTED_iteration$i\_reference2annotation_query.outfmt6

## LOAD BLAST

sqlite3 $db 'CREATE TABLE iteration$i\_annotation2reference_query(query, target, percent_id, length, mismatches, gap_opens, query_start, query_end, target_start, target_end, evalue NUM, bit_score);'
sqlite3 $db 'CREATE TABLE iteration$i\_reference2annotation_query(query, target, percent_id, length, mismatches, gap_opens, query_start, query_end, target_start, target_end, evalue NUM, bit_score);'

sqlite3 $db '.import SORTED_iteration$i\_annotation2reference_query.outfmt6 iteration$i\_annotation2reference_query'
sqlite3 $db '.import SORTED_iteration$i\_reference2annotation_query.outfmt6 iteration$i\_reference2annotation_query'

## RBH

sqlite3 $db 'SELECT iteration$i\_annotation2reference_query.query, iteration$i\_annotation2reference_query.target, iteration$i\_annotation2reference_query.evalue, iteration$i\_annotation2reference_query.bit_score, iteration$i\_reference2annotation_query.query, iteration$i\_reference2annotation_query.target, iteration$i\_reference2annotation_query.evalue, iteration$i\_reference2annotation_query.bit_score
	     FROM iteration$i\_annotation2reference_query
	     INNER JOIN iteration$i\_reference2annotation_query
	     ON iteration$i\_annotation2reference_query.target = iteration$i\_reference2annotation_query.query
	     WHERE iteration$i\_annotation2reference_query.query = iteration$i\_reference2annotation_query.target;' > iteration$i\_annotation_to_v2.1_transcriptome.rbh

sed -i 's/|/\t/g' iteration$i\_annotation_to_v2.1_transcriptome.rbh

## GET ANNOTATION, TARGET

ugh=$(( $i + 1 ))
#why=\`printf '%02d' \$ugh\`

cut -f1 iteration$i\_annotation_to_v2.1_transcriptome.rbh > iteration$i\_annotation
cut -f2 iteration$i\_annotation_to_v2.1_transcriptome.rbh > iteration$i\_reference

## FILTER BLAST

sqlite3 $db 'CREATE TABLE iteration$i\_annotation(id);'
sqlite3 $db 'CREATE TABLE iteration$i\_reference(id);'

sqlite3 $db '.import iteration$i\_annotation iteration$i\_annotation'
sqlite3 $db '.import iteration$i\_reference iteration$i\_reference'

sqlite3 $db 'SELECT * FROM iteration$i\_TOTAL_annotation2reference_query
             WHERE target NOT IN (SELECT * FROM iteration$i\_annotation)
             AND target NOT IN (SELECT * FROM iteration$i\_reference)
             AND query NOT IN (SELECT * FROM iteration$i\_reference)
             AND query NOT IN (SELECT * FROM iteration$i\_annotation);' > iteration\$ugh\_annotation2reference_query.outfmt6

sqlite3 $db 'SELECT * FROM iteration$i\_TOTAL_reference2annotation_query
             WHERE target NOT IN (SELECT * FROM iteration$i\_annotation)
             AND target NOT IN (SELECT * FROM iteration$i\_reference)
             AND query NOT IN (SELECT * FROM iteration$i\_reference)
             AND query NOT IN (SELECT * FROM iteration$i\_annotation);' > iteration\$ugh\_reference2annotation_query.outfmt6

sed -i 's/|/\t/g' iteration\$ugh\_annotation2reference_query.outfmt6
sed -i 's/|/\t/g' iteration\$ugh\_reference2annotation_query.outfmt6"

done
