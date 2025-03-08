
```{}
while read -r line; do grep -A1 "$line" /virtak_db/VMR_MSL38_v1_complete_genomes.fasta > "$line.fasta"; done < DNA_accessions.txt
```

