DNA or RNA accessions used for testing VirTaK performance can be easily searched for in the VirTaK database:

```{bash, eval=FALSE, echo=TRUE}
while read -r line; do grep -A1 "$line" /virtak_db/VMR_MSL38_v1_complete_genomes.fasta > "$line.fasta"; done < DNA_accessions.txt
```

