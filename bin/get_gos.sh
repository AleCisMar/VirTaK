#!/bin/bash

# Script: get_gos.sh
# Description: Retrieve Gene Ontologies based on PFAM accessions in domain_info.tsv
# Usage: get_gos.sh <domain_file> <pfam_go_file> <go_file>
# Example: get_gos.sh domain_info.tsv ~/db/virtak_db/pfam_gomf_most_specific.txt ~/db/virtak_db/go.txt

# Display help message
display_help() {
    echo "################"
    echo "## get_gos.sh ##"
    echo "################"
    echo "Description: 
    Retrieve Gene Ontologies based on PFAM accessions in domain_info.tsv"
    echo "Usage: 
    get_gos.sh <domain_file> <pfam_go_file> <go_file>"
    echo "Example: 
    get_gos.sh domain_info.tsv ~/db/virtak_db/pfam_gomf_most_specific.txt ~/db/virtak_db/go.txt"
    exit 1
}

# Check for the -h (help) option
if [ "$#" -eq 1 ] && [ "$1" == "-h" ]; then
    display_help
fi

# Check the number of command-line arguments
if [ "$#" -ne 3 ]; then
    display_help
fi

# Process command-line arguments
domain_file="$1"
pfam_go_file="$2"
go_file="$3"

# Main loop
echo "Getting PFAM IDs..."
cut -f1,2 "$domain_file" | sed '1d' > 1.tmp

echo "Relating PFAM IDs to GO IDs..."
cat 1.tmp | while read line; do 
    id=$(echo $line | cut -d' ' -f1)
    acc=$(echo $line | cut -d' ' -f2)
    if [[ $(grep "$acc" "$pfam_go_file") ]]; then 
        gos=$(grep "$acc" "$pfam_go_file" | cut -d';' -f1 | cut -d':' -f2 | tr "\n" ";" | sed 's/;$//');
        echo $id,$acc,$gos
    else 
        echo $id,$acc,-
    fi
done > 2.tmp

echo "Getting GO names..."
cat 2.tmp | while read line; do
    id=$(echo $line | cut -d',' -f1)
    acc=$(echo $line | cut -d',' -f2)
    gos=$(echo $line | cut -d',' -f3)
    names=()
    # Use process substitution to avoid subshell
    while IFS= read -r row; do
        if [[ -n $(grep -w "id: GO:$row" "$go_file") ]]; then
            name=$(grep -w -A1 "id: GO:$row" "$go_file" | grep name | sed 's/name: //')
            names+=("$name")
        else
            names+=("$row")
        fi
    done < <(echo $gos | tr ";" "\n")
    # Join the array elements with ';'
    joined_names=$(IFS=";"; echo "${names[*]}")
    echo -e "$id\t$acc\t$gos\t$joined_names"
done > domain_gos.tsv
rm 1.tmp 2.tmp
echo "DONE!"
echo "Domain Gene Ontologies written to domain_gos.tsv"