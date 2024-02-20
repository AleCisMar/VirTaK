#!/bin/bash

# Script: get_domain_info.sh
# Description: Retrieve domain information from Pfam-A.hmm.dat file based on domain names.
# Usage: get_domain_info.sh <domain_file> <pfam_file> <pfam_go_file> <go_file>
# Example: get_domain_info.sh domains.txt ~/db/virtak_db/Pfam-A.hmm.dat

# Display help message
display_help() {
    echo "########################"
    echo "## get_domain_info.sh ##"
    echo "########################"
    echo "Description: 
    Retrieve domain information from Pfam-A.hmm.dat file based on domain names."
    echo "Usage: 
    get_domain_info.sh <domain_file> <pfam_file> <pfam_go_file> <go_file>"
    echo "Example: 
    get_domain_info.sh domains.txt ~/db/virtak_db/Pfam-A.hmm.dat ~/db/virtak_db/pfam_gomf_most_specific.txt ~/db/virtak_db/go.txt"
    exit 1
}

# Check for the -h (help) option
if [ "$#" -eq 1 ] && [ "$1" == "-h" ]; then
    display_help
fi

# Check the number of command-line arguments
if [ "$#" -ne 4 ]; then
    display_help
fi

# Process command-line arguments
domain_file="$1"
pfam_file="$2"
pfam_go_file="$3"
go_file="$4"

# Main loop
echo "Getting PFAM IDs..."
cat "$domain_file" | while read line; do
    grep -A2 "^#=GF ID   $line\$" "$pfam_file" | 
    sed 's/   /#/' | 
    cut -d'#' -f3 | 
    paste -sd',' -
done > 1.tmp

echo "Relating PFAM IDs to GO IDs..."
cat 1.tmp | while read line; do 
    domain=$(echo $line | cut -d',' -f1)
    id=$(echo $line | cut -d',' -f2 | cut -d'.' -f1)
    description=$(echo $line | cut -d',' -f3)
    if [[ $(grep "$id" "$pfam_go_file") ]]; then 
        gos=$(grep "$id" "$pfam_go_file" | cut -d';' -f1 | cut -d':' -f2 | tr "\n" ";" | sed 's/;$//');
        echo $domain,$id,$description,$gos
    else 
        echo $domain,$id,$description,-
    fi
done > 2.tmp

echo "Getting GO names..."
cat 2.tmp | while read line; do
    domain=$(echo $line | cut -d',' -f1)
    id=$(echo $line | cut -d',' -f2)
    description=$(echo $line | cut -d',' -f3)
    gos=$(echo $line | cut -d',' -f4)
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
    echo $domain,$id,$description,$gos,$joined_names
done > domain_info.csv
echo "DONE!"
echo "Domain info written to domain_info.csv"