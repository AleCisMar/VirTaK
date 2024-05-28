# Script used to extract some metadata for network representation
cp reordered_distance_matrix.csv matrix_network.csv
cut -d',' -f1 reordered_distance_matrix.csv | cut -d'|' -f1 > 1.tmp
cut -d',' -f1 reordered_distance_matrix.csv | cut -d'-' -f3,4 > 2.tmp
paste 1.tmp 2.tmp | sed 's/\t/|/' > 3.tmp
cat 3.tmp | while read line; do 
    if [[ $(echo $line | grep "-") ]]; then 
        #echo -e "Other\t$line";
        if [[ $(echo $line | grep "viridae") ]]; then
            echo $line | cut -d'|' -f1 > 1.temp
            echo $line | cut -d'-' -f2 > 2.temp
            paste 1.temp 2.temp | sed 's/\t/|/' > 3.temp
            a=$(cat 3.temp)
            echo -e "Other,$a"
            rm *.temp
        else
            b=$(echo $line | sed 's/-//');
            echo -e "Other,$b";
        fi;
    else 
        c=$(echo $line | cut -d'|' -f1);
        echo -e "Metagenome,$c"; 
    fi; 
done > 4.tmp
sed -i '1d' 4.tmp
sed -i '1s/^/Source,Virus\n/' 4.tmp
cut -d',' -f2 4.tmp > 1.temp
cut -d',' -f1 4.tmp > 2.temp
paste 1.temp 2.temp > metadata.txt
rm *.tmp *.temp
