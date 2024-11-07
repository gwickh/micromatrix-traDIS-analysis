#!/usr/bin/env bash
#author:    :Gregory Wickham
#date:      :20241107
#version    :1.0.0
#desc       :misc scripts
#usage		:bash misc.sh
#===========================================================================================================

#add filename to column of csv file
for csv_file in *.csv; do
    # Create a temporary file to store the updated content
    temp_file="temp_$csv_file"
    # Add the 'name' column header and include the filename in each row
    awk -v filename="$csv_file" 'BEGIN{FS=OFS=","} 
        NR==1 {$0 = $0 FS "name"}  # Add 'name' column to header row
        NR>1 {$0 = $0 FS filename}  # Add filename to each data row
        1' "$csv_file" > "$temp_file"
    # Replace the original CSV with the updated file
    mv "$temp_file" "$csv_file"
done
