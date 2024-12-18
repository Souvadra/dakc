#!/bin/bash

# iterate over the .fq files in the current directory
for file in ./*.fq; do
    # Extract the filename without the extension
    filename=$(basename "$file" .fq)
            
    # Create the output file path
    output_file="$filename.txt"
                        
    # remove first line, third line, fourth line, and repeat every 4th line
    awk 'NR%4==2' "$file" > "$output_file" 

    # Find the maximum line length
    max_length=$(awk '{ if (length > max) max = length } END { print max }' $output_file)

    # Pad the shorter lines with "M" characters
    awk -v max_len="$max_length" '{ printf "%-"max_len"s\n", $0 }' $output_file | sed 's/ /M/g' > output.txt
    mv output.txt $output_file
done
