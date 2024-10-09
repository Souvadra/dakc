#!/bin/bash

# iterate over the .fq files in the current directory
for file in ./*.fq; do
    # Extract the filename without the extension
    filename=$(basename "$file" .fq)
            
    # Create the output file path
    output_file="$filename.txt"
                        
    # remove first line, third line, fourth line, and repeat every 4th line
    awk 'NR%4==2' "$file" > "$output_file" 
done
