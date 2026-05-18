#!/usr/bin/env python3

import argparse

def convert_blast_to_gff(input_file, output_file):
    """
    Converts BLAST tabular output (format 6) to GFF3 format.
    
    Args:
        input_file (str): Path to the BLAST tabular output file.
        output_file (str): Path to the output GFF3 file.
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            cols = line.strip().split("\t")
            query_id = cols[0]
            subject_id = cols[1]
            perc_identity = cols[2]
            start = int(cols[8])
            end = int(cols[9])
            evalue = cols[10]
            
            # Determine strand and ensure start < end
            strand = "+" if start < end else "-"
            start, end = min(start, end), max(start, end)
            
            # Write GFF3 formatted line
            outfile.write(f"{query_id}\tBLAST\texon\t{start}\t{end}\t{perc_identity}\t{strand}\t.\tID={subject_id};evalue={evalue}\n")

def main():
    """
    Main function to parse arguments and run the conversion.
    """
    parser = argparse.ArgumentParser(
        description="Convert BLAST tabular output (format 6) to GFF3 format for use in tools like AUGUSTUS."
    )
    parser.add_argument(
        "-i", "--input", 
        required=True, 
        help="Path to the BLAST tabular output file (format 6)."
    )
    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Path to the output GFF3 file."
    )
    
    args = parser.parse_args()
    
    # Run the conversion
    convert_blast_to_gff(args.input, args.output)
    print(f"Conversion complete. GFF3 file saved to: {args.output}")

if __name__ == "__main__":
    main()

