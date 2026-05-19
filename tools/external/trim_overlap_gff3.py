#!/usr/bin/env python3

def parse_gff3(file):
    """
    Parse a GFF3 file and yield features as dictionaries.
    """
    with open(file) as f:
        for line in f:
            if line.startswith("#"):
                continue  # Skip headers
            fields = line.strip().split("\t")
            if len(fields) == 9:  # Ensure valid GFF3 format
                yield {
                    "seqid": fields[0],
                    "source": fields[1],
                    "type": fields[2],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "attributes": fields[8],
                    "line": line.strip(),  # Store the full line for output
                }


def extract_gene_records(file):
    """
    Extract gene features and their child features (e.g., mRNA, CDS, exon) from a GFF3 file.
    Returns a dictionary where keys are tuples of gene coordinates and values are lists of lines (gene + children).
    """
    genes = {}
    current_gene_key = None
    current_gene_lines = []

    for feature in parse_gff3(file):
        if feature["type"] == "gene":
            # Save the previous gene and its children
            if current_gene_key:
                genes[current_gene_key] = current_gene_lines

            # Start a new gene
            current_gene_key = (feature["seqid"], feature["start"], feature["end"], feature["strand"])
            current_gene_lines = [feature["line"]]

        elif current_gene_key:
            # Add child features (e.g., mRNA, CDS) to the current gene
            current_gene_lines.append(feature["line"])

    # Add the last gene
    if current_gene_key:
        genes[current_gene_key] = current_gene_lines

    return genes


def is_overlapping_same_strand(gene1, gene2):
    """
    Check if two genes overlap and are on the same strand.
    Overlap is defined as sharing the same seqid, overlapping start/end positions,
    and being on the same strand.
    """
    return (
        gene1[0] == gene2[0] and  # seqid matches
        max(gene1[1], gene2[1]) <= min(gene1[2], gene2[2]) and  # start/end overlap
        gene1[3] == gene2[3]  # strand matches
    )


def compare_gff3(file1, file2, output):
    """
    Compare two GFF3 files and exclude overlapping genes (on the same strand) from file1.
    Write non-overlapping genes (including their children) to the output file.
    """
    print(f"Parsing {file1}...")
    genes_file1 = extract_gene_records(file1)
    print(f"Found {len(genes_file1)} gene records in {file1}.")

    print(f"Parsing {file2}...")
    genes_file2 = extract_gene_records(file2)
    print(f"Found {len(genes_file2)} gene records in {file2}.")

    print("Identifying non-overlapping genes...")
    non_overlapping_genes = {}
    
    for gene1_coords, gene1_lines in genes_file1.items():
        overlap_found = False
        for gene2_coords in genes_file2.keys():
            if is_overlapping_same_strand(gene1_coords, gene2_coords):  # Check for overlap on the same strand
                overlap_found = True
                break
        if not overlap_found:
            non_overlapping_genes[gene1_coords] = gene1_lines

    print(f"Found {len(non_overlapping_genes)} non-overlapping gene records.")

    print(f"Writing non-overlapping genes to {output}...")
    with open(output, "w") as out:
        out.write("##gff-version 3\n")
        for lines in non_overlapping_genes.values():
            out.write("\n".join(lines) + "\n\n")  # Write all lines for each unique gene

    print(f"Non-overlapping genes written to {output}.")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Compare two GFF3 files and exclude overlapping genes (on the same strand) from file1."
    )
    parser.add_argument("-f1", "--file1", required=True, help="Path to the first GFF3 file (query file).")
    parser.add_argument("-f2", "--file2", required=True, help="Path to the second GFF3 file (reference file).")
    parser.add_argument("-o", "--output", required=True, help="Path to save non-overlapping gene records.")

    args = parser.parse_args()
    
    compare_gff3(args.file1, args.file2, args.output)
