import pysam
import sys
from collections import defaultdict

def read_annotations(annotation_file):
    annotations = defaultdict(list)
    max_position = float('-inf')
    with open(annotation_file, 'r') as f:
        for line in f:
            if not line.startswith('#'):  # Ignore header lines
                fields = line.strip().split('\t')
                if fields[2] == 'gene':  # Only process gene features
                    start, end = int(fields[3]), int(fields[4])
                    max_position = max(max_position, end)
                    attributes = dict(item.split('=') for item in fields[8].split(';'))
                    gene_name = attributes.get('Name', '').strip('"')  # Strip quotation marks
                    locus_tag = attributes.get('locus_tag', '').strip('"')  # Strip quotation marks
                    for position in range(start, end+1):
                        annotations[position].append((gene_name, locus_tag))
    return annotations, max_position

def create_wig(bam_file, annotation_file):
    # Initialize a defaultdict
    coverage = defaultdict(int)

    # Read the annotations
    annotations, max_position = read_annotations(annotation_file)

    # Open the BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_reverse:
                # For reverse strand reads, use the end coordinate
                coverage[read.reference_end] += 1
            else:
                # For forward strand reads, use the start coordinate
                coverage[read.reference_start] += 1

    # Write the coverage to a WIG file
    with open(bam_file.replace('.bam', '.wig'), 'w') as wig:
        for position in range(0, max_position+1):
            genes = ', '.join([f'{gene_name} {locus_tag}' for gene_name, locus_tag in annotations[position]])
            wig.write(f"{position}\t{coverage.get(position, 0)}\t{genes}\n")

if __name__ == "__main__":
    create_wig(sys.argv[1], sys.argv[2])