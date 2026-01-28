#!/usr/bin/env python3
# v1 2025-12-28

"""
Usage:
reorder_jf_fasta_dump_by_counts.py input.fasta [output.fasta]

  from jellyfish k-mer counter, e.g.

~/sratoolkit.3.3.0-ubuntu64/bin/fastq-dump --split-files --gzip SRR25073705
gzip -dc SRR25073705_?.fastq.gz | jellyfish count -m 25 -s 2G -C -o SRR25073705.counts -t 4 /dev/fd/0
jellyfish dump -L 1000 SRR25073705.l25.counts > SRR25073705.l25.fasta

  example for SRA record, SRX20827684, 24M spots, 4.9G bases
  https://www.ncbi.nlm.nih.gov/sra/SRX20827684[accn]
"""

import sys
from Bio import SeqIO

def main(fasta_in, fasta_out=None):
    seq_dict = {}

    # Read FASTA and populate dict
    for record in SeqIO.parse(fasta_in, "fasta"):
        try:
            key = int(record.id)
        except ValueError:
            raise ValueError(f"FASTA header is not an integer: {record.id}")

        seq_dict.setdefault(key, []).append(str(record.seq))

    # Sort keys descending
    sorted_keys = sorted(seq_dict.keys(), reverse=True)

    # Output
    out_handle = open(fasta_out, "w") if fasta_out else sys.stdout
    seq_counter = 1

    for key in sorted_keys:
        for seq in seq_dict[key]:
            out_handle.write(f">{key}|{seq_counter}\n")
            out_handle.write(f"{seq}\n")
            seq_counter += 1

    if fasta_out:
        out_handle.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit( __doc__ )

    fasta_in = sys.argv[1]
    fasta_out = sys.argv[2] if len(sys.argv) > 2 else None

    main(fasta_in, fasta_out)

