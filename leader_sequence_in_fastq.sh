#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 motifs.fasta reads.fastq.gz" >&2
  exit 1
fi

#leader_sequence_in_fastq.sh leader_from_reads.txt SRR25073705_2.fastq.gz > reads_w_leader.fastq
#Reads with motif: 64505


motifs_fa="$1"
fastq_gz="$2"

gzip -dc -- "$fastq_gz" \
| awk -v M="$motifs_fa" '
  BEGIN {
    # Read FASTA motifs (ignore headers, concatenate sequences)
    seq = ""
    while ((getline line < M) > 0) {
      sub(/\r$/, "", line)
      if (line ~ /^>/) {
        if (seq != "") motifs[++n] = seq
        seq = ""
      } else {
        gsub(/[[:space:]]/, "", line)
        seq = seq line
      }
    }
    if (seq != "") motifs[++n] = seq
    close(M)

    count = 0
  }

  # FASTQ records are 4 lines
  {
    rec[(NR-1)%4 + 1] = $0
    if (NR % 4 == 0) {
      readseq = rec[2]
      hit = 0
      for (i = 1; i <= n; i++) {
        if (index(readseq, motifs[i]) > 0) {
          hit = 1
          break
        }
      }
      if (hit) {
        print rec[1]
        print rec[2]
        print rec[3]
        print rec[4]
        count++
      }
    }
  }

  END {
    # report to stderr so stdout remains valid FASTQ
    print "Reads with motif:", count > "/dev/stderr"
  }
'

