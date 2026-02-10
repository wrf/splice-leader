# splice-leader
analysis and filtering of trans-spliced leader sequences in IsoSeq data

## Finding the motifs in long reads using a genome

Using [minimap2](https://github.com/lh3/minimap2), map the long reads to the genome
```
~/minimap2-2.30_x64-linux/minimap2 -a -x splice --secondary=no GCF_026151205.1_MBARI_Bmic_1.0_genomic.fna CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz | samtools sort - -o BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.bam
samtools flagstat BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.bam
```


### [get_read_skip_from_bam.py](https://github.com/wrf/splice-leader/blob/main/get_read_skip_from_bam.py)
I had used the original script [`get_read_skip_from_bam.py`](https://github.com/wrf/splice-leader/blob/main/get_read_skip_from_bam.py) in both the [Hormiphora genome](https://pmc.ncbi.nlm.nih.gov/articles/PMC8527503/) with IsoSeq reads and the [Aphrocallistes genome](https://pmc.ncbi.nlm.nih.gov/articles/PMC10282587/) with Nanopore reads. Leader sequences were detectable in both species, each from a different phylum.

The script requires either a `.sam` file, or [samtools](https://github.com/samtools/samtools) to convert the `.bam` file to text for the standard input of the script.

```
samtools view BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.bam | get_read_skip_from_bam.py -n 35 -N 45 - > BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.read_skip.tab
Rscript ~/git/genome-reannotations/bam_read_skip_histogram.R BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.read_skip.tab

```

## Finding motifs in long reads de novo
As the reads are in 2-line FASTA format, the first 40-45bp would either contain a leader, or would be random since they are sheared. Using a pipeline with `cut`, these are sent as standard input to [`leader_sequence_uniqc_to_fasta.py`](https://github.com/wrf/splice-leader/blob/main/leader_sequence_uniqc_to_fasta.py). This makes a fasta file that is used as option `-m` in the filtering script [`retain_transcripts_w_leader_motif.py`](https://github.com/wrf/splice-leader/blob/main/retain_transcripts_w_leader_motif.py).

```
gzip -dc CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz | grep -v ">" | cut -c 1-45 | sort | uniq -c | sort -nr | head -n 40 | leader_sequence_uniqc_to_fasta.py --prefix BolinopsisDeep-Iso-V4656
retain_transcripts_w_leader_motif.py -i CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz -m BolinopsisDeep_leader.fasta -o CTE_BolinopsisDeep-Iso-V4656-D10+s7.w_leader.fasta
```

The filtering retains only about 1/3 of the long reads. These, in theory, are complete mRNAs. Most of the others will be sheared randomly, as evident in the genome.

```
# Reading motifs from BolinopsisDeep_leader.fasta
# Counted 18 motifs from BolinopsisDeep_leader.fasta
# Reading sequences from CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz
# Read sequences from 270742, wrote 95311 (35.20%)
```


## Extracting motifs from short reads
In the above examples having a genome and long reads, the motifs were already identified. Here I examined whether those motifs are identifiable in short reads, as they should be present, but immediately removed by short read assemblers due to the high connectivity.

Kmer counts are found with [jellyfish](https://github.com/gmarcais/Jellyfish). These are counted, then exported as fasta sequences. The option `-L 10000` ignores those with a count below 10000, since in theory, all mRNAs from this species should have the leader sequence.

The data are from [SRR25073705](https://www.ncbi.nlm.nih.gov/sra/SRX20827684), from *Mnemiopsis leidyi* regeneration RNA-Seq of whole animal, to check against the genome.

```
~/sratoolkit.3.3.0-ubuntu64/bin/fastq-dump --split-files --gzip -v SRR25073705
gzip -dc SRR25073705_?.fastq.gz | jellyfish count -m 25 -s 2G -C -o SRR25073705.k25.counts -t 4 /dev/fd/0
jellyfish dump -L 10000 SRR25073705.k25.counts > SRR25073705.k25.L10000.fasta
reorder_jf_fasta_dump_by_counts.py SRR25073705.k25.L10000.fasta SRR25073705.k25.L10000.sort.fasta
```

Viewing with `cat SRR25073705.k25.L10000.sort.fasta`, the actual consensus leader is quite low in frequency. PolyA tails are much higher frequency in the kmers, and the other of the top 5 are all connected.

```
>610193|1
CGAAAAGGACGACGCGGCAAGAAAG
>601845|2
CCGAAAAGGACGACGCGGCAAGAAA
>589456|3
GCCGAAAAGGACGACGCGGCAAGAA
>515402|4
AAAAAAAAAAAAAAAAAAAAAAAAA
>462342|5
CGGCCGAAAAGGACGACGCGGCAAG
...
>65418|163
ACTACTATTATACAAATAATTTGAG
```




