# splice-leader
analysis and filtering of trans-spliced leader sequences in IsoSeq data

## Finding the motifs in long reads using a genome
Using [minimap2](https://github.com/lh3/minimap2), map the long reads to the genome:

```
~/minimap2-2.30_x64-linux/minimap2 -a -x splice --secondary=no GCF_026151205.1_MBARI_Bmic_1.0_genomic.fna CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz | samtools sort - -o BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.bam
samtools flagstat BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.bam
```

or for the [Hcal_v1]() genome:

```
~/minimap2-2.30_x64-linux/minimap2 -a -x splice --secondary=no UCSC_Hcal_v1.fa SRR10403849.fasta.gz | samtools sort - -o SRR10403849_vs_UCSC_Hcal_v1.bam
[M::mm_idx_gen::1.936*1.29] collected minimizers
[M::mm_idx_gen::2.660*1.75] sorted minimizers
[M::main::2.660*1.75] loaded/built the index for 45 target sequence(s)
[M::mm_mapopt_update::2.892*1.69] mid_occ = 189
[M::mm_idx_stat] kmer size: 15; skip: 5; is_hpc: 0; #seq: 45
[M::mm_idx_stat::3.020*1.66] distinct minimizers: 22892626 (83.21% are singletons); average occurrences: 1.655; average spacing: 2.922; total length: 110691255
[M::worker_pipeline::705.371*2.99] mapped 210618 sequences
[M::worker_pipeline::1430.972*3.00] mapped 203948 sequences
[M::worker_pipeline::2169.129*3.00] mapped 201757 sequences
[M::worker_pipeline::2915.309*3.00] mapped 201680 sequences
[M::worker_pipeline::3643.484*3.00] mapped 201987 sequences
[M::worker_pipeline::4377.691*3.01] mapped 201348 sequences
[M::worker_pipeline::5111.227*3.01] mapped 201930 sequences
[M::worker_pipeline::5845.491*3.01] mapped 203233 sequences
[M::worker_pipeline::6608.030*3.01] mapped 204604 sequences
[M::worker_pipeline::7376.407*3.01] mapped 204526 sequences
[M::worker_pipeline::8134.723*3.01] mapped 205346 sequences
[M::worker_pipeline::8876.873*3.01] mapped 213068 sequences
[M::worker_pipeline::8900.834*3.01] mapped 8607 sequences
[M::main] Version: 2.30-r1287
[M::main] CMD: /home/wrf/minimap2-2.30_x64-linux/minimap2 -a -x splice --secondary=no UCSC_Hcal_v1.fa SRR10403849.fasta.gz
[M::main] Real time: 8900.890 sec; CPU: 26751.496 sec; Peak RSS: 4.110 GB
[bam_sort_core] merging from 12 files and 1 in-memory blocks...

~/minimap2-2.30_x64-linux/minimap2 -a -x splice --secondary=no UCSC_Hcal_v1.fa SRR10403849.w_leader.fasta | samtools sort - -o SRR10403849_vs_UCSC_Hcal_v1.w_leader.bam
[M::mm_idx_gen::1.945*1.29] collected minimizers
[M::mm_idx_gen::2.662*1.75] sorted minimizers
[M::main::2.662*1.75] loaded/built the index for 45 target sequence(s)
[M::mm_mapopt_update::2.893*1.69] mid_occ = 189
[M::mm_idx_stat] kmer size: 15; skip: 5; is_hpc: 0; #seq: 45
[M::mm_idx_stat::3.023*1.66] distinct minimizers: 22892626 (83.21% are singletons); average occurrences: 1.655; average spacing: 2.922; total length: 110691255
[M::worker_pipeline::826.629*2.98] mapped 209095 sequences
[M::worker_pipeline::1735.226*2.99] mapped 203100 sequences
[M::worker_pipeline::2657.437*2.99] mapped 202617 sequences
[M::worker_pipeline::3563.101*2.99] mapped 202403 sequences
[M::worker_pipeline::4426.804*2.99] mapped 203632 sequences
[M::worker_pipeline::5312.028*2.99] mapped 205295 sequences
[M::worker_pipeline::6227.500*2.99] mapped 205900 sequences
[M::worker_pipeline::6924.467*2.99] mapped 159506 sequences
[M::main] Version: 2.30-r1287
[M::main] CMD: /home/wrf/minimap2-2.30_x64-linux/minimap2 -a -x splice --secondary=no UCSC_Hcal_v1.fa SRR10403849.w_leader.fasta
[M::main] Real time: 6924.555 sec; CPU: 20721.943 sec; Peak RSS: 3.728 GB
[bam_sort_core] merging from 7 files and 1 in-memory blocks...

```

### [get_read_skip_from_bam.py](https://github.com/wrf/splice-leader/blob/main/get_read_skip_from_bam.py)
This script [`get_read_skip_from_bam.py`](https://github.com/wrf/splice-leader/blob/main/get_read_skip_from_bam.py) uses the `.sam` standard input and extracts any *S* from both the beginning and the end of each sequenced that mapped to the genome.

The leader sequence is evident in the `.bam` file because many of the reads are soft-masked at one end. This is written in the CIGAR string with the letter *S*. That is, a large fraction of the reads either *start* or *end* with a soft-masked portion of somewhere from 35-45bp. Several examples of CIGAR strings from [`Bmic_1.0`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_026151205.1/) are shown below, starting with *42S* (since the gene is on the `+` strand). One of them has *45S*, since it appears that there are multiple leader sequences, and any given transcript/gene randomly uses one. From the current data, it is not known if there is a preference of some genes for some leader sequences.

```
42S12M5D91M1I134M265N123M4076N177M4114N307M588N342M196N204M144N488M468N597M3I78M106N46M1I302M
42S13M5D87M1I137M265N123M4076N177M4114N307M588N342M196N204M144N488M468N675M106N305M1S
42S12M5D91M1I134M265N123M4076N177M4114N307M588N342M196N204M144N488M468N597M3I78M106N46M1I251M1S
45S12M5D91M1I134M265N123M4076N177M4114N307M588N342M196N204M144N488M468N597M3I78M106N46M1I302M
```

or all ending with *42S* (`-` strand):

```
84M1I179M4D109M396N191M16I46M1D2M1I13M1I202M1081N78M194N115M117N170M1117N95M156N124M1030N104M221N129M235N156M42S
84M1I12M1I167M4D109M374N213M16I46M1D2M1I13M1I202M1081N78M194N115M117N170M1117N12M1I83M156N124M1030N104M221N109M275N136M42S
84M1I179M4D109M374N213M16I46M1D2M1I13M1I202M1081N78M194N115M117N170M1117N95M156N124M1030N104M221N129M235N156M42S
83M1I95M120N63M374N227M16I46M1D2M1I13M1I202M1081N78M194N115M117N170M1117N95M156N124M1030N104M221N129M235N156M42S
```

I had used the original script in both the [Hormiphora genome](https://pmc.ncbi.nlm.nih.gov/articles/PMC8527503/) with IsoSeq reads and the [Aphrocallistes genome](https://pmc.ncbi.nlm.nih.gov/articles/PMC10282587/) with Nanopore reads. Leader sequences were detectable in both species, each from a different phylum.

The script requires either a `.sam` file, or [samtools](https://github.com/samtools/samtools) to convert the `.bam` file to text for the standard input of the script.

```
samtools view BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.bam | get_read_skip_from_bam.py -n 35 -N 45 - > BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.read_skip.tab
Rscript ~/git/genome-reannotations/bam_read_skip_histogram.R BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.read_skip.tab
```

Redoing for Hormiphora, `Counted 2491938 lines for 2467676 reads; 24262 had no CIGAR string; 1391932 reads had either end matched to target leader length`

```
samtools view SRR10403849_vs_UCSC_Hcal_v1.bam | get_read_skip_from_bam.py -n 35 -N 45 - > SRR10403849_vs_UCSC_Hcal_v1.read_skip.tab
Rscript ~/git/splice-leader/bam_read_skip_histogram.R SRR10403849_vs_UCSC_Hcal_v1.read_skip.tab
cut -f 4 SRR10403849_vs_UCSC_Hcal_v1.read_skip.tab | sort | uniq -c | sort -nr | ~/git/splice-leader/leader_sequence_uniqc_to_fasta.py --prefix Hcal_SRR10403849 > SRR10403849_vs_UCSC_Hcal_v1.leaders.fasta
```


## Finding motifs in long reads de novo
From the above example in a genome, it was also clear that the reads themselves would contain the leader. This can be extracted even if one does not have a genome.

As the reads are in 2-line FASTA format, and are in 5' to 3' direction, the first 40-45bp would either contain a leader, or would be random since they are sheared. 

### [leader_sequence_uniqc_to_fasta.py](https://github.com/wrf/splice-leader/blob/main/leader_sequence_uniqc_to_fasta.py)
Using a pipeline with `cut`, these are sent as standard input to [`leader_sequence_uniqc_to_fasta.py`](https://github.com/wrf/splice-leader/blob/main/leader_sequence_uniqc_to_fasta.py).

```
gzip -dc CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz | grep -v ">" | cut -c 1-45 | sort | uniq -c | sort -nr | head -n 40 | leader_sequence_uniqc_to_fasta.py --prefix BolinopsisDeep-Iso-V4656
```

```
gzip -dc SRR10403849.fasta.gz | grep -v ">" | cut -c 1-45 | sort | uniq -c | sort -nr | head -n 50 | ~/git/splice-leader/leader_sequence_uniqc_to_fasta.py --prefix Hcal_SRR10403849
gzip -dc SRR10403849.fasta.gz | grep -v ">" | cut -c 1-40 | sort | uniq -c | sort -nr | head -n 40 | ~/git/splice-leader/leader_sequence_uniqc_to_fasta.py --prefix Hcal_SRR10403849
```

This makes a fasta file that is used as option `-m` in the filtering script [`retain_transcripts_w_leader_motif.py`](https://github.com/wrf/splice-leader/blob/main/retain_transcripts_w_leader_motif.py).

### [retain_transcripts_w_leader_motif.py](https://github.com/wrf/splice-leader/blob/main/retain_transcripts_w_leader_motif.py)
Once a list of leader sequences is obtained, these can be used to keep only the long reads with the leader. In the example, the filtering retains only about 1/3 of the long reads. These, in theory, are complete mRNAs. Most of the others will be sheared randomly, as is evident in the genome by the "cascade" pattern of the mapping. 

**IMPORTANT**: Full-length mRNAs will generally all have the *exact* same 5' starting position in the genome. Sheared mRNAs will have *random* starting positions, such that it is very unlikely that any two will have the same starting position for the same gene. These sheared mRNAs are **NOT** alternate start sites in some kind of alternative splicing scheme; that *does* happen, but alternate start sites will nearly always have more than one long read with the **exact same** starting position, and *still have* a soft-masked portion of 35-45bp. 

```
retain_transcripts_w_leader_motif.py -i CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz -m BolinopsisDeep_leader.fasta -o CTE_BolinopsisDeep-Iso-V4656-D10+s7.w_leader.fasta
# Reading motifs from BolinopsisDeep_leader.fasta
# Counted 18 motifs from BolinopsisDeep_leader.fasta
# Reading sequences from CTE_BolinopsisDeep-Iso-V4656-D10+s7.fasta.gz
# Read sequences from 270742, wrote 95311 (35.20%)
```


```
retain_transcripts_w_leader_motif.py -i SRR10403849.fasta.gz -m SRR10403849_vs_UCSC_Hcal_v1.leaders.fasta -o SRR10403849.w_leader.fasta
# Reading motifs from SRR10403849_vs_UCSC_Hcal_v1.leaders.fasta
# Counted 178 motifs from SRR10403849_vs_UCSC_Hcal_v1.leaders.fasta
# Reading sequences from SRR10403849.fasta.gz
# Read sequences from 2462652, wrote 1591548 (64.63%)
Hcal_SRR10403849|4|189393|rc	GGAGTTTCAAACTTTTCAACACTACTTTAAACAAATTAATTTG	487327
Hcal_SRR10403849|2|354689|rc	GGGAGTTTCAAACTTTTCAACACTACTTTAAACAAATTAATTTG	915176
...
```

## Extracting motifs from short reads
In the above examples having a genome and long reads, the motifs were easily identified. Here I examined whether those motifs are identifiable in short reads, as they should be present, but immediately removed by short read assemblers due to the high connectivity.

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


# datasets #

[Derelle 2010 Convergent origins and rapid evolution of spliced leader trans-splicing in Metazoa: Insights from the Ctenophora and Hydrozoa](https://pmc.ncbi.nlm.nih.gov/articles/PMC2844618/)
[Douris 2010 Evidence for Multiple Independent Origins of trans-Splicing in Metazoa](https://academic.oup.com/mbe/article/27/3/684/1002183)


