#!/usr/bin/env python3

"""match_noleader_orthologs.py  last modified 2026-05-18

~/git/splice-leader/match_noleader_orthologs.py  --table1 ~/genomes/hormiphora/SRR10403849_vs_UCSC_Hcal_v1.w_leader.compared_coverage.tab  --table2 ~/genomes/bolinopsis/BolinopsisDeep-Iso-V4656_vs_MBARI_Bmic_1.0.w_leader.compared_coverage.tab  --gff1 ~/genomes/hormiphora/Hcv1av93.leader_removed.gff  --gff2 ~/genomes/bolinopsis/GCF_026151205.1_MBARI_Bmic_1.0_genomic.gff  --blast Hcv1av93_vs_Bmic_1.0.blastp.tab.gz  -o Hcal_v1_vs_Bmic_1.0.no_leader_sequence_ortholog_pairs.tsv
"""
#!/usr/bin/env python3

import sys
import argparse
import gzip

def open_text(path):
    """
    Transparently open plain text or gzipped files.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_attributes(attr_string):
    """
    Parse GFF3 attribute column into dict.
    """

    attrs = {}

    for item in attr_string.strip().split(";"):

        if "=" not in item:
            continue

        k, v = item.split("=", 1)
        attrs[k] = v

    return attrs


def normalize_id(x):
    """
    Normalize IDs for flexible matching.

    Handles:
      rna-XM_...
      cds-XP_...
      version stripping
      whitespace trimming
    """

    if x is None:
        return None

    x = x.strip()
    x = x.split()[0]

    for prefix in (
        "rna-",
        "cds-",
        "gene-",
        "transcript-",
        "mrna-",
    ):
        if x.startswith(prefix):
            x = x[len(prefix):]

    return x


def add_alias(mapping, alias, gene_id):
    """
    Add multiple alias forms to lookup table.
    """

    if not alias:
        return

    alias = alias.strip()

    mapping[alias] = gene_id

    norm = normalize_id(alias)

    if norm:
        mapping[norm] = gene_id

        if "." in norm:
            mapping[norm.split(".")[0]] = gene_id


def lookup_gene(blast_id, id_to_gene):
    """
    Flexible BLAST ID -> gene lookup.
    """

    if blast_id is None:
        return None

    raw = blast_id.strip().split()[0]

    candidates = [raw]

    norm = normalize_id(raw)

    if norm:
        candidates.append(norm)

        if "." in norm:
            candidates.append(norm.split(".")[0])

    for c in candidates:
        if c in id_to_gene:
            return id_to_gene[c]

    return None


def read_gene_table(path):
    """
    Read table:
      gene_id chrom start end

    Returns:
      dict gene_id -> first 4 columns
    """

    genes = {}

    with open_text(path) as f:

        for line in f:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 4:
                continue

            leader_flag = fields[11]
            if leader_flag !="nls":
                continue

            gene_id = fields[0]
            genes[gene_id] = fields[:4]


    return genes


def build_blast_id_to_gene_map(gff_file):
    """
    Build flexible mapping:
      BLAST-compatible ID -> gene ID

    Supports:
      transcript IDs
      protein IDs
      CDS IDs
      Name=
      Parent=
      accession-like IDs
    """

    transcript_to_gene = {}
    blast_id_to_gene = {}

    ###########################################################
    # PASS 1
    # transcript/mRNA -> gene
    ###########################################################

    with open_text(gff_file) as f:

        for line in f:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            feature = fields[2]

            if feature not in (
                "mRNA",
                "transcript",
                "lnc_RNA",
                "tRNA",
                "rRNA",
            ):
                continue

            attrs = parse_attributes(fields[8])

            tx_id = attrs.get("ID")
            parent_gene = attrs.get("Parent")

            if not tx_id or not parent_gene:
                continue

            gene_id = parent_gene.split(",")[0]

            transcript_to_gene[tx_id] = gene_id
            transcript_to_gene[normalize_id(tx_id)] = gene_id

            # Add all transcript aliases
            for key in (
                "ID",
                "Name",
                "transcript_id",
                "transcript",
                "locus_tag",
                "gene",
            ):
                add_alias(
                    blast_id_to_gene,
                    attrs.get(key),
                    gene_id
                )

            # Add gene aliases too
            add_alias(blast_id_to_gene, gene_id, gene_id)

    ###########################################################
    # PASS 2
    # CDS/protein -> gene
    ###########################################################

    with open_text(gff_file) as f:

        for line in f:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            feature = fields[2]

            if feature != "CDS":
                continue

            attrs = parse_attributes(fields[8])

            parent = attrs.get("Parent")

            if not parent:
                continue

            parent_tx = parent.split(",")[0]

            gene_id = (
                transcript_to_gene.get(parent_tx)
                or transcript_to_gene.get(normalize_id(parent_tx))
            )

            if not gene_id:
                continue

            # Add CDS/protein aliases
            for key in (
                "protein_id",
                "ID",
                "Name",
                "Parent",
                "Derives_from",
                "product_accession",
                "locus_tag",
                "gene",
            ):
                add_alias(
                    blast_id_to_gene,
                    attrs.get(key),
                    gene_id
                )

            # Explicitly add transcript parent
            add_alias(
                blast_id_to_gene,
                parent_tx,
                gene_id
            )

    return blast_id_to_gene


def read_best_hits(
    blast_file,
    query_id_to_gene,
    subject_id_to_gene,
    query_gene_set
):
    """
    Read BLAST outfmt6.

    Keeps:
      best subject hit per query gene

    Best hit:
      highest bitscore
    """

    best_hits = {}

    unmapped_queries = 0
    unmapped_subjects = 0

    with open_text(blast_file) as f:

        for line in f:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 12:
                continue

            qid = fields[0]
            sid = fields[1]

            bitscore = float(fields[11])

            qgene = lookup_gene(
                qid,
                query_id_to_gene
            )

            sgene = lookup_gene(
                sid,
                subject_id_to_gene
            )

            if qgene is None:
                unmapped_queries += 1
                continue

            if sgene is None:
                unmapped_subjects += 1
                continue

            if qgene not in query_gene_set:
                continue

            if (
                qgene not in best_hits
                or bitscore > best_hits[qgene]["bitscore"]
            ):

                best_hits[qgene] = {
                    "subject_gene": sgene,
                    "bitscore": bitscore,
                }

    print(f"Unmapped query proteins: {unmapped_queries}", file=sys.stderr)
    print(f"Unmapped subject proteins: {unmapped_subjects}", file=sys.stderr)

    return {
        qgene: x["subject_gene"]
        for qgene, x in best_hits.items()
    }


def read_gene_descriptions_from_gff(gff_file, default="NA"):
    """
    Extract gene feature description tags from GFF.

    Returns:
      dict gene_id -> description
    """

    descriptions = {}

    with open_text(gff_file) as f:

        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            feature = fields[2]

            if feature != "gene":
                continue

            attrs = parse_attributes(fields[8])

            gene_id = attrs.get("ID")

            if not gene_id:
                continue

            desc = attrs.get("description", default)

            descriptions[gene_id] = desc

            # also store normalized forms, in case bgene is LOC... not gene-LOC...
            norm_gene_id = normalize_id(gene_id)

            if norm_gene_id:
                descriptions[norm_gene_id] = desc

    return descriptions


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Identify orthologs between two datasets "
            "using BLAST and GFF-derived ID mappings."
        )
    )

    parser.add_argument(
        "--table1",
        required=True,
        help="Species A table: gene chrom start end"
    )

    parser.add_argument(
        "--table2",
        required=True,
        help="Species B table: gene chrom start end"
    )

    parser.add_argument(
        "--gff1",
        required=True,
        help="Species A GFF/GFF3"
    )

    parser.add_argument(
        "--gff2",
        required=True,
        help="Species B GFF/GFF3"
    )

    parser.add_argument(
        "--blast",
        required=True,
        help="BLAST outfmt6 file (gzipped or plain)"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV"
    )

    args = parser.parse_args()

    ###########################################################
    # Read coordinate tables
    ###########################################################

    print("Reading coordinate tables...", file=sys.stderr)

    table1 = read_gene_table(args.table1)
    table2 = read_gene_table(args.table2)

    print(f"Species A genes: {len(table1)}", file=sys.stderr)
    print(f"Species B genes: {len(table2)}", file=sys.stderr)

    ###########################################################
    # Build GFF mappings
    ###########################################################

    print("Building GFF ID mappings...", file=sys.stderr)

    idmap1 = build_blast_id_to_gene_map(args.gff1)
    idmap2 = build_blast_id_to_gene_map(args.gff2)

    print(f"Species A mapped IDs: {len(idmap1)}", file=sys.stderr)
    print(f"Species B mapped IDs: {len(idmap2)}", file=sys.stderr)

    print("Reading species B gene descriptions...")

    desc2 = read_gene_descriptions_from_gff( args.gff2, default="NA" )

    ###########################################################
    # Read BLAST
    ###########################################################

    print("Reading BLAST hits...", file=sys.stderr)

    best_hits = read_best_hits(
        args.blast,
        idmap1,
        idmap2,
        set(table1)
    )

    print(f"Best-hit orthologs found: {len(best_hits)}", file=sys.stderr)

    ###########################################################
    # Write output
    ###########################################################

    kept = 0

    with open(args.output, "w") as out:

        out.write(
            "\t".join([
                "A_gene_id",
                "A_chr",
                "A_start",
                "A_end",
                "B_gene_id",
                "B_chr",
                "B_start",
                "B_end",
                "B_description",
            ]) + "\n"
        )

        for agene, bgene in best_hits.items():

            if agene not in table1:
                continue

            if bgene not in table2:
                continue

            b_description = (
                desc2.get(bgene)
                or desc2.get(normalize_id(bgene))
                or "NA"
            )

            out.write(
                "\t".join(
                    table1[agene] +
                    table2[bgene] +
                    [b_description]
                ) + "\n"
            )

            kept += 1

    print(f"Ortholog pairs written: {kept}", file=sys.stderr)


if __name__ == "__main__":
    main()
