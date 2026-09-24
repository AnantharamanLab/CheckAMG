#!/usr/bin/env python3

import os
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor
import pyfastx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

DELIM = "~"
VALID_MGE_TYPES = {"Phage", "IS_Tn", "Conjugative_element"}

def build_id(contig, contig_type, start=None, end=None):
    if not start and not end:
        return f"{contig}{DELIM}{contig_type}"
    elif start is not None and end is not None:
        return f"{contig}{DELIM}{start}-{end}{DELIM}{contig_type}"
    else:
        raise ValueError(f"Both start and end must be provided together for contig {contig} with type {contig_type}")

def load_fasta_headers(header_file):
    headers = set()
    with open(header_file) as f:
        for line in f:
            if line.startswith('>'):
                headers.add(line[1:].split()[0])
    return headers

def build_needed_contigs(*dfs):
    contigs = set()
    for df in dfs:
        if hasattr(df, 'iter_rows'):
            for row in df.iter_rows(named=True):
                contigs.add(row['contig'])
        else:
            for row in df:
                contigs.add(row['contig'])
    return contigs

def build_seq_dict_from_pyfastx(fasta_path, needed_contigs):
    fx = pyfastx.Fasta(fasta_path, build_index=True)
    seq_dict = {}
    for contig in needed_contigs:
        if contig in fx:
            seq_dict[contig] = str(fx[contig][:])
    return seq_dict

def partition_two_bins(mge_df):
    allowed = mge_df.filter(pl.col('mge_category').is_in(VALID_MGE_TYPES))
    mge_contigs = allowed['contig'].unique().to_list()
    mge_contigs = sorted(mge_contigs)
    half = len(mge_contigs) // 2
    mge_contigs_1 = set(mge_contigs[:half])
    mge_contigs_2 = set(mge_contigs[half:])
    mge1 = allowed.filter(pl.col('contig').is_in(mge_contigs_1))
    mge2 = allowed.filter(pl.col('contig').is_in(mge_contigs_2))
    return mge1, mge2

def extract_mge_only(mge_rows, seq_dict, threads):
    records, regions = [], []
    def worker(row):
        contig, start, end, cat = row['contig'], row['start'], row['end'], row['mge_category']
        if cat == "Phage":
            contig_type = "standalone_virus"
            region_type = "Viral"
            region_subtype = "Provirus"
        else:
            contig_type = "standalone_mge"
            region_type = "MGE"
            region_subtype = cat
        rid = build_id(contig, contig_type, start, end)
        seq = seq_dict[contig][start:end]
        rec = SeqRecord(Seq(seq), id=rid, description="")
        # start always 0 for standalone, region is whole thing
        return rec, (rid, 0, len(seq), contig_type, region_type, region_subtype)
    with ThreadPoolExecutor(max_workers=threads) as exe:
        results = exe.map(worker, mge_rows)
        for rec, annot in results:
            records.append(rec)
            regions.append(annot)
    return records, regions

def extract_chrom_mge_only(mge_df, seq_dict, threads):
    records, regions = [], []
    mge_map = {}
    for row in mge_df.iter_rows(named=True):
        mge_map.setdefault(row['contig'], []).append(row)
    def worker(contig):
        seq = seq_dict[contig]
        start, end = 1, len(seq)
        rid = build_id(contig, "chromosome_mixed", start, end)
        rec = SeqRecord(Seq(seq), id=rid, description="")
        regs, curr = [], 0
        rows = sorted(mge_map[contig], key=lambda r: r['start'])
        for row in rows:
            s, e, cat = row['start'], row['end'], row['mge_category']
            # Host region before MGE
            if curr < s:
                regs.append((rid, curr, s, "chromosome_mixed", "Host", "Host"))
            # MGE region
            if cat == "Phage":
                regs.append((rid, s, e, "chromosome_mixed", "Viral", "Provirus"))
            else:
                regs.append((rid, s, e, "chromosome_mixed", "MGE", cat))
            curr = e
        if curr < len(seq):
            regs.append((rid, curr, len(seq), "chromosome_mixed", "Host", "Host"))
        return rec, regs
    with ThreadPoolExecutor(max_workers=threads) as exe:
        results = exe.map(worker, mge_map.keys())
        for rec, regs in results:
            records.append(rec)
            regions.extend(regs)
    return records, regions

def extract_chromosome_only(contigs, seq_dict, mge_df, threads):
    mge_regions = {}
    for row in mge_df.iter_rows(named=True):
        mge_regions.setdefault(row['contig'], []).append((row['start'], row['end']))
    records, regions = [], []
    def worker(contig):
        seq = seq_dict[contig]
        cut_regions = sorted(mge_regions.get(contig, []), key=lambda x: x[0])
        kept_segs = []
        region_coords = []
        prev = 0
        for s, e in cut_regions:
            if prev < s:
                kept_segs.append((prev, s))
            prev = e
        if prev < len(seq):
            kept_segs.append((prev, len(seq)))
        pieces = []
        for start, end in kept_segs:
            seg = seq[start:end]
            if seg:
                pieces.append(seg)
                region_coords.append((start, end))
        stitched_seq = ''.join(pieces)
        start, end = 1, len(stitched_seq)
        rid = build_id(contig, "chromosome_only", start, end)
        rec = SeqRecord(Seq(stitched_seq), id=rid, description="")
        region_list = []
        offset = 0
        for (start, end) in region_coords:
            seg_len = end - start
            region_list.append((rid, offset, offset+seg_len, "chromosome_only", "Host", "Host"))
            offset += seg_len
        return rec, region_list
    with ThreadPoolExecutor(max_workers=threads) as exe:
        results = exe.map(worker, contigs)
        for rec, regs in results:
            records.append(rec)
            regions.extend(regs)
    return records, regions

def write_outputs(bin_name, records, regions, outdir):
    fa = os.path.join(outdir, f"{bin_name}.fasta")
    ts = os.path.join(outdir, f"{bin_name}_regions.tsv")
    seen_seqs = {}
    unique_records = []
    for rec in records:
        seq = str(rec.seq)
        if rec.id in seen_seqs:
            if seen_seqs[rec.id] != seq:
                raise ValueError(
                    f"Duplicate header '{rec.id}' in bin '{bin_name}' maps to different sequences"
                )
            continue
        seen_seqs[rec.id] = seq
        unique_records.append(rec)
    seen_regions = set()
    unique_regions = []
    for tup in regions:
        if tup in seen_regions:
            continue
        seen_regions.add(tup)
        unique_regions.append(tup)
    SeqIO.write(unique_records, str(fa), 'fasta')
    with open(ts, 'w') as f:
        f.write('contig\tstart\tend\tcontig_type\tregion_type\tregion_subtype\n')
        for tup in unique_regions:
            f.write("\t".join(map(str, tup)) + "\n")
    logging.info(f"Wrote {bin_name}: {fa} / {ts} ({len(unique_records):,} records)")

def main(args):
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(levelname)s:%(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)
    os.makedirs(args.outdir, exist_ok=True)

    logger.info(f"Reading MGE table: {args.mge_table}")
    mge = pl.read_csv(args.mge_table, separator='\t')
    logger.info(f"Initial MGE table: {mge.height:,} rows")
    logger.info(f"Filtering MGEs by minimum length: {args.min_len:,} bp")
    mge = mge.filter(pl.col('mge_length') >= args.min_len)
    logger.info(f"After length filter: {mge.height:,} rows")
    mge = mge.with_columns([
        pl.col('mge_genome_position').str.split_exact(':',1).struct.field('field_0').alias('contig'),
        pl.col('mge_genome_position').str.split_exact(':',1).struct.field('field_1').alias('pos')
    ]).with_columns([
        pl.col('pos').str.split_exact('-',1).struct.field('field_0').cast(pl.Int64).alias('start'),
        pl.col('pos').str.split_exact('-',1).struct.field('field_1').cast(pl.Int64).alias('end')
    ]).drop(['mge_genome_position','pos'])
    logger.info(f"Filtering by headers from: {args.headers}")
    hdrs = load_fasta_headers(args.headers)
    logger.info(f"{len(hdrs):,} headers loaded.")
    mge = mge.filter(pl.col('contig').is_in(hdrs))
    logger.info(f"After header filter: {mge.shape[0]:,} rows")
    mge = mge.filter(pl.col('mge_category').is_in(VALID_MGE_TYPES))
    logger.info(f"After MGE type filter: {mge.shape[0]:,} rows")

    logger.info(f"Strictly partitioning MGEs into two non-overlapping bins of contigs...")
    mge1, mge2 = partition_two_bins(mge)
    logger.info(f"Standalone MGE bin contigs: {mge1['contig'].unique().to_list().__len__():,}, chrom+MGE bin contigs: {mge2['contig'].unique().to_list().__len__():,}")

    needed_contigs = build_needed_contigs(mge1, mge2)
    logger.info(f"Reading {len(needed_contigs):,} required contigs from FASTA...")
    seq_dict = build_seq_dict_from_pyfastx(args.fasta, needed_contigs)

    # Standalone bins (split for chromosome_only, write only after)
    rows1 = [row for row in mge1.iter_rows(named=True)]
    recs1, regs1 = extract_mge_only(rows1, seq_dict, args.threads)

    # Split mge_records by contig to ensure no overlap in chromosome_only and standalone_mge
    provirus_records = []
    provirus_regions = []
    mge_records = []
    mge_regions = []

    for rec, (rid, s, e, contig_type, region_type, region_subtype) in zip(recs1, regs1):
        if contig_type == "standalone_virus":
            provirus_records.append(rec)
            provirus_regions.append((rid, s, e, contig_type, region_type, region_subtype))
        else:
            mge_records.append(rec)
            mge_regions.append((rid, s, e, contig_type, region_type, region_subtype))

    # chromosome_only/standalone_mge split by contig
    if len(mge_records) > 2:
        contigs = [rec.id.split(DELIM)[0] for rec in mge_records]
        unique_contigs = sorted(set(contigs))
        half = len(unique_contigs) // 2
        chr_only_contigs = set(unique_contigs[half:])
        mge_only_contigs = set(unique_contigs[:half])

        # Chromosome only: excise all MGEs from these contigs
        chr_seq_dict = build_seq_dict_from_pyfastx(args.fasta, chr_only_contigs)
        chr_only_records, chr_only_regions = extract_chromosome_only(chr_only_contigs, chr_seq_dict, mge1, args.threads)
        write_outputs('chromosome_only', chr_only_records, chr_only_regions, args.outdir)
        logger.info(f"Wrote chromosome_only bin: {os.path.join(args.outdir, 'chromosome_only.fasta')} / {os.path.join(args.outdir, 'chromosome_only_regions.tsv')}")

        # Standalone MGE: only write MGE records for mge_only_contigs
        mge_final_records = [rec for rec in mge_records if rec.id.split(DELIM)[0] in mge_only_contigs]
        mge_final_regions = [r for i, r in enumerate(mge_regions) if mge_records[i].id.split(DELIM)[0] in mge_only_contigs]
        write_outputs('standalone_mge', mge_final_records, mge_final_regions, args.outdir)
        logger.info(f"Wrote standalone_mge bin: {os.path.join(args.outdir, 'standalone_mge.fasta')} / {os.path.join(args.outdir, 'standalone_mge_regions.tsv')}")
    else:
        write_outputs('standalone_mge', mge_records, mge_regions, args.outdir)
        logger.info(f"Wrote standalone_mge bin: {os.path.join(args.outdir, 'standalone_mge.fasta')} / {os.path.join(args.outdir, 'standalone_mge_regions.tsv')}")

    # Write standalone provirus last
    write_outputs('standalone_virus', provirus_records, provirus_regions, args.outdir)
    logger.info(f"Wrote standalone_virus bin: {os.path.join(args.outdir, 'standalone_virus.fasta')} / {os.path.join(args.outdir, 'standalone_provirus_regions.tsv')}")

    # Chromosome with any MGEs (phage, non-phage, or both)
    recs2, regs2 = extract_chrom_mge_only(mge2, seq_dict, args.threads)
    write_outputs('chromosome_mixed', recs2, regs2, args.outdir)
    logger.info(f"Wrote chromosome_mixed bin: {os.path.join(args.outdir, 'chromosome_mixed.fasta')} / {os.path.join(args.outdir, 'chromosome_mixed_regions.tsv')}")

    logger.info("Done.")

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mge_table', required=True)
    parser.add_argument('--headers', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--min_len', required=False, default=1000, type=int,
                        help="Minimum length of MGE to consider (default: 1000)")
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--threads', type=int, default=50)
    args = parser.parse_args()
    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    import polars as pl
    main(args)
