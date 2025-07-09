#!/usr/bin/env python3

import argparse
import textwrap
import sys
import psutil
from CheckAMG.scripts import CheckAMG_annotate, download_dbs
from CheckAMG.scripts.checkAMG_ASCII import ASCII
from importlib.metadata import version

__version__ = version("checkamg")

available_memory_gb = psutil.virtual_memory().available / (1024 ** 3) # Get available memory in GB

class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
     def _fill_text(self, text, width, indent):
          text = textwrap.dedent(text).strip()
          return textwrap.fill(text, width, initial_indent=indent, subsequent_indent=indent)

def main():
     parser = argparse.ArgumentParser(description="CheckAMG: automated identification and curation of Auxiliary Metabolic Genes (AMGs),"
                                                  " Auxiliary Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs)"
                                                  " in viral genomes.", formatter_class=CustomHelpFormatter)
     parser.add_argument("-v", "--version", action="version", version=f"CheckAMG {__version__}")
     subparsers = parser.add_subparsers(help="CheckAMG modules", dest="command")
          
     download_parser = subparsers.add_parser(
          "download",
          help="Download the databases required by CheckAMG.",
          description="Download the databases required by CheckAMG. This requires ~40 GB of disk space (or ~21 GB finally, if the '--remove' argument is provided).",
          formatter_class=CustomHelpFormatter)
     download_parser.add_argument(
          "-d", "--db_dir", type=str, required=True,
          help="Path to the directory where CheckAMG databases will be downloaded (Required).")
     download_parser.add_argument(
          "-f", "--force", action="store_true", default=False,
          help="Force re-download of databases even if they already exist (default: %(default)s).")
     download_parser.add_argument(
          "-r", "--rm_hmm", action="store_true", default=False,
          help="Remove human-readable HMM files from the database directory after downloading, to save space. CheckAMG only needs the binary files (default: %(default)s).")
     download_parser.add_argument(
          "-t", "--threads", type=int, default=10,
          help="Number of threads to use for downloading and compiling databases (default: %(default)s).")
     
     annotate_parser = subparsers.add_parser("annotate",
                                        help="Predict and curate auxiliary genes in viral genomes based on functional annotations and genomic context.",
                                        description="Predict and curate auxiliary genes in viral genomes based on functional annotations and genomic context.",
                                        formatter_class=CustomHelpFormatter)

     annotate_required = annotate_parser.add_argument_group('required arguments')
     annotate_required.add_argument("-d", "--db_dir", type=str, required=True,
                              help="Path to CheckAMG database files (Required).")
     annotate_required.add_argument("-o", "--output", type=str, required=True,
                              help="Output directory for all generated files and folders (Required).")
     annotate_required.add_argument("-g", "--genomes", type=str, required=False,
                              help="Input viral genome(s) in nucleotide fasta format (.fna or .fasta). Expectation is that "
                                   "individual virus genomes are single contigs.")
     annotate_required.add_argument("-vg", "--vmags", type=str, required=False,
                              help="Path to folder containing vMAGs (multiple contigs) rather than single-contig viral genomes. "
                                   "Expectation is that the folder contains one .fna or .fasta file "
                                   "per virus genome and that each genome contains multiple contigs.")
     annotate_required.add_argument("-p", "--proteins", type=str, required=False,
                              help="Input viral genome(s) in amino-acid fasta format (.faa or .fasta). Required if --input_type is prot. "
                                   "Expectations are that the amino-acid sequence headers are in Prodigal format (>[CONTIG NAME]_[CDS NUMBER] # START # END # FRAME # ...) "
                                   "and that each contig encoding proteins represents a single virus genome.")
     annotate_required.add_argument("-vp", "--vmag_proteins", type=str, required=False,
                              help="Path to folder containing vMAGs (multiple contigs) in amino-acid fasta format (.faa or .fasta) "
                                   "rather than single-contig viral genomes. Expectation is that the folder contains one .faa or .fasta file "
                                   "per virus genome and that each genome file contains amino-acid sequences encoded on multiple contigs. "
                                   "Required if --input_type is 'prot'.")
          
     annotate_parser.add_argument("--input_type", type=str, required=False, default="nucl",
                              help="Specifies whether the input files are nucleotide genomes (nucl) or translated amino-acid genomes (prot). "
                                   "Providing proteins as input will skip the pyrodigal-gv step, "
                                   "but it will be unable to tell whether viral genomes are circular, potentially losing additional evidence "
                                   "for verifying the viral origin of putative auxiliary genes. (default: %(default)s).")
     annotate_parser.add_argument("-l", "--min_len", type=int, required=False, default = 5000,
                              help="Minimum length in base pairs for input sequences (default: %(default)s).")
     annotate_parser.add_argument("-f", "--min_orf", type=int, required=False, default = 4,
                              help="Minimum number of open reading frames (proteins) inferred by pyrodigal-gv for input sequences (default: %(default)s).")
     annotate_parser.add_argument("-n", "--min_annot", type=float, required=False, default=0.20,
                              help="Minimum percentage (0.0-1.0) of genes in a genome/contig required to have been assigned a "
                                   "functional annotation using the CheckAMG database to be considered for contextual analysis. "
                                   "(default: %(default)s).")
     annotate_parser.add_argument("-c", "--cov_fraction", type=float, required=False, default=0.5,
                              help="Minimum covered fraction (of the user viral protein) for HMM searches (default: %(default)s).")
     annotate_parser.add_argument("-e", "--evalue", type=float, required=False, default=1e-5,
                              help="Maximum fallback E-value for HMM searches when database-provided cutoffs are not available (default: %(default)s).")
     annotate_parser.add_argument("-b", "--bit_score", type=int, required=False, default=50,
                              help="Minimum fallback bit score for HMM searches when database-provided cutoffs are not available (default: %(default)s).")
     annotate_parser.add_argument("-bh", "--bitscore_fraction_heuristic", type=float, required=False, default=0.5,
                              help="Retain HMM hits scoring at least this fraction of the database-provided threshold under heuristic filtering (default: %(default)s).")
     annotate_parser.add_argument("-s", "--scaling_factor", type=float, required=False, default=1.5, min=1.0,
                              help="Scaling factor used to multiply the minimum bit score and minimum covered fraction provided by the '-b' and '-c' arguments to come up with"
                                   " a heuristic, stricter threshold for HMM hits (this is ONLY used when curating gene annotations that match to 'soft' filter keywords; default: %(default)s)."
                                   " This is to avoid making functional predictions to viral genes that are often misannotated with auxiliary functions, requiring stronger "
                                   " HMM search results to be assigned a function that is often incorrect if less-strict thresholds were used. If database-provided, trusted bit score cutoffs are available "
                                   " for the matching HMMs, those are used instead of the heuristic threshold for minimum bit score, but the scaling factor is still applied to the "
                                   " minimum covered fraction. This is not used for gene annotations that match to 'hard' filter keywords, which are always excluded from the final auxiliary gene predictions.")
     annotate_parser.add_argument("-Z", "--window_size", type=int, required=False, default=25000,
                              help="Size in base pairs of the window used to calculate the average VL-score of genes on a contig (default: %(default)s).")
     annotate_parser.add_argument("-F", "--max_flank", type=int, required=False, default=5000,
                              help="Maximum length in base pairs to check on the left/right flanks of potentially auxiliary genes when checking for virus-like genes and non-virus-like genes (default: %(default)s).")
     annotate_parser.add_argument("-V", "--min_flank_Vscore", type=float, required=False, default=10.0,
                              help="Minimum V-score of genes in flanking regions required to verify a potential auxiliary gene as viral and not host sequence contamination (0.0-10.0) (default: %(default)s).")
     annotate_parser.add_argument("-H", "--use_hallmark", required=False, default=False, action=argparse.BooleanOptionalAction,
                              help="Use viral hallmark gene annotations instead of V-scores when checking flanking regions of potential auxiliary genes for viral verification (default: %(default)s).")
     annotate_parser.add_argument("-t", "--threads", type=int, required=False, default=10,
                              help="Number of threads to use for pyrodigal-gv and pyhmmer (default: %(default)s).")
     annotate_parser.add_argument("-m", "--mem", type=int, required=False, default=round(available_memory_gb*0.80), # 80% of available memory
                              help="Maximum amount of memory allowed to be allocated in GB (default: 80%% of available).")
     annotate_parser.add_argument("--debug", required=False, default=False, action=argparse.BooleanOptionalAction,
                              help="Log CheckAMG with debug-level detail (default: %(default)s).")

     de_novo_parser = subparsers.add_parser(
          "de-novo",
          help="(Not yet implemented) Predict auxiliary genes with an annotation-independent method using Protein Set Transformer.",
          description="Not yet implemented.",
          formatter_class=CustomHelpFormatter)

     aggregate_parser = subparsers.add_parser(
          "aggregate",
          help="(Not yet implemented) Aggregate the results of the CheckAMG annotate and de-novo modules to produce a final report of auxiliary gene predictions.",
          description="Not yet implemented.",
          formatter_class=CustomHelpFormatter)
          
     end_to_end_parser = subparsers.add_parser(
          "end-to-end",
          help="(Not yet implemented) Executes CheckAMG annotate, de-novo, and aggregate in tandem.",
          description="Not yet implemented.",
          formatter_class=CustomHelpFormatter)
     
     if "--version" not in sys.argv and "-v" not in sys.argv:
          print(ASCII)
          sys.stdout.flush()
          
     # Validate that only one subcommand is given
     subcommands = {"download", "annotate", "de-novo", "aggregate", "end-to-end"}
     used_subcommands = [arg for arg in sys.argv[1:] if arg in subcommands]
     if len(used_subcommands) > 1:
          parser.error(f"Too many arguments provided ({', '.join(used_subcommands)}). Please specify only one CheckAMG module to run.")
          
     args = parser.parse_args()
               
     if args.command == "download":
          download_dbs.download_all(dest=args.db_dir, force=args.force, threads=args.threads)
          if args.rm_hmm:
               download_dbs.remove_human_readable_files(dest=args.db_dir)
     elif args.command == "annotate":
          if args.input_type == "nucl" and not args.genomes and not args.vmags:
               parser.error("At least one of --genomes or --vmags is required when --input_type is 'nucl'.")
          if args.input_type == "nucl" and (args.proteins or args.vmag_proteins):
               parser.error("Cannot provide --proteins or --vmag_proteins when --input_type is 'nucl'.")
          if args.input_type == "prot" and not args.proteins and not args.vmag_proteins:
               parser.error("At least one of --proteins or --vmag_proteins is required when --input_type is 'prot'.")
          if args.input_type == "prot" and (args.genomes or args.vmags):
               parser.error("Cannot provide --genomes or --vmags when --input_type is 'prot'.")
          if (args.genomes and args.proteins) or (args.vmags and args.vmag_proteins):
               parser.error("Cannot provide both --genomes/--vmags and --proteins/--vmag_proteins.")
          CheckAMG_annotate.create_output_dir(args.output)
          config_path = CheckAMG_annotate.generate_config(args)
          CheckAMG_annotate.run_snakemake(config_path, args)
     elif args.command == "de-novo":
          print("CheckAMG de novo functionality will be implemented here.")
     elif args.command == "aggregate":
          print("CheckAMG aggregate functionality will be implemented here.")
     elif args.command == "end-to-end":
          print("CheckAMG end-to-end pipeline functionality will be implemented here.")
     else:
          print("Error: Please specify a CheckAMG module to run.", file=sys.stderr)
          parser.print_help()
          sys.exit(1)

if __name__ == "__main__":
     main()

