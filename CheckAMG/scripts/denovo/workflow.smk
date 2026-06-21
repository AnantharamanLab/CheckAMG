import os
from pathlib import Path

rule all:
    input:
        os.path.join(config["paths"]["output_dir"], "snakemake", "run_inference.done")

# If input nucleotide sequences are provided, run pyrodigal-gv to predict and translate proteins, then compute ESM2 and PST embeddings, and run model inference
if config["query_input"] == "fna":
    rule run_pyrodigal_gv:
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "run_pyrodigal_gv.done"))
        params:
            input_single_contigs = Path(config["paths"]["input_single_contig_fna"]),
            input_bins = config["paths"]["input_bin_fna_dir"],
            wdir = os.path.join(config["paths"]["output_dir"]),
            output_dir = os.path.join(config["paths"]["output_dir"], "pyrodigal-gv"),
            single_contig_prots = os.path.join(config["paths"]["output_dir"], "pyrodigal-gv", "single_contig_proteins.faa"),
            bin_proteins_subdir = directory(os.path.join(config["paths"]["output_dir"], "pyrodigal-gv", "bin_proteins")),
            gene_to_genome = os.path.join(config["paths"]["output_dir"], "gene_to_genome.txt"),
            starting_point = True,
            debug = bool(config["debug"]),
            log = config["log"],
            step = (1,5),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Predicting genes in input contigs with pyrodigal-gv & translating"
        script:
            os.path.join(config["paths"]["common_dir"], "run_pyrodigal.py")

    # Filter translated pyrodigal-gv sequences by AA length and number of CDS per contig for PST
    rule filter_prots:
        input:
            os.path.join(config["paths"]["output_dir"], "snakemake", "run_pyrodigal_gv.done")
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "filter_prots.done"))
        params:
            input_single_contig_prots = Path(os.path.join(config["paths"]["output_dir"], "pyrodigal-gv", "single_contig_proteins.faa")),
            input_bin_prots = Path(os.path.join(config["paths"]["output_dir"], "pyrodigal-gv", "bin_proteins")),
            output_filtered_faa = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa")),
            contig_to_genome = Path(Path(config["paths"]["output_dir"]).joinpath("contig_to_genome.tsv")),
            input_type = "nucl",
            min_cds = 5,
            max_cds = 2048,
            max_protein_length = 20_000,
            debug = bool(config["debug"]),
            log = config["log"],
            step = (2,5),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Filtering translated pyrodigal-gv sequences using a minimum number of CDS of {params.min_cds}"
        script:
            os.path.join(config["paths"]["steps_dir"], "filter_prots.py")

    # Compute ESM2 embeddings and graph from input protein FASTA
    rule compute_query_esm2_graph:
        input:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "filter_prots.done"))
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_esm2.done"))
        params:
            query_fasta_file = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa")),
            model_ckpt = config["paths"]["model_ckpt"],
            esm2_ckpt_dir = config["paths"]["esm2_ckpt_dir"],
            output_dir = config["paths"]["output_dir"],
            esm2_graph_file = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".graphfmt.h5")),
            accelerator = config["accelerator"],
            devices = config["devices"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (3,5),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Compute ESM2 embeddings from input FASTA"
        script:
            os.path.join(config["paths"]["steps_dir"], "fasta_to_esm2.py")

    # Compute PST embeddings from ESM2 graphs
    rule compute_query_pst_embeddings:
        input:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_esm2.done"))
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_pst.done"))
        params:
            esm2_graph_file = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".graphfmt.h5")),
            query_file_pst = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".PST-EMBED.h5")),
            model_ckpt = config["paths"]["model_ckpt"],
            batch_size = config["batch_size"],
            accelerator = config["accelerator"],
            devices = config["devices"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (4,5),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Compute PST embeddings from ESM2"
        script:
            os.path.join(config["paths"]["steps_dir"], "esm2_to_pst_embed.py")

    # Run model inference using PST embeddings (all input modes)
    rule run_inference:
        input:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_pst.done"))
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "run_inference.done"))
        params:
            query_fasta_file = Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa"),
            query_file_pst = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".PST-EMBED.h5")),
            train_data_file = config["paths"]["train_data_file"],
            train_embed_file = config["paths"]["train_embed_file"],
            train_index_file = config["paths"]["train_index_file"],
            train_labels_file = config["paths"]["train_labels_file"],
            model_ckpt = config["paths"]["model_ckpt"],
            accelerator = config["accelerator"],
            devices = config["devices"],
            num_index_cells = config["num_index_cells"],
            num_probe_cells = config["num_probe_cells"],
            knn_k = config["knn"],
            batch_size = config["batch_size"],
            nn_batch_size = config["nn_batch_size"],
            output_file = Path(Path(config["paths"]["output_dir"]).joinpath("predictions.tsv")),
            contig_to_genome_file = Path(Path(config["paths"]["output_dir"]).joinpath("contig_to_genome.tsv")),
            pst_thresholds = config["pst_thresholds"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (5,5),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "PST model inference"
        script:
            os.path.join(config["paths"]["steps_dir"], "pst_embed_to_inference.py")

# If input proteins are provided as FASTA, filter them, compute ESM2 and PST embeddings, and run model inference
elif config["query_input"] == "faa":
    # Filter sequences by AA length and number of CDS per contig for PST
    rule filter_prots:
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "filter_prots.done"))
        params:
            input_single_contig_prots = Path(config["paths"]["input_single_contig_faa"]),
            input_bin_prots = config["paths"]["input_bin_faa"],
            output_filtered_faa = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa")),
            contig_to_genome = Path(Path(config["paths"]["output_dir"]).joinpath("contig_to_genome.tsv")),
            input_type = "prot",
            min_cds = 5,
            max_cds = 2048,
            max_protein_length = 20_000,
            debug = bool(config["debug"]),
            log = config["log"],
            step = (1,4),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Filtering translated pyrodigal-gv sequences using a minimum number of CDS of {params.min_cds}"
        script:
            os.path.join(config["paths"]["steps_dir"], "filter_prots.py")

    # Compute ESM2 embeddings and graph from input protein FASTA
    rule compute_query_esm2_graph:
        input:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "filter_prots.done"))
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_esm2.done"))
        params:
            query_fasta_file = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa")),
            model_ckpt = config["paths"]["model_ckpt"],
            esm2_ckpt_dir = config["paths"]["esm2_ckpt_dir"],
            output_dir = config["paths"]["output_dir"],
            esm2_graph_file = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".graphfmt.h5")),
            accelerator = config["accelerator"],
            devices = config["devices"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (2,4),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Compute ESM2 embeddings from input FASTA"
        script:
            os.path.join(config["paths"]["steps_dir"], "fasta_to_esm2.py")

    # Compute PST embeddings from ESM2 graphs
    rule compute_query_pst_embeddings:
        input:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_esm2.done"))
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_pst.done"))
        params:
            esm2_graph_file = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".graphfmt.h5")),
            query_file_pst = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".PST-EMBED.h5")),
            model_ckpt = config["paths"]["model_ckpt"],
            batch_size = config["batch_size"],
            accelerator = config["accelerator"],
            devices = config["devices"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (3,4),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Compute PST embeddings from ESM2"
        script:
            os.path.join(config["paths"]["steps_dir"], "esm2_to_pst_embed.py")

    # Run model inference using PST embeddings (all input modes)
    rule run_inference:
        input:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_pst.done"))
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "run_inference.done"))
        params:
            query_fasta_file = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa")),
            query_file_pst = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".PST-EMBED.h5")),
            train_data_file = config["paths"]["train_data_file"],
            train_embed_file = config["paths"]["train_embed_file"],
            train_index_file = config["paths"]["train_index_file"],
            train_labels_file = config["paths"]["train_labels_file"],
            model_ckpt = config["paths"]["model_ckpt"],
            accelerator = config["accelerator"],
            devices = config["devices"],
            num_index_cells = config["num_index_cells"],
            num_probe_cells = config["num_probe_cells"],
            knn_k = config["knn"],
            batch_size = config["batch_size"],
            nn_batch_size = config["nn_batch_size"],
            output_file = Path(Path(config["paths"]["output_dir"]).joinpath("predictions.tsv")),
            contig_to_genome_file = Path(Path(config["paths"]["output_dir"]).joinpath("contig_to_genome.tsv")),
            pst_thresholds = config["pst_thresholds"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (4,4)
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "PST model inference"
        script:
            os.path.join(config["paths"]["steps_dir"], "pst_embed_to_inference.py")

# If ESM2 graph formatted file is provided as input, skip directly to PST embedding and inference steps
elif config["query_input"] == "esm2":
    # Compute PST embeddings from ESM2 graphs
    rule compute_query_pst_embeddings:
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_pst.done"))
        params:
            esm2_graph_file = Path(config["paths"]["query_esm2_file"]),
            query_file_pst = Path(Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".PST-EMBED.h5")),
            model_ckpt = config["paths"]["model_ckpt"],
            batch_size = config["batch_size"],
            accelerator = config["accelerator"],
            devices = config["devices"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (1,2),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Compute PST embeddings from ESM2"
        script:
            os.path.join(config["paths"]["steps_dir"], "esm2_to_pst_embed.py")

    # Run model inference using PST embeddings (all input modes)
    rule run_inference:
        input:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "compute_query_pst.done"))
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "run_inference.done"))
        params:
            query_fasta_file = config["paths"]["query_fasta"],
            query_file_pst = Path(config["paths"]["output_dir"]).joinpath("combined_proteins.filtered.faa").with_suffix(".PST-EMBED.h5"),
            train_data_file = config["paths"]["train_data_file"],
            train_embed_file = config["paths"]["train_embed_file"],
            train_index_file = config["paths"]["train_index_file"],
            train_labels_file = config["paths"]["train_labels_file"],
            model_ckpt = config["paths"]["model_ckpt"],
            accelerator = config["accelerator"],
            devices = config["devices"],
            num_index_cells = config["num_index_cells"],
            num_probe_cells = config["num_probe_cells"],
            knn_k = config["knn"],
            batch_size = config["batch_size"],
            nn_batch_size = config["nn_batch_size"],
            output_file = Path(Path(config["paths"]["output_dir"]).joinpath("predictions.tsv")),
            contig_to_genome_file = Path(Path(config["paths"]["output_dir"]).joinpath("contig_to_genome.tsv")),
            pst_thresholds = config["pst_thresholds"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (2,2),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "PST model inference"
        script:
            os.path.join(config["paths"]["steps_dir"], "pst_embed_to_inference.py")

# If PST embeddings are provided as input, skip directly to inference step
elif config["query_input"] == "pst_embed":
    rule run_inference:
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "run_inference.done"))
        params:
            query_fasta_file = config["paths"]["query_fasta"],
            query_file_pst = config["paths"]["query_file_pst"],
            train_data_file = config["paths"]["train_data_file"],
            train_embed_file = config["paths"]["train_embed_file"],
            train_index_file = config["paths"]["train_index_file"],
            train_labels_file = config["paths"]["train_labels_file"],
            model_ckpt = config["paths"]["model_ckpt"],
            accelerator = config["accelerator"],
            devices = config["devices"],
            num_index_cells = config["num_index_cells"],
            num_probe_cells = config["num_probe_cells"],
            knn_k = config["knn"],
            batch_size = config["batch_size"],
            nn_batch_size = config["nn_batch_size"],
            output_file = Path(Path(config["paths"]["output_dir"]).joinpath("predictions.tsv")),
            contig_to_genome_file = Path(Path(config["paths"]["output_dir"]).joinpath("contig_to_genome.tsv")),
            pst_thresholds = config["pst_thresholds"],
            debug = bool(config["debug"]),
            log = config["log"],
            step = (1,1),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "PST model inference"
        script:
            os.path.join(config["paths"]["steps_dir"], "pst_embed_to_inference.py")

else:
    raise ValueError(f"Unsupported query input mode: {config['query_input']}")
