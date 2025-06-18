import os

rule all:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "make_final_table.done")

# Filter input sequences by length
rule filter_by_length:
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "filter_by_length.done"))
    params:
        build_or_annotate = "build",
        input_single_contig_genomes = config["input_single_contig_genomes"],
        input_vmag_fastas = config["input_vmag_fastas"],
        min_len = config["min_len"],
        output = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_fna_by_length"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Filtering input sequences using a minimum length of {params.min_len} bp"
    script:
        os.path.join(config["paths"]["scripts_dir"], "filter_by_length.py")

# Check circulatiry of user genomes
rule check_circular:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "filter_by_length.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "check_circular.done"))
    params:
        build_or_annotate = "build",
        input_single_contig_genomes = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_fna_by_length", "single_contig_genomes.fna"),
        input_vmag_fastas = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_fna_by_length", "vMAG_fna"),
        tr_min_len = 20,
        tr_max_len = 1000,
        tr_max_count = 8,
        tr_max_ambig = 0.2,
        tr_max_basefreq = 0.70,
        kmer_max_freq = 1.5,
        k = 15,
        circularity_tbl = os.path.join(config["paths"]["output_dir"], "wdir", "circular_contigs.tsv"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Checking the circularity of input sequences and writing results to {params.circularity_tbl}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "check_circular.py")

# Annotate user genomes
rule run_pyrodigal_gv:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "check_circular.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "run_pyrodigal_gv.done"))
    params:
        build_or_annotate = "build",
        input_single_contig_genomes = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_fna_by_length", "single_contig_genomes.fna"),
        input_vmag_fastas = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_fna_by_length", "vMAG_fna"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir"),
        output_dir = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv"),
        single_contig_prots = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv", "single_contig_proteins.faa"),
        vmag_proteins_subdir = directory(os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv", "vMAG_proteins")),
        gene_to_genome = os.path.join(config["paths"]["output_dir"], "wdir", "gene_to_genome.txt"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Predicting genes in input genomes with pyrodigal-gv & translating"
    script:
        os.path.join(config["paths"]["scripts_dir"], "run_pyrodigal.py")

# Filter translated pyrodigal-gv sequences by minimum number of CDS
rule filter_by_cds:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "run_pyrodigal_gv.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "filter_by_cds.done"))
    params:
        build_or_annotate = "build",
        input_type = "nucl",
        input_prot_subdir = os.path.join(config["paths"]["output_dir"], "wdir", "pyrodigal-gv"),
        min_cds = config["min_cds"],
        output = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Filtering translated pyrodigal-gv sequences using a minimum number of CDS of {params.min_cds}"
    script:
        os.path.join(config["paths"]["scripts_dir"], "filter_by_cds.py")

# Assign annotations to the proteins in the database
rule assign_annots:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "filter_by_cds.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "annotate_hmm.done"))
    params:
        build_or_annotate = "build",
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        hmm_vscores = os.path.join(config["paths"]["files_dir"], "vscores.csv"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir"),
        cov_fraction = config["mmseqs_params"]["cov_fraction"],
        min_bitscore = 50,
        max_evalue = 1e-5,
        vscores = os.path.join(config["paths"]["output_dir"], "wdir", "vscores.tsv"),
        all_hmm_results = os.path.join(config["paths"]["output_dir"], "wdir", "hmm_results.tsv"),
        db_dir = config["paths"]["db_dir"],
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Assigning functions input proteins using an HMMsearch with the databases in {params.db_dir}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "annotate_hmm.py")

# Obtain gene information from input (prodigal-formatted) .faa and genome information from...
rule index_genes:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "annotate_hmm.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "index_genes.done"))
    params:
        build_or_annotate = "build",
        cluster_taxa_levels = config["cluster_ranks"],
        gene_index = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index.tsv"),
        out_parent = os.path.join(config["paths"]["output_dir"], "wdir"),
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        vmag_proteins_subdir = directory(os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds", "vMAG_proteins")),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Writing gene- and genome-level data to {params.gene_index}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "map_protein_data.py")

# Merge the annotations with the protein database
rule add_annots:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "index_genes.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "add_annots.done"))
    params:
        build_or_annotate = "build",
        gene_index = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index.tsv"),
        gene_index_annotated = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        vscores = os.path.join(config["paths"]["output_dir"], "wdir", "vscores.tsv"),
        all_hmm_results = os.path.join(config["paths"]["output_dir"], "wdir", "hmm_results.tsv"),
        db_dir = config["paths"]["db_dir"],
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Adding annotations to {params.gene_index} and writing to {params.gene_index_annotated}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "add_annots.py")

# Analyze the genomic context of annotations
rule genome_context:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "add_annots.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "genome_context.done"))
    params:
        build_or_annotate = "build",
        outparent = os.path.join(config["paths"]["output_dir"], "results"),
        context_table = os.path.join(config["paths"]["output_dir"], "results", "genes_genomic_context.tsv"),
        gene_index_annotated = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        circular_contigs = os.path.join(config["paths"]["output_dir"], "wdir", "circular_contigs.tsv"),
        annotation_percent_threshold = config["annotation_percent_threshold"],
        min_window_avg_lscores = config["min_window_avg_lscores"],
        window_size = config["window_size"],
        minimum_flank_vscore = config["minimum_flank_vscore"],
        max_flank_length = config["max_flank_length"],
        use_hallmark = config["use_hallmark"],
        hallmark_path = os.path.join(config["paths"]["files_dir"], "viral_hallmark_genes.csv"),
        mobile_genes_path = os.path.join(config["paths"]["files_dir"], "mobile_genes.csv"),
        vscore_ref = os.path.join(config["paths"]["files_dir"], "vscores.csv"),
        tmp_dir = os.path.join(config["paths"]["output_dir"], "wdir"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Analyzing the genomic context of V-scores and L-scores in {params.gene_index_annotated} and writing results to {params.context_table}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "genome_context.py")

# Curate the predicted functions based on their genomic context
rule curate_annots:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "genome_context.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "curate_results.done"))
    params:
        build_or_annotate = "build",
        context_table = os.path.join(config["paths"]["output_dir"], "results", "genes_genomic_context.tsv"),
        metabolism_table = os.path.join(config["paths"]["files_dir"], "AMGs.tsv"),
        physiology_table = os.path.join(config["paths"]["files_dir"], "APGs.tsv"),
        regulation_table = os.path.join(config["paths"]["files_dir"], "AReGs.tsv"),
        metabolism_table_out = os.path.join(config["paths"]["output_dir"], "results", "metabolic_genes_curated.tsv"),
        physiology_table_out = os.path.join(config["paths"]["output_dir"], "results", "physiology_genes_curated.tsv"),
        regulation_table_out = os.path.join(config["paths"]["output_dir"], "results", "regulation_genes_curated.tsv"),
        bypass_min_bitscore = 80,
        bypass_min_cov = 0.8,
        all_annot_out_table = os.path.join(config["paths"]["output_dir"], "results", "gene_annotations.tsv"),
        hmm_ref = os.path.join(config["paths"]["files_dir"], "hmm_id_to_name.csv"),
        false_amgs = os.path.join(config["paths"]["files_dir"], "false_amgs.txt"),
        false_apgs = os.path.join(config["paths"]["files_dir"], "false_apgs.txt"),
        false_aregs = os.path.join(config["paths"]["files_dir"], "false_aregs.txt"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Writing the curated metabolic/regulatory/physiology protein results."
    script:
        os.path.join(config["paths"]["scripts_dir"], "curate_annots.py")

# Create MMseqs2 database from input protein files and run mmseqs search and mmseqs convertalis
rule mmseqs_filtered_proteins:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "curate_results.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "mmseqs_filtered_proteins.done"))
    params:
        protein_db = os.path.join(config["paths"]["output_dir"], "wdir", "protein_db", "protein_db"),
        input_prot_subdir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        out_db_dir = os.path.join(config["paths"]["output_dir"], "wdir", "protein_db"),
        wdir_parent = os.path.join(config["paths"]["output_dir"], "wdir"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "search"),
        search_parent = os.path.join(config["paths"]["output_dir"], "wdir", "protein_search_db"),
        search_db = os.path.join(config["paths"]["output_dir"], "wdir", "protein_search_db", "protein_search_db"),
        search_sensitivity = config["mmseqs_params"]["sensitivity"],
        min_seq_id = 0, # Keep minimum sequence identity at 0.0 since we want to sample the entire range of sequence identity
        cov_fraction = config["mmseqs_params"]["cov_fraction"],
        search_output = os.path.join(config["paths"]["output_dir"], "wdir", "search", "protein_search_output.tsv"),
        debug = bool(config["debug"]),
        log = config["log"]
    resources:
        mem = config["mem_limit"]
    threads:
        config["threads"]
    message:
        "Creating MMseqs2 database from input protein files and running mmseqs search"
    script:
        os.path.join(config["paths"]["scripts_dir"], "cluster_filtered_proteins.py")

# Preprocess search output for AAI calculations
rule preprocess_aai:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "mmseqs_filtered_proteins.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "preprocess_aai.done"))
    threads:
        config["threads"]
    params:
        search_output = os.path.join(config["paths"]["output_dir"], "wdir", "search", "protein_search_output.tsv"),
        input_prot_subdir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        gene_map_file = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "gene_map.tsv"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "aai"),
        debug = bool(config["debug"]),
        log = config["log"]
    message:
        "Preparing gene-to-genome map for AAI calculations."
    resources:
        mem = config["mem_limit"]
    script:
        os.path.join(config["paths"]["scripts_dir"], "preprocess_genome_aai.py")

# Compute the amino-acid identity of all proteins in the database
rule compute_aai:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "preprocess_aai.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "compute_aai.done"))
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    params:
        search_output = os.path.join(config["paths"]["output_dir"], "wdir", "search", "protein_search_output.processed.tsv"),
        gene_map_file = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "gene_map.tsv"),
        cov_fraction = config["mmseqs_params"]["cov_fraction"],
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "aai"),
        debug = bool(config["debug"]),
        log = config["log"]
    message:
        "Computing pairwise amino-acid identities of all proteins in the database."
    script:
        os.path.join(config["paths"]["scripts_dir"], "compute_genome_aai.py")

# Cluster genomes using amino-acid identities
rule cluster_aai:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "compute_aai.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "cluster_aai.done"))
    threads:
        config["threads"]
    params:
        gene_map = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "gene_map.tsv"),
        cluster_taxa_levels = config["cluster_ranks"],
        duckdb = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "alignments_and_clusters.duckdb"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "aai"),
        debug = bool(config["debug"]),
        log = config["log"]
    resources:
        mem = config["mem_limit"]
    message:
        "Clustering genomes based on AAI results at the taxonomic rank(s): {params.cluster_taxa_levels}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "cluster_genome_aai.py")

# Protein clustering for each genome cluster
rule cluster_proteins_within_genome_cluster:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "cluster_aai.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "protein_clustering.done"))
    threads:
        config["threads"]
    params:
        duckdb = os.path.join(config["paths"]["output_dir"], "wdir", "aai", "alignments_and_clusters.duckdb"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir"),
        protein_clusters_tsv = os.path.join(config["paths"]["output_dir"], "wdir", "genome_clusters_protein_clusters.tsv"),
        cluster_taxa_levels = config["cluster_ranks"],
        min_seq_id = config["mmseqs_params"]["min_seq_id"],
        cov_fraction = config["mmseqs_params"]["cov_fraction"],
        debug = bool(config["debug"]),
        log = config["log"]
    resources:
        mem = config["mem_limit"]
    message:
        "Clustering proteins within each genome cluster based on amino-acid sequence identity."
    script:
        os.path.join(config["paths"]["scripts_dir"], "cluster_proteins_within_genome_clusters.py")

# Add protein and genome cluster information to the gene index
rule add_cluster_info_to_gene_index:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "protein_clustering.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "add_cluster_info_to_gene_index.done"))
    threads:
        config["threads"]
    params:
        gene_index_annotated = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        protein_clusters_tsv = os.path.join(config["paths"]["output_dir"], "wdir", "genome_clusters_protein_clusters.tsv"),
        out_parent = os.path.join(config["paths"]["output_dir"], "wdir"),
        debug = bool(config["debug"]),
        log = config["log"]
    resources:
        mem = config["mem_limit"]
    message:
        "Adding genome and protein cluster information to the gene index."
    script:
        os.path.join(config["paths"]["scripts_dir"], "add_cluster_info_to_gene_index.py")

# Calculate auxiliary scores of protein clusters and generate organized output tables
rule aux_scores:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "add_cluster_info_to_gene_index.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "aux_scores.done"))
    params:
        gene_index = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        cluster_taxa_levels = config["cluster_ranks"],
        aux_scores = os.path.join(config["paths"]["output_dir"], "results", "aux_scores.tsv"),
        exclude_singletons = bool(config["exclude_singletons"]),
        outdir = os.path.join(config["paths"]["output_dir"], "results"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        """
        Calculating auxiliary scores of proteins in {params.gene_index} and wiring to {params.aux_scores}.
        """
    script:
        os.path.join(config["paths"]["scripts_dir"], "aux_scores.py")

# Extract the protein clusters and their proteins that remained after filtering by aux_scores.py
rule extract_filtered_clusters:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "aux_scores.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "extract_filtered_clusters.done"))
    params:
        aux_scores = os.path.join(config["paths"]["output_dir"], "results", "aux_scores.tsv"),
        outdir = os.path.join(config["paths"]["output_dir"], "results", "protein_sequences"),
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        vmag_proteins_subdir = directory(os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds", "vMAG_proteins")),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        """
        Extracting filtered protein clusters.
        """
    script:
        os.path.join(config["paths"]["scripts_dir"], "extract_filtered_clusters.py")

# Organize proteins by viral origin confidence
rule organize_proteins:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "extract_filtered_clusters.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "organize_proteins.done"))
    params:
        build_or_annotate = "build",
        metabolism_table = os.path.join(config["paths"]["output_dir"], "results", "metabolic_genes_curated.tsv"),
        physiology_table = os.path.join(config["paths"]["output_dir"], "results", "physiology_genes_curated.tsv"),
        regulation_table = os.path.join(config["paths"]["output_dir"], "results", "regulation_genes_curated.tsv"),
        all_genes_annotated = os.path.join(config["paths"]["output_dir"], "results", "gene_annotations.tsv"),
        aux_scores = os.path.join(config["paths"]["output_dir"], "results", "aux_scores.tsv"),
        outdir = os.path.join(config["paths"]["output_dir"], "results", "protein_sequences"),
        ranks = config["cluster_ranks"],
        input_prot_subdir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Organizing proteins into auxiliary & metabolic, auxiliary not metabolic, and metabolic not auxiliary categories."
    script:
        os.path.join(config["paths"]["scripts_dir"], "organize_proteins.py")

rule cluster_filtered_genes:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "organize_proteins.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "cluster_quality_filtered_prots.done"))
    params:
        input_prot_dir = os.path.join(config["paths"]["output_dir"], "results", "protein_sequences"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "mmseqs_filtered_clusters"),
        cluster_taxa_levels = config["cluster_ranks"],
        confidence_levels = config["confidence_levels"],
        min_seq_id = 0,
        cov_fraction = 0.8,
        min_score = 0.3,
        sensitivity = config["mmseqs_params"]["sensitivity"],
        debug = bool(config["debug"]),
        log = config["log"]
    resources:
        mem = config["mem_limit"]
    threads:
        config["threads"]
    message:
        "Clustering auxiliary proteins with MMseqs"
    script:
        os.path.join(config["paths"]["scripts_dir"], "cluster_quality_filtered_prots.py")

rule align_protein_clusters:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "cluster_quality_filtered_prots.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "align_protein_clusters.done"))
    params:
        cluster_file = os.path.join(config["paths"]["output_dir"], "wdir", "mmseqs_filtered_clusters", "filtered_clusters.tsv"),
        dblookup = os.path.join(config["paths"]["output_dir"], "wdir", "mmseqs_filtered_clusters", "filtered_prots.db.lookup"),
        input_prot_dir = os.path.join(config["paths"]["output_dir"], "results", "protein_sequences"),
        confidence_levels = config["confidence_levels"],
        acc_prefix = config["hmm_acc_prefix"],
        prot_clust_to_accession = os.path.join(config["paths"]["output_dir"], "wdir", "mmseqs_filtered_clusters", "cluster_names.tsv"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir", "protein_cluster_msa"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        """
        Aligning proteins in each auxiliary protein cluster.
        """
    script:
        os.path.join(config["paths"]["scripts_dir"], "align_protein_clusters.py")

# Build HMM profiles for each protein cluster
rule build_hmm_profiles:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "align_protein_clusters.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "build_hmm_profiles.done"))
    params:
        outdir = os.path.join(config["paths"]["output_dir"], "results", "hmm_profiles"),
        alignments_dir = os.path.join(config["paths"]["output_dir"], "wdir", "protein_cluster_msa", "aligned"),
        alignments_dir_parent = os.path.join(config["paths"]["output_dir"], "wdir", "protein_cluster_msa"),
        prot_clust_to_accession = os.path.join(config["paths"]["output_dir"], "wdir", "mmseqs_filtered_clusters", "cluster_names.tsv"),
        output_table = os.path.join(config["paths"]["output_dir"], "results", "hmms.tsv"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Building HMM profiles for each scored protein cluster."
    script:
        os.path.join(config["paths"]["scripts_dir"], "build_hmm_profiles.py")

# Make the final summarized table with annotations, classifications, genomic context, auxiliary status, and HMM profile IDs
rule make_final_table:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "build_hmm_profiles.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "make_final_table.done"))
    params:
        build_or_annotate = "build",
        aux_scores = os.path.join(config["paths"]["output_dir"], "results", "aux_scores.tsv"),
        all_genes_annotated = os.path.join(config["paths"]["output_dir"], "results", "gene_annotations.tsv"),
        gene_index = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        hmm_table = os.path.join(config["paths"]["output_dir"], "results", "hmms.tsv"),
        metabolism_table = os.path.join(config["paths"]["output_dir"], "results", "metabolic_genes_curated.tsv"),
        physiology_table = os.path.join(config["paths"]["output_dir"], "results", "physiology_genes_curated.tsv"),
        regulation_table = os.path.join(config["paths"]["output_dir"], "results", "regulation_genes_curated.tsv"),
        final_table = os.path.join(config["paths"]["output_dir"], "results", "final_results.tsv"),
        hmm_results = os.path.join(config["paths"]["output_dir"], "wdir", "hmm_results.tsv"),
        acc_prefix = config["hmm_acc_prefix"],
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Creating the final summarized table with annotations, classifications, genomic context, auxiliary status, and HMM profile IDs."
    script:
        os.path.join(config["paths"]["scripts_dir"], "make_final_table.py")
