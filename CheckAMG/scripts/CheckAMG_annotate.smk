import os

input_type = config["input_type"]

rule all:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "make_final_table.done")

if input_type == "nucl":
    # Filter input sequences by length
    rule filter_by_length:
        output:
            touch(os.path.join(config["paths"]["output_dir"], "wdir", "filter_by_length.done"))
        params:
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
            input_type = config["input_type"],
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

elif input_type == "prot":
    # Filter input amino-acid sequences by minimum number of CDS
    rule filter_by_cds:
        output:
            touch(os.path.join(config["paths"]["output_dir"], "wdir", "filter_by_cds.done"))
        params:
            input_type = config["input_type"],
            input_single_contig_prots = config["input_single_contig_prots"],
            input_vmag_prots = config["input_vmag_prots"],
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

else:
    raise ValueError("Invalid input_type: {input_type}")

# Assign functional annotations to the proteins in the database
rule assign_annots:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "filter_by_cds.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "annotate_hmm.done"))
    params:
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        hmm_vscores = os.path.join(config["paths"]["files_dir"], "vscores.csv"),
        vscores = os.path.join(config["paths"]["output_dir"], "wdir", "vscores.tsv"),
        all_hmm_results = os.path.join(config["paths"]["output_dir"], "wdir", "hmm_results.tsv"),
        wdir = os.path.join(config["paths"]["output_dir"], "wdir"),
        db_dir = config["paths"]["db_dir"],
        cov_fraction = config["cov_fraction"],
        min_bitscore = config["min_bitscore"],
        min_bitscore_fraction_heuristic = config["min_bitscore_fraction_heuristic"],
        max_evalue = config["max_evalue"],
        kegg_cutoff_file = os.path.join(config["paths"]["files_dir"], "KEGG_cutoffs.csv"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Assigning V-scores and L-scores to proteins in {params.protein_dir} using an HMMsearch of the annotations in {params.hmm_vscores}"
    script:
        os.path.join(config["paths"]["scripts_dir"], "annotate_hmm.py")


# Obtain gene information from input (prodigal-formatted) .faa and genome information from...
rule index_genes:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "annotate_hmm.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "index_genes.done"))
    params:
        cluster_taxa_levels = None,
        gene_index = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index.tsv"),
        vscores = os.path.join(config["paths"]["output_dir"], "wdir", "vscores.tsv"),
        all_hmm_results = os.path.join(config["paths"]["output_dir"], "wdir", "hmm_results.tsv"),
        out_parent = os.path.join(config["paths"]["output_dir"], "wdir"),
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        vmag_proteins_subdir = directory(os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds", "vMAG_proteins"),),
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
        "Adding calculated v-scores in {params.vscores} to {params.gene_index} and writing to {params.gene_index_annotated}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "add_annots.py")

# Analyze the genomic context of annotations
rule genome_context:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "add_annots.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "genome_context.done"))
    params:
        outparent = os.path.join(config["paths"]["output_dir"], "results"),
        context_table = os.path.join(config["paths"]["output_dir"], "results", "genes_genomic_context.tsv"),
        gene_index_annotated = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index_annotated.tsv"),
        circular_contigs = os.path.join(config["paths"]["output_dir"], "wdir", "circular_contigs.tsv"),
        annotation_percent_threshold = config["annotation_percent_threshold"],
        window_size = config["window_size"],
        minimum_flank_vscore = config["minimum_flank_vscore"],
        max_flank_length = config["max_flank_length"],
        use_hallmark = config["use_hallmark"],
        hallmark_path = os.path.join(config["paths"]["files_dir"], "viral_hallmark_genes.csv"),
        mobile_genes_path = os.path.join(config["paths"]["files_dir"], "mobile_genes.csv"),
        vscore_ref = os.path.join(config["paths"]["files_dir"], "vscores.csv"),
        lgbm_model = os.path.join(config["paths"]["files_dir"], "lgbm_model.joblib"),
        feature_names = os.path.join(config["paths"]["files_dir"], "lgbm_feature_names.joblib"),
        thresholds = os.path.join(config["paths"]["files_dir"], "lgbm_thresholds.joblib"),
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
        context_table = os.path.join(config["paths"]["output_dir"], "results", "genes_genomic_context.tsv"),
        metabolism_table = os.path.join(config["paths"]["files_dir"], "AMGs.tsv"),
        physiology_table = os.path.join(config["paths"]["files_dir"], "APGs.tsv"),
        regulation_table = os.path.join(config["paths"]["files_dir"], "AReGs.tsv"),
        metabolism_table_out = os.path.join(config["paths"]["output_dir"], "results", "metabolic_genes_curated.tsv"),
        physiology_table_out = os.path.join(config["paths"]["output_dir"], "results", "physiology_genes_curated.tsv"),
        regulation_table_out = os.path.join(config["paths"]["output_dir"], "results", "regulation_genes_curated.tsv"),
        soft_keyword_bypass_scaling_factor = config["soft_keyword_bypass_scaling_factor"],
        cov_fraction = config["cov_fraction"],
        min_bitscore = config["min_bitscore"],
        kegg_cutoff_file = os.path.join(config["paths"]["files_dir"], "KEGG_cutoffs.csv"),
        foam_cutoff_file = os.path.join(config["paths"]["files_dir"], "FOAM_cutoffs.csv"),
        all_annot_out_table = os.path.join(config["paths"]["output_dir"], "results", "gene_annotations.tsv"),
        hmm_ref = os.path.join(config["paths"]["files_dir"], "hmm_id_to_name.csv"),
        false_amgs = os.path.join(config["paths"]["files_dir"], "false_amgs.csv"),
        false_apgs = os.path.join(config["paths"]["files_dir"], "false_apgs.csv"),
        false_aregs = os.path.join(config["paths"]["files_dir"], "false_aregs.csv"),
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
        
# Organize proteins into auxiliary & metabolic, auxiliary not metabolic, and metabolic not auxiliary categories
rule organize_proteins:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "curate_results.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "organize_proteins.done"))
    params:
        metabolism_table = os.path.join(config["paths"]["output_dir"], "results", "metabolic_genes_curated.tsv"),
        physiology_table = os.path.join(config["paths"]["output_dir"], "results", "physiology_genes_curated.tsv"),
        regulation_table = os.path.join(config["paths"]["output_dir"], "results", "regulation_genes_curated.tsv"),
        all_genes_annotated = os.path.join(config["paths"]["output_dir"], "results", "gene_annotations.tsv"),
        protein_dir = os.path.join(config["paths"]["output_dir"], "wdir", "filtered_input", "filtered_faa_by_cds"),
        aux_fasta_dir = os.path.join(config["paths"]["output_dir"], "results"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Writing the final AMG results to {params.aux_fasta_dir}."
    script:
        os.path.join(config["paths"]["scripts_dir"], "organize_proteins.py")

# Make the final summarized table with annotations, genomic context, and classifications
rule make_final_table:
    input:
        os.path.join(config["paths"]["output_dir"], "wdir", "organize_proteins.done")
    output:
        touch(os.path.join(config["paths"]["output_dir"], "wdir", "make_final_table.done"))
    params:
        all_genes_annotated = os.path.join(config["paths"]["output_dir"], "results", "gene_annotations.tsv"),
        gene_index = os.path.join(config["paths"]["output_dir"], "wdir", "gene_index.tsv"),
        metabolism_table = os.path.join(config["paths"]["output_dir"], "results", "metabolic_genes_curated.tsv"),
        physiology_table = os.path.join(config["paths"]["output_dir"], "results", "physiology_genes_curated.tsv"),
        regulation_table = os.path.join(config["paths"]["output_dir"], "results", "regulation_genes_curated.tsv"),
        final_table = os.path.join(config["paths"]["output_dir"], "results", "final_results.tsv"),
        debug = bool(config["debug"]),
        log = config["log"]
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Creating the final summarized table with annotations, genomic context, and classifications."
    script:
        os.path.join(config["paths"]["scripts_dir"], "make_final_table.py")