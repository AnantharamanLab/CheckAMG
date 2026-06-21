import os
from pathlib import Path

# Embedding the training data is an optional final step. When requested, it runs
# as a separate rule so that a disruption after model training does not force a
# retrain: training output (train_pst_model.done) is preserved and only the
# embedding step is rerun.
train_embed_output = config["paths"]["train_embed_output"]
n_steps = 2 if train_embed_output else 1

rule all:
    input:
        [
            os.path.join(config["paths"]["output_dir"], "snakemake", "train_pst_model.done"),
            *(
                [os.path.join(config["paths"]["output_dir"], "snakemake", "embed_train_data.done")]
                if train_embed_output
                else []
            ),
        ]

rule train_pst_model:
    output:
        touch(os.path.join(config["paths"]["output_dir"], "snakemake", "train_pst_model.done"))
    params:
        # Inputs
        train_file = Path(config["paths"]["train_file"]) if config["paths"]["train_file"] else None,
        test_file = Path(config["paths"]["test_file"]) if config["paths"]["test_file"] else None,
        model_ckpt = Path(config["paths"]["model_ckpt"]) if config["paths"]["model_ckpt"] else None,
        # Outputs
        logdir = Path(config["paths"]["logdir"]) if config["paths"]["logdir"] else None,
        # Parameters
        batch_size = config["batch_size"],
        learning_rate = config["learning_rate"],
        max_epochs = config["max_epochs"],
        dropout = config["dropout"],
        margin = config["margin"],
        knn = config["knn"],
        max_ap_pairs = config["max_ap_pairs"],
        positive_mining_strategy = config["positive_mining_strategy"],
        negative_mining_strategy = config["negative_mining_strategy"],
        opposite_class_negatives = bool(config["opposite_class_negatives"]),
        class_weighting = config["class_weighting"],
        context_size = config["context_size"],
        # Resources
        accelerator = config["accelerator"],
        devices = config["devices"],
        # Other
        debug = bool(config["debug"]),
        seed = config["seed"],
        log = config["log"],
        step = (1, n_steps),
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Train a CheckAMG-PST model for de-novo mode"
    script:
        os.path.join(config["paths"]["steps_dir"], "train_pst_model.py")

# Embed the training proteins with the newly trained best model. Runs only when
# --save-train-embed is set (train_embed_output is configured) and depends on the
# completed training step so it can be restarted independently.
if train_embed_output:
    rule embed_train_data:
        input:
            os.path.join(config["paths"]["output_dir"], "snakemake", "train_pst_model.done")
        output:
            touch(os.path.join(config["paths"]["output_dir"], "snakemake", "embed_train_data.done"))
        params:
            # Inputs
            train_file = Path(config["paths"]["train_file"]) if config["paths"]["train_file"] else None,
            logdir = Path(config["paths"]["logdir"]) if config["paths"]["logdir"] else None,
            # Outputs
            train_embed_output = Path(train_embed_output),
            # Parameters
            batch_size = config["batch_size"],
            class_weighting = config["class_weighting"],
            # Resources
            accelerator = config["accelerator"],
            devices = config["devices"],
            # Other
            debug = bool(config["debug"]),
            log = config["log"],
            step = (2, 2),
        threads:
            config["threads"]
        resources:
            mem = config["mem_limit"]
        message:
            "Embed training proteins with the trained CheckAMG-PST model"
        script:
            os.path.join(config["paths"]["steps_dir"], "embed_train_data.py")
