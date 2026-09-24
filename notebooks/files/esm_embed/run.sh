#!/bin/bash

FAA=$1
OUTDIR=$2
MODELESM=$3
ACCEL=$4
NDEVICES=$5

echo 'Date: ' `date`
echo 'Host: ' `hostname`
echo 'System: ' `uname -spo`
echo 'CUDA_VISIBLE_DEVICES: ' $CUDA_VISIBLE_DEVICES
echo 'Gpu: ' `nvidia-smi -L | grep $CUDA_VISIBLE_DEVICES`

CHECKPOINTSTAR="checkpoints.tar.gz"

cat <<EOF
Train FAA file:           $FAA
Output directory:         $OUTDIR
ESM model:                $MODELESM
Accelerator:              $ACCEL
Number of devices:        $NDEVICES
ESM2 checkpoints archive: $CHECKPOINTSTAR
EOF

set -e
ENVNAME="pst" # conda-packed PST environment name
TARBALL="${ENVNAME}.tar.gz"
ENVDIR=$ENVNAME

# Set DDP debug info
export NCCL_DEBUG="INFO"
export TORCH_CPP_LOG_LEVEL="INFO"
export TORCH_DISTRIBUTED_DEBUG="INFO"
export CUDA_LAUNCH_BLOCKING="1"

# move data

echo "Moving $FAA over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$FAA .

echo "Moving $CHECKPOINTSTAR over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$CHECKPOINTSTAR .
echo "Untarring $CHECKPOINTSTAR"
tar -xzf $CHECKPOINTSTAR

CHECKPOINTS="checkpoints"

# CONDA
echo "Moving $TARBALL over from $STAGING/$USER/conda_envs" 
cp $STAGING/$USER/conda_envs/$TARBALL .
export PATH
mkdir $ENVDIR
tar -xzf $TARBALL -C $ENVDIR
. $ENVDIR/bin/activate
rm $TARBALL

# untar checkpoints folder, pass the folder to --model
run_pst_embed () {
    pst embed \
        --input $FAA \
        --outdir $OUTDIR \
        --model.esm $MODELESM \
        --model.torch_hub . \
        --trainer.accelerator $ACCEL \
        --trainer.devices $NDEVICES
}

clean_up () {
    rm -rf $FAA $ENVDIR $CHECKPOINTS $CHECKPOINTSTAR

    if [ -d "$OUTDIR" ]
    then
	tar -czf $OUTDIR.tar.gz $OUTDIR && mv $OUTDIR.tar.gz $GRPSTAG/$USER
	rm -rf $OUTDIR
    fi
}

# run tool and ALWAYS cleanup
run_pst_embed && clean_up || clean_up