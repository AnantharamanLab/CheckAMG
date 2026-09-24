#!/bin/bash

INPUT_FAA=$1
TRAIN_INDEX=$2
TRAIN_LABELS=$3
CKPT=$4
OUTDIR=$5
KNN=$6

# simple executable for CHTC
echo 'Date: ' `date`
echo 'Host: ' `hostname`
echo 'System: ' `uname -spo`
echo 'CUDA_VISIBLE_DEVICES: ' $CUDA_VISIBLE_DEVICES
echo 'Gpu: ' `nvidia-smi -L | grep $CUDA_VISIBLE_DEVICES`

CHECKPOINTSTAR="checkpoints.tar.gz"

cat <<EOF
Input data file:            $INPUT_FAA
Train index file:           $TRAIN_INDEX
Train labels file:          $TRAIN_LABELS
Model checkpoint:           $CKPT
Output directory:           $OUTDIR
ESM2 checkpoints archive:   $CHECKPOINTSTAR
KNN:                        $KNN
EOF

set -e
ENVNAME="checkamg_1.0_gpu"
TARBALL="${ENVNAME}.tar.gz"
ENVDIR=$ENVNAME

# Set DDP debug info
export NCCL_DEBUG="INFO"
export TORCH_CPP_LOG_LEVEL="INFO"
export TORCH_DISTRIBUTED_DEBUG="INFO"
export CUDA_LAUNCH_BLOCKING="1"

# move data
echo "Moving $INPUT_FAA over from $GRPSTAG/$USER/denovo_splits"
cp $GRPSTAG/$USER/denovo_splits/$INPUT_FAA .

echo "Moving $CKPT over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$CKPT .

echo "Moving $TRAIN_INDEX over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$TRAIN_INDEX .

echo "Moving $TRAIN_LABELS over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$TRAIN_LABELS .

echo "Moving $CHECKPOINTSTAR over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$CHECKPOINTSTAR .

echo "Untarring $CHECKPOINTSTAR"
tar -xzf $CHECKPOINTSTAR

CHECKPOINTS="checkpoints"

# conda
echo "Moving $TARBALL over from $STAGING/$USER/conda_envs"
cp $STAGING/$USER/conda_envs/$TARBALL .
export PATH
mkdir $ENVDIR
tar -xzf $TARBALL -C $ENVDIR
. $ENVDIR/bin/activate
rm $TARBALL

run_checkamg_denovo () {
    checkamg de-novo \
        --query-proteins "$INPUT_FAA" \
        --train-index-file "$TRAIN_INDEX" \
        --train-labels-file "$TRAIN_LABELS" \
        --model-ckpt "$CKPT" \
        --output "$OUTDIR" \
        --esm2-ckpt-dir $CHECKPOINTS \
        --knn $KNN \
        -a gpu \
        --devices 1 \
        --mem 1000
}

clean_up () {
    rm -rf $INPUT_FAA $TRAIN_INDEX $TRAIN_LABELS $CKPT $ENVDIR $CHECKPOINTS $CHECKPOINTSTAR

    if [ -d "$OUTDIR" ]
    then
	tar -czf $OUTDIR.tar.gz $OUTDIR && mv $OUTDIR.tar.gz $GRPSTAG/$USER
	rm -rf $OUTDIR
    
    fi
}

# run tool and ALWAYS cleanup
echo "Running checkamg with command:"
echo "checkamg de-novo \
        --query-proteins "$INPUT_FAA" \
        --train-index-file "$TRAIN_INDEX" \
        --train-labels-file "$TRAIN_LABELS" \
        --model-ckpt "$CKPT" \
        --output "$OUTDIR" \
        --esm2-ckpt-dir $CHECKPOINTS \
        --knn $KNN \
        -a gpu \
        --devices 1 \
        --mem 1000"
    
run_checkamg_denovo && clean_up || clean_up