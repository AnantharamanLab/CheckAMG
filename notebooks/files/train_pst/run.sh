#!/bin/bash

TRAINDATA=$1
OUTDIR=$2
MARGINS=$3
LR=$4
CKPT=$5
EPOCHS=$6
BATCHSIZE=$7
APPAIRS=$8
POSMINE=$9
NEGMINE=${10}
OPPNEGS=${11}
CLASSWEIGHT=${12}
CONTEXTSIZE=${13}

# simple executable for CHTC
echo 'Date: ' `date`
echo 'Host: ' `hostname`
echo 'System: ' `uname -spo`
echo 'CUDA_VISIBLE_DEVICES: ' $CUDA_VISIBLE_DEVICES
echo 'Gpu: ' `nvidia-smi -L | grep $CUDA_VISIBLE_DEVICES`

cat <<EOF
Train data file:            $TRAINDATA
Output directory:           $OUTDIR
Margins:                    $MARGINS
Learning rate:              $LR
Model checkpoint:           $CKPT
Max number of epochs:       $EPOCHS
Batch size:                 $BATCHSIZE
Max number of AP pairs:     $APPAIRS
Positive mining strategy:   $POSMINE
Negative mining strategy:   $NEGMINE
Opposite class negatives:   $OPPNEGS
Class weighing strategy:    $CLASSWEIGHT
Context size:               $CONTEXTSIZE
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
echo "Moving $TRAINDATA over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$TRAINDATA .
echo "Moving $CKPT over from $GRPSTAG/$USER"
cp $GRPSTAG/$USER/$CKPT .

# conda
echo "Moving $TARBALL over from $STAGING/$USER/conda_envs"
cp $STAGING/$USER/conda_envs/$TARBALL .
export PATH
mkdir $ENVDIR
tar -xzf $TARBALL -C $ENVDIR
. $ENVDIR/bin/activate
rm $TARBALL

run_checkamg_train () {
    checkamg train \
        --train-file "$TRAINDATA" \
        --output "$OUTDIR" \
        --margin "$MARGINS" \
        --lr "$LR" \
        --model-ckpt "$CKPT" \
        --max-epochs "$EPOCHS" \
        --batch-size "$BATCHSIZE" \
        --max-ap-pairs "$APPAIRS" \
        --positive-mining-strategy "$POSMINE" \
        --negative-mining-strategy "$NEGMINE" \
        --"$OPPNEGS" \
        --class-weighting "$CLASSWEIGHT" \
        --context-size "$CONTEXTSIZE" \
        --mem 1000 \
        -a gpu \
        --save-train-embed \
        --verbose
}

clean_up () {
    rm -rf $DATA $CKPT $TRAINDATA $ENVDIR

    if [ -d "$OUTDIR" ]
    then
	tar -czf $OUTDIR.tar.gz $OUTDIR && mv $OUTDIR.tar.gz $GRPSTAG/$USER
	rm -rf $OUTDIR
    
    fi
}

# run tool and ALWAYS cleanup
echo "Running checkamg with command:"
echo "checkamg train \
        --train-file "$TRAINDATA" \
        --output "$OUTDIR" \
        --margin "$MARGINS" \
        --lr "$LR" \
        --model-ckpt "$CKPT" \
        --max-epochs "$EPOCHS" \
        --batch-size "$BATCHSIZE" \
        --max-ap-pairs "$APPAIRS" \
        --positive-mining-strategy "$POSMINE" \
        --negative-mining-strategy "$NEGMINE" \
        --"$OPPNEGS" \
        --class-weighting "$CLASSWEIGHT" \
        --context-size "$CONTEXTSIZE" \
        --mem 1000 \
        -a gpu \
        --save-train-embed \
        --verbose"
    
run_checkamg_train && clean_up || clean_up