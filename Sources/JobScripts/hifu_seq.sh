#!/bin/bash

cd ..

sonications=20

for i in {10..16}
do
    digits=12       # number of digits in seed
    a=$(date +%s)
    seed=$((a*RANDOM))

    while [ ${#seed} -lt 12 ]; do
        seed="${seed}$RANDOM"
    done

    mkdir Outputs/HIFU/SEQ/${sonications}/${i}
    LOG_DIR="Outputs/HIFU/SEQ/${sonications}/${i}/"
    ./hifuSeq ${seed} ${LOG_DIR}
done
