#!/bin/bash

cd ..

for i in {10..16}
do
    digits=12       # number of digits in seed
    a=$(date +%s)
    seed=$((a*RANDOM))

    while [ ${#seed} -lt 12 ]; do
        seed="${seed}$RANDOM"
    done

    mkdir Outputs/Bench/SEQ/${i}
    LOG_DIR="Outputs/Bench/SEQ/${i}/"
    ./benchSeq ${seed} ${LOG_DIR}
done
