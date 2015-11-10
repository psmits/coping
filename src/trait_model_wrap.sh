#!/bin/bash
FILES=../data/data_dump/coping_info*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/coping sample num_samples=6000 num_warmup=6000 thin=6 \
      init=0 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/coping_out_${i}.csv &
  done
  wait
done
