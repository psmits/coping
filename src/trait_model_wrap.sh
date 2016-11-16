#!/bin/bash
FILES=../data/data_dump/trait_w_gaps_revamp.data.R
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/turnover_revamp sample num_samples=5000 num_warmup=5000 thin=5 \
      id=$i \
      init=0 \
      data file=$f \
      output file=../data/mcmc_out/turnover_revamp_${i}.csv &
  done
  wait
done
