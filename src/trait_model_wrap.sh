#!/bin/bash
FILES=../data/data_dump/trait_info*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/turnover_horseshoe sample num_samples=5000 num_warmup=5000 thin=5 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=../data/mcmc_out/turnover_${i}.csv &
  done
  wait
done
