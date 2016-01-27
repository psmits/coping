#!/bin/bash
FILES=../data/data_dump/trait_info*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/trait_model sample num_samples=2000 num_warmup=2000 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/trait_${i}.csv &
  done
  wait
done
