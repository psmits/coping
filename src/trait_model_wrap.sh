#!/bin/bash
FILES=../data/data_dump/trait_info*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/trait_categorical sample num_samples=1000 num_warmup=1000 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/trait_categorical_out_${i}.csv &
  done
  wait
done
