#!/bin/bash
for i in `seq 1 4`;
do
  ../stan/turnover_rwprior sample \
    adapt delta=0.95 \
    num_samples=5000 num_warmup=5000 thin=5 \
    algorithm=hmc engine=nuts max_depth=15 \
    id=$i \
    init=0 \
    data file=../data/data_dump/trait_w_gaps_NALMA.data.R \
    output file=../data/mcmc_out/turnover_rwprior_${i}.csv &
done
