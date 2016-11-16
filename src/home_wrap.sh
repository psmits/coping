../stan/turnover_revamp variational \
  id=1 \
  init=0 \
  data file=../data/data_dump/trait_w_gaps_revamp.data.R \
  output file=../data/mcmc_out/turnover_revamp_advi.csv &
wait
