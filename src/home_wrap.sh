#../stan/turnover_revamp variational \
#  id=1 \
#  init=0 \
#  data file=../data/data_dump/trait_w_gaps.data.R \
#  output file=../data/mcmc_out/turnover_revamp_1_advi.csv &
#wait
../stan/turnover_rwprior variational \
  id=1 \
  init=0 \
  data file=../data/data_dump/trait_w_gaps_NALMA.data.R \
  output file=../data/mcmc_out/turnover_rwprior_advi_NALMA.csv &
wait
