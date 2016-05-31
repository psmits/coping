#!/bin/bash
FILES=/peter/coping/data/data_dump/trait_info*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    /peter/coping/stan/turnover_simple sample num_samples=5000 num_warmup=5000 thin=5 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=/peter/coping/data/mcmc_out/turnover_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    /peter/coping/stan/turnover_noenv sample num_samples=5000 num_warmup=5000 thin=5 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=/peter/coping/data/mcmc_out/turnover_noenv_${i}.csv &
  done
  wait
done
#!/bin/bash
FILES=/peter/coping/data/data_dump/trait_w_gaps*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    /peter/coping/stan/turnover_simple sample num_samples=5000 num_warmup=5000 thin=5 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=/peter/coping/data/mcmc_out/turnover_gap_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    /peter/coping/stan/turnover_noenv sample num_samples=5000 num_warmup=5000 thin=5 \
      id=$i \
      random seed=420 \
      data file=$f \
      output file=/peter/coping/data/mcmc_out/turnover_noenv_gap_${i}.csv &
  done
  wait
done
