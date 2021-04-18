#!/bin/bash

all_ars=(5 10 15 20 25 30)
all_sigmas=(0.88)
all_betas=(0.1 0.3 0.5)

# array from 0.1 to 0.5 with step 0.2, it will be equal to (0.1 0.3 0.5)
#all_betas=`seq 0.1 0.2 0.5`

# array from 1 to 10, inclusive
all_reps=`seq 1 50`


for ar in ${all_ars[@]}; do
  for sigma in ${all_sigmas[@]}; do
    for beta in ${all_betas[@]}; do  
      for rep in ${all_reps[@]}; do
        echo $ar $sigma $beta $rep
        time nice ./runsim -ar $ar -sigma $sigma -beta $beta -rep $rep -boost -boostfactor 7.26
      done
    done
  done
done
