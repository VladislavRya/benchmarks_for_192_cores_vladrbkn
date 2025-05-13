#!/bin/bash

benchmarks=(cg jac3d bt_mz)

for benchmark in "${benchmarks[@]}"
do
  for prog_dir in "$benchmark"/*/
  do
    cd "$prog_dir"
    ./run_prog.sh
    cd -
  done
done