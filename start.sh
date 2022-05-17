#!/bin/sh

for var in 125 200 500
do
    value=$(echo "scale=2;$var / 100" | bc)
    eval "export PARAM1=$value"
    if [ $var = 125 ]; then
      ID=$(sbatch --parsable --array=1-1000 launch.sh)
      ID=$(sbatch --parsable --dependency=afterany:${ID} recomb.sh)
      ID=$(sbatch --parsable --dependency=after:${ID} clean.sh)
    else
      ID=$(sbatch --parsable --dependency=after:${ID} --array=1-1000 launch.sh)
      ID=$(sbatch --parsable --dependency=afterany:${ID} recomb.sh)
      ID=$(sbatch --parsable --dependency=after:${ID} clean.sh)
    fi
done
