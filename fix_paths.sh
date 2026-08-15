#!/usr/bin/bash
# fix_paths.sh - point the run scripts at a local conda environment and MOPAC-PI build.
#
# Run once after cloning.
#
# Usage:  CONDA_ENV=myenv MOPACPI=/path/to/mopacpi.x ./fix_paths.sh
set -eu

CONDA_ENV=${CONDA_ENV:-gtc_mopac}
MOPACPI=${MOPACPI:-/home/rliang/mopacpi/mopacpi.x}
CONDA_SH='$(conda info --base)/etc/profile.d/conda.sh'

cd "$(dirname "$0")"
n=0
for f in $(grep -rl 'python3p8\|mopacpi-master\|anaconda3/etc/profile.d/conda.sh' . \
             --include='*.sh' --include='*.py' 2>/dev/null); do
    sed -i "s#conda activate python3p8#conda activate ${CONDA_ENV}#g; \
            s#/home/pan60047/mopacpi-master/mopacpi.x#${MOPACPI}#g; \
            s#source /home/pan60047/anaconda3/etc/profile.d/conda.sh#source ${CONDA_SH}#g" "$f"
    n=$((n + 1))
done
echo "repointed $n scripts  (env=${CONDA_ENV}, mopacpi=${MOPACPI})"
