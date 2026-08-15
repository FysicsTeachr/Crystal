source ~/.bashrc
export CUDA_VISIBLE_DEVICES=1
export TeraChem=/home/ruibinliang/terachem_06_2023/terachem/TeraChem
export LD_LIBRARY_PATH=/home/ruibinliang/terachem_06_2023/terachem/TeraChem/lib:/home/ruibinliang/anaconda3/lib:$LD_LIBRARY_PATH
export NBOEXE=/home/ruibinliang/terachem_06_2023/terachem/TeraChem/bin/nbo6.i4.exe
export OMP_NUM_THREADS=32 
export PATH=/home/ruibinliang/terachem_06_2023/terachem/TeraChem/bin:$PATH

source /home/ruibinliang/terachem_06_2023/terachem/TeraChem/SetTCVars.sh
export OpenMM=/home/ruibinliang/anaconda3/
export OPENMM_PLUGIN_DIR=$TeraChem/lib/plugins
$TeraChem/bin/terachem ttc.in > ttc.out 
