#!/usr/bin/bash
#SBATCH --time=120:00:00
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --account=liang
#SBATCH --qos=liang_noc
#SBATCH --job-name pm6fomo
#SBATCH --output=slurm-out
#SBATCH --error=slurm-error
#SBATCH --mem=32GB
cd $SLURM_SUBMIT_DIR

##./fix_delay.sh

rm -r failed
rm -r dissociated
mkdir failed
mkdir dissociated
rm -r allfolder
mkdir allfolder
rm -r yesfolder
mkdir yesfolder
cp 1.xyz yesfolder/
rm -r nofolder 
mkdir nofolder
rm -r newlat
mkdir newlat
rm -r climbinglat
mkdir climbinglat
rm -r nofoldertmp
mkdir nofoldertmp
rm -r prefolder
mkdir prefolder

rm filesinf*.txt
rm shell*
rm set0*
rm foundp*
rm failed/*
rm pre*
rm allfiles*
rm converted*
rm dissociated/*
rm climbing*
rm new_lattice*
rm dissociated.xyz

./a.out

