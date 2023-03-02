#!/bin/bash 
    
#SBATCH --time=144:00:00                #Requested walltime 
#SBATCH -c1
#SBATCH -A snic2020-3-4
#SBATCH --job-name="tb__u3__d0_63__dimeric__0_5"

cd /proj/polymer/users/x_ihosa/tb/u3__d0_63/dimeric/0_5/
#make
chmod +x ./run
./run

