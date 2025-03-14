#!/bin/bash
#
# Example data staging job script that copies a directory from DataStore to Eddie with rsync
# 
# Job will restart from where it left off if it runs out of time 
# (so setting an accurate hard runtime limit is less important)

# Grid Engine options start with a #$
# Name job and set to use current working directory
#$ -N stagein 
#$ -cwd
# Choose the staging environment
#$ -q staging

# Hard runtime limit
#$ -l h_rt=12:00:00 

# Make the job resubmit itself if it runs out of time: rsync will start where it left off
#$ -r yes
#$ -notify
trap 'exit 99' sigusr1 sigusr2 sigterm

# Source and destination directories
#
# Source path on DataStore in the staging environment
# Note: these paths are only available on the staging nodes
# It should start with one of /exports/csce/datastore, /exports/chss/datastore, /exports/cmvm/datastore or /exports/igmm/datastore
#
SOURCE=/exports/cmvm/datastore/eb/groups/macqueen_lab/JS/2023_chicken_organoid
#
# Destination path on Eddie. It should be on the fast Eddie HPC filesystem, starting with one of:
# /exports/csce/eddie, /exports/chss/eddie, /exports/cmvm/eddie, /exports/igmm/eddie or /exports/eddie/scratch, 
#
DESTINATION=/exports/eddie/scratch/s1964742

# Perform copy with rsync
# Note: do not use -p or -a (implies -p) as this can break file ACLs at the destination
rsync -rl ${SOURCE} ${DESTINATION}

cd /exports/eddie/scratch/s1964742/2023_chicken_organoid
source /etc/profile
tar -xf "/exports/eddie/scratch/s1964742/2023_chicken_organoid/LV_230323.tar.gz"
tar -xf "/exports/eddie/scratch/s1964742/2023_chicken_organoid/LV_240323.tar.gz"
