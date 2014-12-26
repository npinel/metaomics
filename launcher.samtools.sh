#! /bin/bash
################################################################################
# Submit this job in passive mode by 
#
#   oarsub [options] -S ./launch.05a.mapping.bwa.sh
#
################################################################################

##########################
#   The OAR  directives  #
##########################
# do NOT remove the # from lines containing OAR directorives
# use these to specify oarsub options

#OAR -l nodes=1/core=12,walltime=0:30:00
#OAR -n sam.index
#OAR --stdout test_launcher
#OAR --stderr test_launcher

#####################################
#   The UL HPC specific directives  #
#####################################
if [ -f  /etc/profile ]; then
    .  /etc/profile
fi

# load shortcuts to capture the current location of $LAOTS
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

# Modules to preload - last updated 2013-10-31
# the current version of 05a.mapping.bwa.sh was written for bwa 0.7
# check compatibility if upgrading versions
MODULE_TO_LOAD=("SAMtools/0.1.18-goolf-1.4.10")

# program/task to run
TASK='samtools'

# Define here a file containing the arguments to pass to the task, one line per 
# expected run.
# ARG_TASK_FILE=$LAOTS/params/20131031.05a.mapping.bwa
ARG_TASK_FILE='/home/users/npinel/lao_time_series/params.samtools.index.txt'

# DIRECTORY WHERE TO RUN
# adjust as necessary
cd $LAOWD

# Characteristics of the reservation: number of cores on the first (and normally
# only one) node
NB_CORES_HEADNODE=`cat ${OAR_NODEFILE} | uniq -c | head -n1 | awk '{print $1}'`
# Default value
: ${NB_CORES_HEADNODE:=1}

# Total number of tasks to be executed
# NB_TASKS can be adjusted manually for parallel runs of the same task if needed
[ -n "${ARG_TASK_FILE}" ] && NB_TASKS=`cat ${ARG_TASK_FILE} | wc -l` || NB_TASKS=1

################# Let's go ###############
# Load the required modules
for m in ${MODULE_TO_LOAD[*]}; do 
    module load $m
done

if [ -z "${ARG_TASK_FILE}" ]; then 
    seq ${NB_TASKS} | parallel -u -j ${NB_CORES_HEADNODE} ${TASK} {}
else 
    cat ${ARG_TASK_FILE} | parallel -u -j ${NB_CORES_HEADNODE} --colsep ' ' ${TASK} {}
    # OR
    # parallel -u -j ${NB_CORES_HEADNODE} --colsep ' ' -a ${ARG_TASK_FILE} ${TASK} {}
fi

