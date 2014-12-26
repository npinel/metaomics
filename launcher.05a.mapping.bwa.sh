#! /bin/bash
################################################################################
# parallel_launcher.sh -  Example of a generic launcher script using
#     [GNU Parallel](http://www.gnu.org/software/parallel/) able to run a
#     program  across reserved nodes. 
#     
# Submit this job in passive mode by 
#
#   oarsub [options] -S ./parallel_launcher.sh
################################################################################

##########################
#   The OAR  directives  #
##########################
# do NOT remove the # from lines containing OAR directives
#
#OAR -l nodes=8/core=12,walltime=8:00:00
#OAR -n bwa_batch
#OAR --stdout launcher-bwa-asm1MGasm2MT-%jobid%.log
#OAR --stderr launcher-bwa-asm1MGasm2MT-%jobid%.log

#####################################
#                                   #
#   The UL HPC specific directives  #
#                                   #
#####################################
if [ -f  /etc/profile ]; then
    .  /etc/profile
fi

# load ESB-specific shortcuts
source /mnt/nfs/projects/ecosystem_biology/esb.shortcuts

# Modules to preload
MODULE_TO_LOAD=('BWA/0.7.4-goolf-1.4.10' 'SAMtools/0.1.18-goolf-1.4.10')

# Characteristics of the reservation
NB_CORES=`cat ${OAR_NODEFILE} | wc -l`
NB_HOSTS=`cat ${OAR_NODEFILE} | uniq | wc -l`

# The [serial] task to be executed
TASK='/home/users/npinel/git/laoPipeline/05a.mapping.bwa.sh'

# Define here a file containing the arguments to pass to the task, one line per 
# exected run
ARG_TASK_FILE='/home/users/npinel/lao_time_series/params.bwa_launcher.txt'

# Total number of tasks to be executed
[ -n "${ARG_TASK_FILE}" ] && NB_TASKS=`cat ${ARG_TASK_FILE} | wc -l` || NB_TASKS=$(( 2*${NB_HOSTS} ))

# Number of concurrent cores that have to be used on a node to perform a single task
NB_CORE_PER_TASK=12
NB_JOBS=$NB_HOSTS

#####################################
#   The GNU parallel directives     # 
#####################################
# File with sshlogins. The file consists of sshlogins on separate lines. 
GP_SSHLOGINFILE=/tmp/gnuparallel_hostfile.${OAR_JOBID}

# Eventually drop here the options you want to pass to GNU parallel
GP_OPTS=

################# Let's go ###############

# DIRECTORY WHERE TO RUN 
cd $LAOWD

# Prepare an sshloginfile for GNU parallel to define connection to remote nodes.
# 3 versions are defined here:
#    1. ${GP_SSHLOGINFILE}.core : each line correspond exactly to 1 core on a node
#    2. ${GP_SSHLOGINFILE}.node : each line correspond exactly to 1 node, thus
#       of the format  
#           <#cores>/oarsh <hostname>
#    3. {GP_SSHLOGINFILE}.task : each line correspond exactly to the
#       resource of one task as defined by ${NB_CORE_PER_TASK}, thus of the format  
#           <#core_per_task>/oarsh <hostname>

cat $OAR_NODEFILE | awk '{printf "oarsh %s\n", $1}'  > ${GP_SSHLOGINFILE}.core
GP_SSHLOGIN_OPT1=`cat $OAR_NODEFILE | awk '{printf "--sshlogin \"oarsh %s\" ", $1}'`

cat $OAR_NODEFILE | uniq -c | awk '{printf "%s/oarsh %s\n", $1, $2}' > ${GP_SSHLOGINFILE}.node
GP_SSHLOGIN_OPT2=`cat $OAR_NODEFILE | uniq -c | awk '{printf "--sshlogin \"%s/oarsh %s\" ", $1, $2}'`

cat $OAR_NODEFILE | uniq -c | while read line; do
    NB_CORE=`echo $line  | awk '{ print $1 }'`
    HOSTNAME=`echo $line | awk '{ print $2 }'`
    n=$(( ${NB_CORE}/${NB_CORE_PER_TASK} ))
    SSHLOGIN="$NB_CORE_PER_TASK/oarsh $HOSTNAME"
    while [ $(( n -= 1 )) -ge 0 ]; do 
        echo "${SSHLOGIN}" >> ${GP_SSHLOGINFILE}.task
        GP_SSHLOGIN_OPT="${GP_SSHLOGIN_OPT} --sshlogin '${SSHLOGIN}'"
    done
done

# Load the required modules
if [ -z "${ARG_TASK_FILE}" ]; then 
    # ============
    #  Example 1:
    # ============
    # use GNU parallel to perform the tasks on the reserved nodes: 
    #    ${TASK} 1
    #    ${TASK} 2
    #    [...]
    #    ${TASK} ${NB_TASKS}
    seq ${NB_TASKS} | parallel --tag -u -j ${NB_JOBS} --sshloginfile ${GP_SSHLOGINFILE}.task ${GP_OPTS} ${TASK} {}
else 
    # ============
    #  Example 2:
    # ============
    # use GNU parallel to perform the tasks on the nodes to run in
    # parallel on ${NB_CORE_PER_TASK} cores for each line of
    # ${ARG_TASK_FILE} : 
    #    ${TASK} <line1>
    #    ${TASK} <line2>
    #    [...]
    #    ${TASK} <lastline>
    for m in ${MODULE_TO_LOAD[*]}; do
	module load $m
    done
    cat ${ARG_TASK_FILE} | parallel --tag -u -j ${NB_JOBS} --sshloginfile ${GP_SSHLOGINFILE}.task --colsep ' ' ${GP_OPTS} ${TASK} {}
fi


# Cleanup (not mandatory as OAR will clean /tmp on job termination)
# [ -f "${GP_SSHLOGINFILE}.core" ] && rm -f ${GP_SSHLOGINFILE}.core
# [ -f "${GP_SSHLOGINFILE}.node" ] && rm -f ${GP_SSHLOGINFILE}.node
# [ -f "${GP_SSHLOGINFILE}.task" ] && rm -f ${GP_SSHLOGINFILE}.task

