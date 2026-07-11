# Functions to set up the molcas environment

# Limit the  number of threads to 4:
export OMP_NUM_THREADS=4
# location of the Chilton molcas scripts
export MOLCAS_ENV=/home/users/mfr24/molcas_env
# location of our extra molcas scripts
export MOLCAS_EXTRA_PATH=/home/users/mfr24/dev/molcas_scripts
export PATH=$MOLCAS_EXTRA_PATH:$PATH
# open-molcas
export PATH=/home/users/mfr24/open-molcas-26.06:$PATH
export MOLCAS=/home/users/mfr24/open-molcas-26.06
# On your own machine you can use /tmp for the workdir:
export MOLCAS_WORKDIR=/tmp
# On servers you can use more memory, such as 32 GB
export MOLCAS_MEM=16000
# On servers use /scratch/$USER for the workdir. 
# You will have to create this yourself on each server. 
# export MOLCAS_WORKDIR=/scratch/$USER
 
# Helper to clean paths 
path_remove() {
    local var_name=$1
    local path_to_remove=$2
    local var_val="${!var_name}"
    var_val=$(echo -n "$var_val" | awk -v RS=: -v ORS=: '$0 != "'$path_to_remove'"' | sed 's/:$//')
    export "$var_name"="$var_val"
}

# Setup for molcas environment
molcas_activate(){
    echo "Set up paths and virtual environment for molcas"
    local old_ps1="${PS1:-}"
    export VIRTUAL_ENV_DISABLE_PROMPT=1
    # load  molcas virtual environment 
    source "$MOLCAS_ENV/bin/activate"
    unset VIRTUAL_ENV_DISABLE_PROMPT
    export _OLD_VIRTUAL_PS1="$old_ps1"
    PS1="(molcas) ${old_ps1}"
    export PS1
    echo "Activated molcas - to deactivate type: 'deactivate'"
    echo "Running other python apps in this environment will give unpredictable results."
}   
 

