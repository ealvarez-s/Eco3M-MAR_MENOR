# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=
ulimit -s unlimited
module purge
#SYMPHONIE
module load intel/18.2 intelmpi/18.2 hdf5/1.10.2-intelmpi netcdf/4.6.1-intelmpi pnetcdf/1.9.0-intelmpi

# User specific aliases and functions
module load nco/4.7.5 

# >>> conda initialize >>>
## !! Contents within this block are managed by 'conda init' !!
#__conda_setup="$('/usr/local/intel/2020.0.015/intelpython3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    eval "$__conda_setup"
#else
#    if [ -f "/usr/local/intel/2020.0.015/intelpython3/etc/profile.d/conda.sh" ]; then
#        . "/usr/local/intel/2020.0.015/intelpython3/etc/profile.d/conda.sh"
#    else
#        export PATH="/usr/local/intel/2020.0.015/intelpython3/bin:$PATH"
#    fi
#fi
#unset __conda_setup
# <<< conda initialize <<<

#Ferret
alias ferret='/users/p13120/duhaut/softs/ferret/bin/ferret_v7.2'
source /users/p13120/duhaut/softs/ferret/ferret_paths.sh

