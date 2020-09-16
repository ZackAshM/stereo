#!/usr/bin/env bash

#################
### ARG PARSE ###
#################

FILENAME=$(basename $0)
ABS_PATH=$(realpath $0)
PROGRAM_DIR=$(dirname $ABS_PATH)
DATA_PATH="$(dirname $PROGRAM_DIR)/data/images"

display_help() {
    arg_sty() {
	retval=$(tput bold)"${1}"$(tput sgr0)
	shift
	if [ "$#" -gt 0 ]; then
	    retval="${retval} "$(tput smul)"${@}"$(tput sgr0)
	fi
	echo $retval
    }

    echo "NAME"
    echo "   $FILENAME - Executes run_sim.py over a number of jobs and"
    echo "   saves the combination of results into one file called $(tput smul)filename$(tput sgr0)."
    echo
    echo "USAGE"
    echo "   $(arg_sty $FILENAME) [options] [$(arg_sty -h)] [$(arg_sty -t total)] [$(arg_sty -j jobs)] $(tput smul)filename$(tput sgr0)"
    echo
    echo "OPTIONS"
    echo "   $(arg_sty --help)"
    echo "   $(arg_sty -h)              Print this help message."
    echo
    echo "   $(arg_sty --total total)"
    echo "   $(arg_sty -t total)        The total number of solve trials to achieve."
    echo "                   Default=10000."
    echo
    echo "   $(arg_sty --jobs jobs)"
    echo "   $(arg_sty -j jobs)         How many jobs to split the sim into. Default=10."
    echo
    echo "   $(arg_sty --test_args)     Only print the arguments."
    exit 1
}

# saner programming env: these switches turn some bugs into errors
set -o errexit -o pipefail -o noclobber -o nounset

# -allow a command to fail with !’s side effect on errexit
# -use return value from ${PIPESTATUS[0]}, because ! hosed $?
! getopt --test > /dev/null 
if [[ ${PIPESTATUS[0]} -ne 4 ]]; then
    echo 'I’m sorry, `getopt --test` failed in this environment.'
    exit 1
fi

# define CL args
OPTIONS=ht:j:
LONGOPTS=help,test_args,total:,j:

# -regarding ! and PIPESTATUS see above
# -temporarily store output to be able to check for errors
# -activate quoting/enhanced mode (e.g. by writing out “--options”)
# -pass arguments only via   -- "$@"   to separate them correctly
! PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@")
if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    # e.g. return value is 1
    #  then getopt has complained about wrong arguments to stdout
    exit 2
fi
# read getopt’s output this way to handle the quoting right:
eval set -- "$PARSED"

# set defaults
TOTAL=10000 JOBS=10 TEST_ARGS=false
# now enjoy the options in order and nicely split until we see --
while true; do
    case "$1" in
	-h|--help)
	    display_help
	    exit 0
	    ;;
	--test_args)
	    TEST_ARGS=true
	    shift
	    ;;
	-t|--total)
	    TOTAL=$2
	    shift 2
	    ;;
	-j|--jobs)
	    JOBS=$2
	    shift 2
	    ;;
        --)
            shift
            break
            ;;
        *)
            echo "Programming error"
            exit 3
            ;;
    esac
done

# handle non-option arguments
if [[ $# -ne 1 ]]; then
    echo "$(basename $0): An output file name is required. Use -h or --help to see usage."
    exit 4
fi

# set required arguments
RUN_NAME=$1

if [ "$TEST_ARGS" = true ]; then
    echo "name: $RUN_NAME, total: $TOTAL, jobs: $JOBS"
    exit 0
fi

### END ARG PARSE ###


##############
### SCRIPT ###
##############

# create the file names
FILENAMES=""
TRIALS=$(( ($TOTAL + $JOBS - 1) / $JOBS ))
VERBOSE_FACTOR=$(( 2 * $TRIALS / 10 ))
for ((i=1; i<=$JOBS; i++)); do
    FILENAMES="${FILENAMES}_${RUN_NAME}${i} "
done
FILENAMES=${FILENAMES:-1}

# define how to finish the program
finish() {
    # combine the files into one
    echo "$FILENAME: Jobs finished. Combining data to $DATA_PATH/$RUN_NAME..."
    COMBINE="${PROGRAM_DIR}/combine.py"
    python3 $COMBINE _$RUN_NAME $RUN_NAME -s 1 -N $JOBS

    # delete the individual files
    FILENAMES=$(echo $FILENAMES | tr " " ,)
    eval "rm -rf $DATA_PATH/{$FILENAMES}.txt"

    echo "$FILENAME: Finished."

    exit 0
}

# call finish even if interrupted
trap finish SIGINT

# run the sim using all cores
RUN_SIM="${PROGRAM_DIR}/run_sim.py"
parallel --termseq INT,5000 -u -j 100% python3 $RUN_SIM -t $TRIALS -v $VERBOSE_FACTOR ::: $FILENAMES

exit 0

### END SCRIPT ###
