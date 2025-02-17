script_raw = """
#!/bin/bash
ulimit -s unlimited

re='^[0-9]+$' # regex for detecting numbers

if [ $# == 0 ] || ! [[ $1 =~ $re ]]; then
  echo -e "\nCRIMSON test suite runner script.\nUsage: ./mytest <numberOfProcessors> <options>\nTry option --gtest_help for info.\n"
  exit 1
fi


crimsonShuffleBool=0
crimsonRandomSeedProvidedOnCommandLine=0
for commandLineArg in "$@"; do
  case $commandLineArg in
    --gtest_shuffle*)
      echo -e "\n--gtest_shuffle is disabled due to random seed generation problems in parallel. Use --crimson_shuffle (and perhaps --crimson_random_seed=<num>) instead.\n"
      exit 1
      ;;
    --gtest_random_seed*)
      echo -e "\n--gtest_random_seed is disabled due to random seed generation problems in parallel. Use --crimson_random_seed instead.\n"
      exit 1
      ;;
    --crimson_shuffle*)
      crimsonShuffleBool=1
      ;;
    --crimson_random_seed=*)
      # Get the random seed from the command line:
      crimsonRandomSeed="${commandLineArg#*=}"
      crimsonRandomSeedProvidedOnCommandLine=1
      ;;
    *)
      ;;
  esac
done

# Generate a single random seed number and pass it to the test suite
# (avoids having each test thread generate its own random order for the tests,
# which would have consequences of great badness)
#
# Initialise empty commands:
shuffleCommand=
seedCommand=
if [ $crimsonShuffleBool == 1 ]; then
  if [ $crimsonRandomSeedProvidedOnCommandLine == 0 ]; then
    # gtest works with a random seed between 0 and 99,999. We approximate this range (generating a random number between 0 and 98301):
    randomSeed=$(($RANDOM*($RANDOM%%3+1)))
  else
    randomSeed=$crimsonRandomSeed
  fi
  shuffleCommand=--gtest_shuffle
  seedCommand=--gtest_random_seed=$randomSeed
fi
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export CRIMSON_FLOWSOLVER_HOME="%s"
export PHASTA_CONFIG="%s/configs"
export LD_LIBRARY_PATH=""
# SET PYTHON PATH HERE IF YOU'RE HAVING PROBLEMS WITH MODULE LOADS IN PYTHON
#export PYTHONPATH=/usr/lib/python2.7/dist-packages

mpirun -np $1 ./test $shuffleCommand $seedCommand $2 $3 $4 $5 $6 $7 $8 $9
"""
