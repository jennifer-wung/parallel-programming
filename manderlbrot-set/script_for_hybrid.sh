#PBS -q multi
#PBS -N TEST
#PBS -l nodes=5:ppn=2
#PBS -j oe

# -- changing to working directory --
cd $PBS_O_WORKDIR

echo "Running on $PBS_NUM_PPN cores."
# number of OpenMP threads must be matched to the number of cores
# -- openmp environment variables -- 
OMP_NUM_THREADS=$PBS_NUM_PPN;
export OMP_NUM_THREADS

# -- mpi environment variables --
export MV2_ENABLE_AFFINITY=0

echo "Running on $PBS_NUM_NODES nodes."
# -- setting num mpi process --
NUM_MPI_PROCESS_PER_NODE=1;


mpirun -ppn $NUM_MPI_PROCESS_PER_NODE ./MS_Hybrid_static2 $OMP_NUM_THREADS -2 2 -2 2 100 100 disable
#./MS_OpenMP_dynamic 12 -2 2 -2 2 2000 2000 disable
