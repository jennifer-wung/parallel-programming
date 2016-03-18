#PBS -q multi
#PBS -N Jen_test
#PBS -l nodes=3:ppn=1
#PBS -l pmem=80gb
#PBS -l walltime=0:30:00
#PBS -j oe

#cd $PBS_O_WORKDIR

#mpirun HW3/MS_MPI_static2 4 -2 2 -2 2 2000 2000 disable
#mpirun HW3/MS_Hybrid_dynamic2 10 -2 2 -2 2 2000 2000 disable
mpirun HW3/MS_MPI_dynamic 1 -2 2 -2 2 2000 2000 disable
#HW3/MS_OpenMP_dynamic 12 -2 2 -2 2 2000 2000 disable
