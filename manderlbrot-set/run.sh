echo "***MPI static***"
mpirun -np 5 ./MS_MPI_static 4 -2 2 -2 2 400 400 enable

echo "***MPI dynamic***"
mpirun -np 5 ./MS_MPI_dynamic 4 -2 2 -2 2 400 400 enable

echo "***OpenMP dynamic***"
./MS_OpenMP_dynamic 4 -2 2 -2 2 400 400 enable

echo "***OpenMP static***"
./MS_OpenMP_static 6 -2 2 -2 2 400 400 enable

echo "***Hybrid static***"
mpirun -np 5 ./MS_Hybrid_static 4 -2 2 -2 2 400 400 enable

echo "***Hybrid dynamic***"
mpirun -np 5 ./MS_Hybrid_dynamic 4 -2 2 -2 2 400 400 enable

