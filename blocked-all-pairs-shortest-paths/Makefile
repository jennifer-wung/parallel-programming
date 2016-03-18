all: shortest_paths_cuda.exe shortest_paths_openmp.exe shortest_paths_mpi.exe
shortest_paths_cuda.exe: shortest_paths_cuda.cu
	nvcc -arch=sm_20 -o shortest_paths_cuda.exe shortest_paths_cuda.cu
shortest_paths_openmp.exe: shortest_paths_openmp.cu
	nvcc -arch=sm_20 -o shortest_paths_openmp.exe -Xcompiler="-fopenmp" ./shortest_paths_openmp.cu
shortest_paths_mpi.exe: shortest_paths_mpi.cu
	nvcc -arch=sm_20 -o shortest_paths_mpi.exe -I/usr/lib64/mpich/include -L/usr/lib64/mpich/lib -lmpich ./shortest_paths_mpi.cu
clean:
	rm -f shortest_paths_cuda.exe shortest_paths_openmp.exe shortest_paths_mpi.exe

