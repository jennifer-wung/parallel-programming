# compiling and linking all examples.
all: MS_Hybrid_dynamic MS_Hybrid_static MS_MPI_dynamic MS_MPI_static MS_OpenMP_dynamic MS_OpenMP_static

# compiling 
MS_Hybrid_dynamic:	MS_Hybrid_dynamic.cpp 
	mpic++ -fopenmp -lX11 MS_Hybrid_dynamic.cpp -o MS_Hybrid_dynamic

MS_Hybrid_static:	MS_Hybrid_static.cpp
	mpic++ -fopenmp -lX11 -o MS_Hybrid_static MS_Hybrid_static.cpp 

MS_MPI_dynamic:		MS_MPI_dynamic.cpp
	mpic++ -o MS_MPI_dynamic -lX11 MS_MPI_dynamic.cpp

MS_MPI_static:	MS_MPI_static
	mpic++ -o MS_MPI_static -lX11 MS_MPI_static.cpp

MS_OpenMP_dynamic:	MS_OpenMP_dynamic.cpp
	g++ -fopenmp -lX11 -o MS_OpenMP_dynamic MS_OpenMP_dynamic.cpp

MS_OpenMP_static:	MS_OpenMP_static.cpp
	g++ -fopenmp -lX11 -o MS_OpenMP_static MS_OpenMP_static.cpp
