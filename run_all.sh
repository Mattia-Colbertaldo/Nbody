./hw2-1/build/serial -n 1000 -s 1 -o serial.txt
./hw2-1/build/openmp -n 1000 -s 1 -t 8 -o omp.txt
mpirun -n 8 ./hw2-2/build/mpi -n 1000 -s 1 -o mpi.txt