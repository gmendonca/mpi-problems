MPICC:=mpicc

program: gauss_mpi.c get_data.c get_data_modified.c

	 $(MPICC) -c gauss_mpi.c

	 $(MPICC) -o gauss_mpi gauss_mpi.o

	 $(MPICC) -c get_data.c

	 $(MPICC) -o get_data get_data.o

	 $(MPICC) -c get_data_modified.c

	 $(MPICC) -o get_data_modified get_data_modified.o

clean:
	rm -rf *.o
