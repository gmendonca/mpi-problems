# mpi-problems

In order to run this project you need to install [mpich2](https://www.mpich.org/).
You can run the [makefile](https://github.com/gmendonca/mpi-problems/blob/master/makefile) to compile.

```sh
$ make
```

Or you can compile separately, like:

```sh
$ mpicc -c get_data.c
$ mpicc -o get_data get_data.o
```

```sh
$ mpicc -c get_data_modified.c
$ mpicc -o get_data_modified get_data_modified.o
```

```sh
$ mpicc -c gauss_mpi.c
$ mpicc -o gauss_mpi gauss_mpi.o
```
To run the codes, do the following:

```sh
$ bash run_get_data.bash
```

```sh
$ bash run_get_datamodified.bash
```

```sh
$ bash run_gauss.bash
```

Or run the directly on the terminal, for example, with 8 process:

```sh
mpirun -n 8 ./get_data
```

```sh
mpirun -n 8 ./get_data_modified
```

```sh
mpirun -n 8 ./gauss_mpi
```
