#!/bin/bash
#by ishida 2022/10


echo "enter the number of nodes!"
read node


echo "enter the number of processes!"
read proc


#begin the interactive job
pjsub --interact -g gy29 -L rscgrp=interactive,node=$node --mpi proc=$proc


#mpiexec.hydra -n ${PJM_MPI_PROC} ./fvcom --casename Tokyo