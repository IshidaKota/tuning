#!/bin/bash
#PJM -L "node=36"
#PJM -L "elapse=10:00"
#PJM -j
#PJM -X

module load netcdf/4.4.1.1
module unload intel/2017
module load oneapi/2021.3

NUM_NODES=1
NUM_CORES=36
NUM_PROCS=36

export I_MPI_PERHOST=$NUM_CORES
export I_MPI_FABRICS=shm:ofi

export I_MPI_HYDRA_BOOTSTRAP=rsh
export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh
export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}

mpiexec.hydra -n $NUM_PROCS ./fvcom --casename tst
