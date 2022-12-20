#!/bin/bash
#by ishida 2022/10

module load hdf5 
module load netcdf-fortran 
module load netcdf 
. /work/gy29/y29007/miniconda/etc/profile.d/conda.sh
export HOME=/work/gy29/y29007

export LD_LIBRARY_PATH="${HOME}/Github/fvcom442/FVCOM_source/libs/install/lib:$LD_LIBRARY_PATH"

#イメージ:1~9,10,11~12月で分ける

mpiexec.hydra -n ${PJM_MPI_PROC} ./fvcom --casename Tokyo2
wait
#Tokyo2は通常通りのファイル、ただしRST_ONにする。
mv ./Tokyo2_restart_0001.nc ../input_testcase2
#restartfileの移動
wait
mpiexec.hydra -n ${PJM_MPI_PROC} ./fvcom --casename Tokyo2_hot
wait
#再びrestartfileを生成しておいて、最後に使う
mv ./Tokyo2_hot_restart_0001.nc ../input_testcase2
wait
mpiexec.hydra -n ${PJM_MPI_PROC} ./fvcom --casename Tokyo_hot3
wait
conda activate nco
echo 'begin merging netcdf output files...'
ncrcat -O Tokyo2_0001.nc Tokyo2_hot_0001.nc Tokyo2_hot2_0001.nc Tokyo_all.nc
conda deactivate
echo 'All done!'