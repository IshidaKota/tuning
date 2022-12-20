#!/bin/bash
#by ishid 2022/12

. /work/gy29/y29007/miniconda/etc/profile.d/conda.sh
conda activate pyfvcom

python tune.py Tokyo5 temperature
python tune.py Tokyo5 salinity

python plot_contour.py Tokyo5 temperature
python plot_contour.py Tokyo5 salinity

python plot_vertical_contour.py Tokyo5 temperature
python plot_vertical_contour.py Tokyo5 salinity

cd doc
platex pic5
dvipdfmx pic5
mv pic5.pdf pic5_w20r10.pdf
cd ../../run
mv Tokyo5_0001.nc Tokyo5w20r05.nc
cd ../tune
