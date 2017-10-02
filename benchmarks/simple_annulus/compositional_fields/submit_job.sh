#!/bin/bash

interpolation=('bilinear' 'cell_average')
grid_resolution=('2' '3' '4' '5' '6' '7')
job_template_filename='job_template.sh'
#density_directory_names=('rho_0' 'rho_1' 'rho_r' 'rho_r^2')
density_directory_names=('rho_r^2')
FEM_types=('Q1_P0' 'Q2_Q1' 'Q2_P-1' 'Q3_Q2' 'Q3_P-2')
top_level_dir=`pwd`

for density_dir in ${density_directory_names[@]}; do

  top_level_dir=`echo $PWD`
  cd $density_dir

  for FEM_type in ${FEM_types[@]}; do
    cd $FEM_type

for interpolation_scheme in ${interpolation[@]}; do
  cd $interpolation_scheme
  for grid_res in ${grid_resolution[@]}; do
    cd $grid_res
    echo $PWD
    sbatch $job_template_filename
    cd ..
  done 
   cd ..
  done
   cd ..
  done
  cd $top_level_dir
done
