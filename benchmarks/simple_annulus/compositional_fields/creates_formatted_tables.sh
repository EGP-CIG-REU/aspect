#!/bin/bash

interpolation=('bilinear' 'cell_average')
grid_resolution=('2' '3' '4' '5' '6')

top_level_dir=`pwd`

for interpolation_scheme in ${interpolation[@]}; do
  cd $interpolation_scheme
  
  for grid in ${grid_resolution[@]}; do 
    for grid_res in ${grid_resolution[@]}; do
      grep -i "Error" $grid/${grid}x${grid}_$((n_p_c))_output/log.txt >> raw_data_tmp
    done
  done 

  i=0
  echo $interpolation_scheme
  echo "Grid_resolution		L2_velocity	L2_pressure" >> formatted_tables_tmp
  cat raw_data_tmp | while read line; do 
    l2_u=`echo $line | cut -d':' -f2 | cut -d',' -f3`
    l2_p=`echo $line | cut -d':' -f2 | cut -d',' -f4`
   
    echo "${grid_resolution[i]}		$l2_u		$l2_p" >> formatted_tables_tmp
    i=$((i+1))
  done

  cd $top_level_dir
done
