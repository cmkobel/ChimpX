#!/bin/bash

matchcount=100
step=1

for id in {1..500}
do
	echo ${id}
	lastz ../00reference/pan_tro_3.0/windows/window_${id}_*.fasta --self --step=$step --notransition --exact=$matchcount --format=rdotplot > 0dotplots/window_${id}.dotplot
done


	

#lastz ../fasta_files_ampliconic/ampliconic_region${id}.fa ../fasta_files_ampliconic/ampliconic_region${id}.fa --step=$step --notransition --exact=$matchcount --format=rdotplot > testing1.dotplot
	