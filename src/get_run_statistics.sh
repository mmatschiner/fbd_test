# m_matschiner Wed Jul 4 22:30:20 CEST 2018

for analysis_dir in ../res/beast/{CladeAge,FBD}*
do
	analysis_id=`basename ${analysis_dir}`
	mkdir -p ${analysis_dir}/summary
	echo "analysis_id replicate_id time_per_iteration min_ess iterations_required" > ${analysis_dir}/summary/run_statistics.txt
	for replicate_dir in ${analysis_dir}/replicates/r???
	do
		# Get the replicate ID.
		replicate_id=`basename ${replicate_dir}`
		
		# Get the time per iterations from the respective snapp_out.txt file.
		time_string=`cat ${replicate_dir}/*_out.txt | grep Msamples | tail -n 1 | tr -s " " | cut -d " " -f 5 | sed 's/\/Msamples//g'`
		time_per_iteration=`ruby convert_time_string.rb ${time_string}`

		# Get the lowest ess value.
		unset R_HOME
		min_ess=`Rscript get_min_ess.r ${replicate_dir}/*.log`

		# Get the number of iterations required for run convergence.
		iterations_required=`Rscript get_number_of_iterations_to_convergence.r ${replicate_dir}/*.log`

		# Print statistics.
		echo "${analysis_id} ${replicate_id} ${time_per_iteration} ${min_ess} ${iterations_required}"

	done >> ${analysis_dir}/summary/run_statistics.txt
done
