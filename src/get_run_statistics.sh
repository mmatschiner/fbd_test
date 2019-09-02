# m_matschiner Wed Jul 4 22:30:20 CEST 2018

for analysis_dir in ../res/beast/{CladeAge,FBD}*
do
	analysis_id=`basename ${analysis_dir}`
	mkdir -p ${analysis_dir}/summary
	echo "analysis_id replicate_id iterations burnin time_per_iteration min_ess iterations_required" > ${analysis_dir}/summary/run_statistics.txt
	for replicate_dir in ${analysis_dir}/replicates/r???
	do
		# Get the replicate ID.
		replicate_id=`basename ${replicate_dir}`
		
		# Get the total number of mcmc iterations.
		iterations=`cat ${replicate_dir}/*.log | tail -n 1 | cut -f 1`

		# Get the time per iterations from the respective snapp_out.txt file.
		time_string=`cat ${replicate_dir}/*_out.txt | grep Msamples | tail -n 1 | tr -s " " | cut -d " " -f 5 | sed 's/\/Msamples//g'`
		time_per_iteration=`ruby convert_time_string.rb ${time_string}`

		# Read the burnin from the burnin.txt file.
		burnin_file=${replicate_dir}/burnin.txt
		burnin=`cat ${burnin_file}`

		# Get the number of iterations required for run convergence.
		unset R_HOME
		if [ ${analysis_dir} == ../res/beast/FBD_Range3 ]
		then
			iterations_required="NA"
		else
			iterations_required=`Rscript get_number_of_iterations_to_convergence.r ${replicate_dir}/*.log ${burnin} | tr -d " "`
		fi

		# Get the lowest ess value.
		min_ess=`Rscript get_min_ess.r ${replicate_dir}/*.log ${burnin} | tr -d " "`

		# Print statistics.
		echo "${analysis_id} ${replicate_id} ${iterations} ${burnin} ${time_per_iteration} ${min_ess} ${iterations_required}"
	done >> ${analysis_dir}/summary/run_statistics.txt
done
