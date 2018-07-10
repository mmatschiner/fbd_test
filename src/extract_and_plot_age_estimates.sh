# m_matschiner Thu Jul 5 15:45:28 CEST 2018

# Extract and plot the true node ages versus the estimated node ages.
for analysis_dir in ../res/beast/{CladeAge,FBD}*
do
	estimates_file=${analysis_dir}/summary/estimates.txt
	plot=${analysis_dir}/summary/estimates.svg
	ruby extract_and_plot_age_estimates.rb ${analysis_dir} ${estimates_file} ${plot}
done
