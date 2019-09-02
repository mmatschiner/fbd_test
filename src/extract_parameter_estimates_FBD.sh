# m_matschiner Thu Jul 5 17:15:45 CEST 2018

# Extract estimates for the parameters diversification rate, turnover, and sampling proportion.
for analysis_dir in ../res/beast/{CladeAge,FBD}*
do
	estimates_file=${analysis_dir}/summary/parameters.txt
	ruby extract_parameter_estimates_FBD.rb ${analysis_dir} ${estimates_file}
done
