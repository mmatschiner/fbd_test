# m_matschiner Thu Jul 5 17:15:38 CEST 2018

# Get the command-line arguments.
analysis_dir = ARGV[0]
estimates_file_name = ARGV[1]
# See whether the analysis directory exists.
raise "The analysis directory #{analysis_dir} could not be found!" unless Dir.exists?(analysis_dir)

# Determine the replicate directories and see whether they contain .log log files
potential_replicate_dirs = []
Dir.entries("#{analysis_dir}/replicates").each do |e|
	potential_replicate_dirs << e if e.match(/r\d\d\d/)
end
replicate_dirs = []
potential_replicate_dirs.each do |p|
	dir_entries = Dir.entries("#{analysis_dir}/replicates/#{p}")
	dir_entries.each do |de|
		if de.match(/.*\.log/)
			replicate_dirs << "#{analysis_dir}/replicates/#{p}"
			break
		end
	end
end
replicate_dirs.sort!

# If no .log files were found, raise an issue.
raise "No files with ending .log were found in directory #{analysis_dir}" if replicate_dirs == []

# Initiate arrays for estimates of the diversification rate, the turnover, and teh sampling proportion.
diversification_rates = [] # diversificationRateFBD.t:tree1
turnovers = [] # turnoverFBD.t:tree1
sampling_proportions = [] # samplingProportionFBD.t:tree1

# Prepare the output string.
estimates_output = "diversification_rate\tturnover\tsampling_proportion\n"

# Repeat for each replicate directory.
replicate_dirs.each do |r|

	# Initiate arrays for this replicate, for estimates of the diversification rate, the turnover, and teh sampling proportion.
	diversification_rates_this_rep = [] # diversificationRateFBD.t:tree1
	turnovers_this_rep = [] # turnoverFBD.t:tree1
	sampling_proportions_this_rep = [] # samplingProportionFBD.t:tree1

	# Read the burnin.txt file for this replicate.
	burnin_file = File.open("#{r}/burnin.txt")
	burnin = burnin_file.read.to_i

	# Read the .log file for this replicate.
	replicate_dir_entries = Dir.entries(r)
	log_file_name = ""
	replicate_dir_entries.each do |e|
		if e.match(/.*\.log/)
			if log_file_name == ""
				log_file_name = e
			else
				raise "Found more than one file with ending \".log\" in \"#{r}\"!"
			end
		end
	end
	print "\rReading file #{r}/#{log_file_name}..."
	log_file = File.open("#{r}/#{log_file_name}")
	log_lines = log_file.readlines
	log_file.close
	x = 0
	log_line = log_lines[x]
	while log_line[0] == "#"
		x += 1
		log_line = log_lines[x]
	end
	header = log_lines[x]
	log_lines_without_header = log_lines[x+1..-1]
	log_lines_without_burnin = []
	log_lines_without_header.each do |l|
		log_lines_without_burnin << l if l.split[0].to_i > burnin	
	end

	# Analyze the .log file for this replicate.
	header_ary = header.split("\t")
	diversification_rate_index = 0
	turnover_index = 0
	sampling_proportion_index = 0
	header_ary.size.times do |x|
		if header_ary[x].match(/diversificationRateFBD/)
			diversification_rate_index = x
		elsif header_ary[x].match(/turnoverFBD/)
			turnover_index = x
		elsif header_ary[x].match(/samplingProportionFBD/)
			sampling_proportion_index = x
		end
	end
	log_lines_without_burnin.each do |l|
		line_ary = l.split("\t")
		diversification_rates_this_rep << line_ary[diversification_rate_index].to_f
		turnovers_this_rep << line_ary[turnover_index].to_f
		sampling_proportions_this_rep << line_ary[sampling_proportion_index].to_f
	end

	# Reduce the number of samples to 1000.
	diversification_rates_this_rep_red = diversification_rates_this_rep.sample(1000)
	turnovers_this_rep_red = turnovers_this_rep.sample(1000)
	sampling_proportions_this_rep_red = sampling_proportions_this_rep.sample(1000)

	# Add posterior samples of this replicate to the arrays for all replicates.
	diversification_rates_this_rep_red.each {|i| diversification_rates << i}
	turnovers_this_rep_red.each {|i| turnovers << i}
	sampling_proportions_this_rep_red.each {|i| sampling_proportions << i}
end

# Add parameter estimates to output.
diversification_rates.size.times do |x|
	estimates_output << "#{'%.5f' % diversification_rates[x]}\t#{'%.5f' % turnovers[x]}\t#{'%.5f' % sampling_proportions[x]}\n"
end

# Feedback.
puts " done."

# Write the estimates file.
estimates_file = File.new(estimates_file_name, "w")
estimates_file.write(estimates_output)
estimates_file.close
