# m_matschiner Thu Jul 5 17:15:38 CEST 2018

# Get the command-line arguments.
analysis_dir = ARGV[0]
estimates_file_name = ARGV[1]
# See whether the analysis directory exists.
raise "The analysis directory #{analysis_dir} could not be found!" unless Dir.exists?(analysis_dir)

# Determine the replicate directories and see whether they contain .mrca log files
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

# Set the burnin proportion.
burnin_proportion = 0.2

# Prepare the output string.
estimates_output = "diversification_rate\tturnover\tsampling_proportion\n"

# Repeat for each replicate directory.
replicate_dirs.each do |r|

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
	burnin = (burnin_proportion * log_lines_without_header.size).round
	log_lines_without_burnin = log_lines_without_header[burnin+1..-1]

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
		diversification_rates << line_ary[diversification_rate_index].to_f
		turnovers << line_ary[turnover_index].to_f
		sampling_proportions << line_ary[sampling_proportion_index].to_f
	end
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
