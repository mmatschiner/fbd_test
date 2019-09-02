class Array
	def sum
		sum = 0
		self.each {|x| sum += x}
		sum
	end
	def mean
		self.sum/self.size.to_f
	end
	def median
		if self.size.modulo(2) == 1
			sorted_array = self.sort
			median = sorted_array[self.size/2]
		else
			sorted_array = self.sort
			median = (sorted_array[self.size/2]+sorted_array[(self.size/2)+1])/2.0
		end
		median
	end
	def hpd_lower(proportion)
		raise "The interval should be between 0 and 1!" if proportion >= 1 or proportion <= 0
		sorted_array = self.sort
		hpd_index = 0
		min_range = sorted_array[-1]
		diff = (proportion*self.size).round
		(self.size-diff).times do |i|
			min_value = sorted_array[i]
			max_value = sorted_array[i+diff-1]
			range = max_value - min_value
			if range < min_range
				min_range = range
				hpd_index = i
			end
		end
		sorted_array[hpd_index]
	end
	def hpd_upper(proportion)
		raise "The interval should be between 0 and 1!" if proportion >= 1 or proportion <= 0
		sorted_array = self.sort
		hpd_index = 0
		min_range = sorted_array[-1]
		diff = (proportion*self.size).round
		(self.size-diff).times do |i|
			min_value = sorted_array[i]
			max_value = sorted_array[i+diff-1]
			range = max_value - min_value
			if range < min_range
				min_range = range
				hpd_index = i
			end
		end
		sorted_array[hpd_index+diff-1]
	end
end

# Get the command-line arguments.
analysis_dir = ARGV[0]
estimates_file_name = ARGV[1]
svg_file_name = ARGV[2]

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
		if de.match(/.*\.mrca/)
			replicate_dirs << "#{analysis_dir}/replicates/#{p}"
			break
		end
	end
end
replicate_dirs.sort!

# If no .mrca files were found, raise an issue.
raise "No files with ending .mrca were found in directory #{analysis_dir}" if replicate_dirs == []

# Prepare arrays for the true age, and for mean, and for the lower and upper boundary of the 95% HPD interval.
true_ages = []
estimated_ages_means = []
estimated_ages_lowers = []
estimated_ages_uppers = []

# Prepare the output string.
estimates_output = "true age\testimated age\n"

# Repeat for each replicate directory.
replicate_dirs.each do |r|

	# Read the burnin.txt file for this replicate.
	burnin_file = File.open("#{r}/burnin.txt")
	burnin = burnin_file.read.to_i

	# Read the .mrca log file for this replicate.
	replicate_dir_entries = Dir.entries(r)
	mrca_log_file_name = ""
	replicate_dir_entries.each do |e|
		if e.match(/.*\.mrca/)
			if mrca_log_file_name == ""
				mrca_log_file_name = e
			else
				raise "Found more than one file with ending \".mrca\" in \"#{r}\"!"
			end
		end
	end
	print "\rReading file #{r}/#{mrca_log_file_name}..."
	mrca_log_file = File.open("#{r}/#{mrca_log_file_name}")
	mrca_log_lines = mrca_log_file.readlines
	mrca_log_file.close
	header = mrca_log_lines[0]
	mrca_log_lines_without_header = mrca_log_lines[1..-1]
	mrca_log_lines_without_burnin = []
	mrca_log_lines_without_header.each do |l|
		mrca_log_lines_without_burnin << l if l.split[0].to_i > burnin	
	end

	# Analyze the .mrca log file for this replicate.
	header_ary = header.split("\t")
	mrca_indices = []
	header_ary.size.times do |x|
		if header_ary[x].match("TreeHeight.t:tree") or header_ary[x].match(/mrcatime/)
			mrca_indices << x
		end
	end
	mrca_indices.each do |x|
		true_age = mrca_log_lines_without_header[0].split[x].to_f
		estimated_ages = []
		mrca_log_lines_without_burnin.each do |l|
			mrca_log_lines_without_burnin_ary = l.split("\t")
			estimated_ages << mrca_log_lines_without_burnin_ary[x].to_f
		end

		# Determine mean, lower, and upper HPD interval of estimates.
		estimated_ages_mean = estimated_ages.mean
		estimated_ages_lower = estimated_ages.hpd_lower(0.95)
		estimated_ages_upper = estimated_ages.hpd_upper(0.95)

		# Add to arrays.
		true_ages << true_age
		estimated_ages_means << estimated_ages_mean
		estimated_ages_lowers << estimated_ages_lower
		estimated_ages_uppers << estimated_ages_upper

		# Add true ages and estimates to output.
		estimates_output << "#{'%.5f' % true_age}\t#{'%.5f' % estimated_ages_mean} (#{'%.5f' % estimated_ages_lower}-#{'%.5f' % estimated_ages_upper})\n"
	end

end

# Feedback.
puts " done."

# Write the estimates file.
estimates_out_file = File.new(estimates_file_name, "w")
estimates_out_file.write(estimates_output)
estimates_out_file.close

# Feedback.
print "Analyzing estimates..."

# Some specifications for the SVG output.
dimX = 600
dimY = 600
max_age = 125
cr = 2
line_width = 2
frame_stroke_width = 2
window_size = 20
dot_color = "grey"
dot_alpha = 1.0
hpd_color = "grey"
hpd_alpha = 0.5
hpd_width = 1.0
font_family = "Helvetica"
font_size = 14
text_x = 20

# Initiate the stats output.
stats_output = ""

# Prepare the header of the SVG string.
svg_output = ""
svg_output << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
svg_output << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n"
svg_output << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"#{dimX}\" height=\"#{dimY}\" viewBox=\"0 0 #{dimX} #{dimY}\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n"
svg_output << "\n"

# Write some style definitions.
svg_output << "  <defs>\n"
svg_output << "    <circle id=\"dot\" cx=\"0\" cy=\"0\" r=\"#{cr}\" fill=\"#{dot_color}\" fill-opacity=\"#{dot_alpha}\" />\n"
svg_output << "    <path id=\"hpd\" cx=\"0\" cy=\"0\" stroke=\"#{hpd_color}\" stroke-width=\"#{hpd_width}\" stroke-opacity=\"#{hpd_alpha}\" />\n"
svg_output << "  </defs>\n"
svg_output << "\n"

# Write the frame.
svg_output << "  <!--Frame-->\n"
svg_output << "  <rect style=\"stroke:black; stroke-width:#{frame_stroke_width/2.0}px; fill:none\" x=\"0\" y=\"0\" width=\"#{dimX}\" height=\"#{dimY}\" />\n"
svg_output << "\n"

# Write bars for each hpd.
svg_output << "  <!--Highest Posterior Density intervals-->\n"
true_ages.size.times do |a|
	cx = ((true_ages[a]/max_age.to_f)*dimX)
	cy_lower = dimY-((estimated_ages_lowers[a]/max_age.to_f)*dimY)
	cy_upper = dimY-((estimated_ages_uppers[a]/max_age.to_f)*dimY)
	svg_output << "  <path id=\"hpd\" style=\"stroke:#{hpd_color}; stroke-width:#{hpd_width}px; stroke-opacity:#{hpd_alpha}\" d=\"M #{cx},#{cy_lower} L #{cx},#{cy_upper}\" />\n"
end
svg_output << "\n"

# Write dots for each mean estimate.
svg_output << "  <!--Dots-->\n"
true_ages.size.times do |a|
	cx = ((true_ages[a]/max_age.to_f)*dimX)
	cy = dimY-((estimated_ages_means[a]/max_age.to_f)*dimY)
	svg_output << "  <use x=\"#{cx}\" y=\"#{cy}\" xlink:href=\"#dot\"/>\n"
end
svg_output << "\n"

# Write the optimum line.
svg_output << "  <!--Optimum-->\n"
svg_output << "  <path id=\"optimum\" style=\"fill:none; stroke:black; stroke-width:#{line_width}px\" stroke-dasharray=\"10,10\" d=\"M 0,#{dimY} L #{dimX},0\" />\n"
svg_output << "\n"

# Write the mean line of mean estimates.
svg_output << "  <!--Mean of mean estimates-->\n"
collection_x_means = []
collection_y_means = []
(window_size/2).upto(max_age-(window_size/2)) do |x|
	collection_x = []
	collection_y = []
	true_ages.size.times do |a|
		if true_ages[a] > x - (window_size/2) and true_ages[a] <= x + (window_size/2)
			collection_x << true_ages[a]
			collection_y << estimated_ages_means[a]
		end
	end
	if collection_x.size >= 5
		collection_x_means << collection_x.mean
		collection_y_means << collection_y.mean
	end
end
cxStart = (collection_x_means[0]/max_age.to_f)*dimX
cyStart = dimY - (collection_y_means[0]/max_age.to_f)*dimY
svg_output << "  <path id=\"mean\" style=\"fill:none; stroke:black; stroke-width:#{line_width}px\" d=\"M #{cxStart},#{cyStart} "
collection_x_means.size.times do |x|
	cx = (collection_x_means[x]/max_age.to_f)*dimX
	cy = dimY - (collection_y_means[x]/max_age.to_f)*dimY
	svg_output << "L #{cx},#{cy} "
end
svg_output << "\"/>\n"
svg_output << "\n"

# Prepare the stats part of the SVG.
svg_output << "  <!--Stats-->\n"

# Calculate the RMSD value between true ages and mean estimates.
sd = 0
true_ages.size.times {|a| sd += (true_ages[a]-estimated_ages_means[a])**2}
msd = sd/true_ages.size.to_f
rmsd = Math.sqrt(msd)
rmsd_y = 30
svg_output << "  <text x=\"#{text_x}\" y=\"#{rmsd_y}\" font-family=\"#{font_family}\" font-size=\"#{font_size}\" fill=\"black\">\n"
svg_output << "    RMSD: #{rmsd.round(2)}\n"
svg_output << "  </text>\n"

# Calculate the overall proportion of true ages in bins of 10 time units that within the 95% HPD interval.
inside_hpd_000_020 = 0
inside_hpd_020_040 = 0
inside_hpd_040_060 = 0
inside_hpd_060_080 = 0
inside_hpd_080_100 = 0
inside_hpd_100_120 = 0
true_ages_000_020 = 0
true_ages_020_040 = 0
true_ages_040_060 = 0
true_ages_060_080 = 0
true_ages_080_100 = 0
true_ages_100_120 = 0
true_ages.size.times do |a|
	if true_ages[a] < 20
		true_ages_000_020 += 1
		if estimated_ages_uppers[a] >= true_ages[a] and estimated_ages_lowers[a] <= true_ages[a]
			inside_hpd_000_020 += 1
		end
	elsif true_ages[a] < 40
		true_ages_020_040 += 1
		if estimated_ages_uppers[a] >= true_ages[a] and estimated_ages_lowers[a] <= true_ages[a]
			inside_hpd_020_040 += 1
		end
	elsif true_ages[a] < 60
		true_ages_040_060 += 1
		if estimated_ages_uppers[a] >= true_ages[a] and estimated_ages_lowers[a] <= true_ages[a]
			inside_hpd_040_060 += 1
		end
	elsif true_ages[a] < 80
		true_ages_060_080 += 1
		if estimated_ages_uppers[a] >= true_ages[a] and estimated_ages_lowers[a] <= true_ages[a]
			inside_hpd_060_080 += 1
		end
	elsif true_ages[a] < 100
		true_ages_080_100 += 1
		if estimated_ages_uppers[a] >= true_ages[a] and estimated_ages_lowers[a] <= true_ages[a]
			inside_hpd_080_100 += 1
		end
	elsif true_ages[a] < 120
		true_ages_100_120 += 1
		if estimated_ages_uppers[a] >= true_ages[a] and estimated_ages_lowers[a] <= true_ages[a]
			inside_hpd_100_120 += 1
		end
	end
end
inside_hpd_proportion_000_020 = inside_hpd_000_020/true_ages_000_020.to_f
inside_hpd_proportion_020_040 = inside_hpd_020_040/true_ages_020_040.to_f
inside_hpd_proportion_040_060 = inside_hpd_040_060/true_ages_040_060.to_f
inside_hpd_proportion_060_080 = inside_hpd_060_080/true_ages_060_080.to_f
inside_hpd_proportion_080_100 = inside_hpd_080_100/true_ages_080_100.to_f
inside_hpd_proportion_100_120 = inside_hpd_100_120/true_ages_100_120.to_f
stats_output << "Number of true ages 000-020: #{true_ages_000_020}\n"
stats_output << "Number of true ages 020-040: #{true_ages_020_040}\n"
stats_output << "Number of true ages 040-060: #{true_ages_040_060}\n"
stats_output << "Number of true ages 060-080: #{true_ages_060_080}\n"
stats_output << "Number of true ages 080-100: #{true_ages_080_100}\n"
stats_output << "Number of true ages 100-120: #{true_ages_100_120}\n"

stats_output << "Proportion within HPD interval 000-020: #{inside_hpd_proportion_000_020}\n"
stats_output << "Proportion within HPD interval 020-040: #{inside_hpd_proportion_020_040}\n"
stats_output << "Proportion within HPD interval 040-060: #{inside_hpd_proportion_040_060}\n"
stats_output << "Proportion within HPD interval 060-080: #{inside_hpd_proportion_060_080}\n"
stats_output << "Proportion within HPD interval 080-100: #{inside_hpd_proportion_080_100}\n"
stats_output << "Proportion within HPD interval 100-120: #{inside_hpd_proportion_100_120}\n"

# Calculate the proportion of true ages that is lower than, within, or greater than the 95% HPD interval.
lower_than_hpd = 0
inside_hpd = 0
greater_than_hpd = 0
true_ages.size.times do |a|
	if estimated_ages_uppers[a] < true_ages[a]
		greater_than_hpd += 1
	elsif estimated_ages_uppers[a] >= true_ages[a] and estimated_ages_lowers[a] <= true_ages[a]
		inside_hpd += 1
	elsif estimated_ages_lowers[a] > true_ages[a]
		lower_than_hpd += 1
	else
		raise "Unexpected relationship of true ages and hpd!"
	end
end
lower_than_hpd_proportion = lower_than_hpd/true_ages.size.to_f
inside_hpd_proportion = inside_hpd/true_ages.size.to_f
greater_than_hpd_proportion = greater_than_hpd/true_ages.size.to_f

# Write the overall proportion of ages below, inside, and above the HPD to the svg string.
lower_than_hpd_proportion_y = 70
svg_output << "  <text x=\"#{text_x}\" y=\"#{lower_than_hpd_proportion_y}\" font-family=\"#{font_family}\" font-size=\"#{font_size}\" fill=\"black\">\n"
svg_output << "    true age &#60; 95% HPD: #{(lower_than_hpd_proportion*100).round(1)}%\n"
svg_output << "  </text>\n"
inside_hpd_proportion_y = 90
svg_output << "  <text x=\"#{text_x}\" y=\"#{inside_hpd_proportion_y}\" font-family=\"#{font_family}\" font-size=\"#{font_size}\" fill=\"black\">\n"
svg_output << "    true age in 95% HPD: #{(inside_hpd_proportion*100).round(1)}%\n"
svg_output << "  </text>\n"
greater_than_hpd_proportion_y = 110
svg_output << "  <text x=\"#{text_x}\" y=\"#{greater_than_hpd_proportion_y}\" font-family=\"#{font_family}\" font-size=\"#{font_size}\" fill=\"black\">\n"
svg_output << "    true age &#62; 95% HPD: #{(greater_than_hpd_proportion*100).round(1)}%\n"
svg_output << "  </text>\n"
svg_output << "\n"
# Calculate the overall mean HPD width.
hpd_widths = []
true_ages.size.times do |a|
	hpd_widths << estimated_ages_uppers[a]-estimated_ages_lowers[a]
end
mean_hpd_width = hpd_widths.mean

# Write the mean HPD widths to the svg string.
mean_hpd_width_y = 50
svg_output << "  <text x=\"#{text_x}\" y=\"#{mean_hpd_width_y}\" font-family=\"#{font_family}\" font-size=\"#{font_size}\" fill=\"black\">\n"
svg_output << "    mean HPD width: #{mean_hpd_width.round(2)}\n"
svg_output << "  </text>\n"

# Finalize the SVG string
svg_output << "</svg>\n"

# Write the SVG string to file.
svg_file = File.new(svg_file_name, 'w')
svg_file.write(svg_output)
svg_file.close

# Feedback.
puts " done."
puts "Wrote plot to file #{svg_file_name}."
