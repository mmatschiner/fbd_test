# m_matschiner Wed Jun 20 12:52:04 CEST 2018

# Load the phylsim package.
$libPath = "./phylsim/"
require "./phylsim/main.rb"

class Array
	def sum
		sum = 0
		self.each {|x| sum += x}
		sum
	end
	def mean
		self.sum/self.size.to_f
	end
end

class Branch2
	attr_accessor :id
	attr_accessor :origin
	attr_accessor :termination
	attr_accessor :end_cause
	attr_accessor :extant
	attr_accessor :parent_id
	attr_accessor :sister_id
	attr_accessor :daughter_id1
	attr_accessor :daughter_id2
	attr_accessor :species_id
	attr_accessor :progeny_id
	attr_accessor :first_occurrence_age
	attr_accessor :first_occurrence_age_in_progeny
	attr_accessor :fulfilled_criteria
	attr_reader :fossil
	def initialize(id)
		@id = id
		@origin = nil
		@termination = nil
		@end_cause = nil
		@extant = nil
		@parent_id = nil
		@sister_id = nil
		@daughter_id1 = nil
		@daughter_id2 = nil
		@species_id = nil
		@progeny_id = []
		@first_occurrence_age = nil
		@first_occurrence_age_in_progeny = nil
		@fulfilled_criteria = []
	end
	def to_s
		output = ""
		output << "ID:                               #{@id}\n"
		output << "Origin:                           #{@origin.round(3)}\n"
		output << "Termination:                      #{@termination.round(3)}\n"
		output << "End cause:                        #{@end_cause}\n"
		output << "Extant:                           #{@extant}\n"
		output << "Parent ID:                        #{@parent_id}\n"
		output << "Sister ID:                        #{@sister_id}\n"
		output << "Daughter ID 1:                    #{@daughter_id1}\n"
		output << "Daughter ID 2:                    #{@daughter_id2}\n"
		output << "Species ID:                       #{@species_id}\n"
		output << "Progeny ID:                       "
		if @progeny_id == []
			output << "none\n"
		else
		    @progeny_id.each {|p| output << "#{p}, "}
		    output.chomp!(", ")
		    output << "\n"
		  end
		if @first_occurrence_age == nil
			output << "First occurrence age:             none\n"
		else
			output << "First occurrence age:             #{@first_occurrence_age}\n"
		end
		if @first_occurrence_age_in_progeny == nil
			output << "First occurrence age in progeny:  none\n"
		else
			output << "First occurrence age in progeny:  #{@first_occurrence_age_in_progeny}\n"
		end
		output << "Fulfills criteria:                "
		@fulfilled_criteria.each do |c|
			output << "#{c}, "
		end
		output.chomp!(", ")
		output << "\n"

		output << "\n"
		output
	end
	def addFossils(fossil)
	  @fossil = fossil # an array with flexible number of items
	end
end

# If no command line arguments are given, write the help text.
if ARGV == [] or ARGV[0] == "-h"
	puts "Available options:"
	puts "  -s scheme name (required)"
	puts "E.g. write 'ruby prepare_BEAST_input.rb -s CladeAge'"
	exit
end

# Define true parameter values and ranges used for their estimation.
true_speciation_rate = 0.12
true_extinction_rate = 0.06
true_net_diversification_rate = true_speciation_rate - true_extinction_rate
true_turnover = true_extinction_rate / true_speciation_rate
true_sampling_rate = 0.01
true_sampling_proportion = true_sampling_rate / (true_sampling_rate + true_extinction_rate)

# Read the scheme id argument.
accepted_scheme_ids = ["CladeAge", "CladeAge_Range", "FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2", "FBD_Range3"]
# CladeAge - place CladeAge priors on all those branches
# FBD - uses the fossilized birth death model with the sampling proportion fixed to the true sampling proportion.
# FBX - as FBD, but estimating the sampling proportion.
# FBY - as FBD, but estimating the diversification rate.
# FBZ - as FBD, but estimating the sampling proportion and the diversification rate.
# CladeAge_Range - as CladeAge but with ranges of uncertainty for all parameter estimates.
# FBD_Range1 - only oldest fossils, estimate diversification rate and turnover from ranges, fix sampling proportion.
# FBD_Range2 - only oldest fossils, estimate sampling proportion, diversification rate, and turnover from ranges.
# FBD_Range3 - all fossils, estimate sampling proportion, diversification rate, and turnover from ranges.
# All schemes ending in RS - as without but trees were sampled with the random sampling scheme instead of the diversified sampling scheme.
if ARGV.include?("-s")
	full_scheme_id = ARGV[ARGV.index("-s")+1]
	scheme_id = full_scheme_id
	random_sampling = false
	if full_scheme_id[-2..-1] == "RS"
		scheme_id = full_scheme_id[0..-3]
		random_sampling = true
	end
	unless accepted_scheme_ids.include?(scheme_id)
		raise_string = "ERROR: The specified calibration scheme id #{scheme_id} is not yet accepted. Please modify the code or use a different scheme (accepted schemes:"
		accepted_scheme_ids.each {|s| raise_string << " #{s},"}
		raise_string.chomp(",")
		raise_string << ")!"
		puts raise_string
		exit 1
	end
	analysis_id = scheme_id
	analysis_dir_name = "../res/beast/#{full_scheme_id}/replicates"
else
	raise "Please specify a scheme name with option \"-s\"!"
end

# Set the number of replicate species trees.
number_of_replicates = 20

# Based on the above arguments, specify directories and ids.
dataset_dir_name = "../res/datasets"

# Repeat for 10 replicates.
number_of_replicates.times do |r|

	# Define the working directory.
	replicate_id = "r#{(r+1).to_s.rjust(3).gsub(" ","0")}"

	# Read the info file and get the number of extant and sampled species.
	if random_sampling
		info_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.random_sampling.reconstructed.seqs.info"
	else
		info_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.reconstructed.seqs.info"
	end
	info_file = File.open(info_file_name)
	info_file_lines = info_file.readlines
	n_extant_species = 0
	n_sampled_species = 0
	info_file_lines.each do |l|
		n_extant_species = l.split.last.to_i if l.strip[0..23] == "Number of extant species"
		n_sampled_species = l.split.last.to_i if l.strip[0..24] == "Number of sampled species"
	end
	if n_extant_species == 0 or n_sampled_species == 0
		puts "ERROR: The number of extant or samples species could not be read from the info file!"
		exit 1
	end

	# Read the branch file and store branch information in array 'branch'
	if random_sampling
		branch_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.random_sampling.reconstructed.seqs.branches.txt"
	else
		branch_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.reconstructed.seqs.branches.txt"
	end
	branch_file = File.open(branch_file_name)
	branch_file_lines = branch_file.readlines
	branch_file.close
	branch = []
	branch_file_lines[3..-1].each do |l|
		l.strip!
		unless l == ""
			line_ary = l.split(":")
			identifier = line_ary[0].strip
			value = line_ary[1].strip
			if identifier == "ID"
				branch << Branch2.new(id = value)
			elsif identifier == "Origin"
				branch.last.origin = value.to_f
			elsif identifier == "Termination"
				branch.last.termination = value.to_f
			elsif identifier == "End cause"
				branch.last.end_cause = value
			elsif identifier == "Extant"
				if value == "true"
					branch.last.extant = true
				elsif value == "false"
					branch.last.extant = false
				end
			elsif identifier == "Parent ID"
				branch.last.parent_id = value
			elsif identifier == "Daughter ID 1"
				branch.last.daughter_id1 = value
			elsif identifier == "Daughter ID 2"
				branch.last.daughter_id2 = value
			elsif identifier == "Species ID"
				branch.last.species_id = value
			elsif identifier == "First occurrence age"
				branch.last.first_occurrence_age = value.to_f
			end
		end
	end

	# Read the dmp file (to add complete fossil information to branches, which is needed for the FBD_ALL and FBX_all schemes)
	if random_sampling
		dump_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.random_sampling.reconstructed.seqs.dmp"
	else
		dump_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.reconstructed.seqs.dmp"
	end
	dmp_tree = Tree.load(fileName = dump_file_name, verbose=false)

	# For each branch, find the correct branch in the dumped tree.
	branch.each do |b|
		branch_candidates_in_dmp_tree = []
		dmp_tree.branch.each do |dtb|
			if b.id == dtb.id
				branch_candidates_in_dmp_tree << dtb
			end
		end
		if branch_candidates_in_dmp_tree.size == 1
			b.addFossils(branch_candidates_in_dmp_tree[0].fossil)
		else
			raise "Branches from files #{branch_file_name} and #{dump_file_name} could not be uniquely assigned to each other."
		end
	end

	# Now add the complete progeny and the sister of all branches.
	branch.each do |b|
		# Collect the complete progeny of this branch.
		progeny_id = []
		progeny_id_this_complete = []
		unless b.daughter_id1 == "none"
			# Add the two daughters to the progeny.
			progeny_id << b.daughter_id1
			progeny_id_this_complete << false
			progeny_id << b.daughter_id2
			progeny_id_this_complete << false
			# Add all the rest until the progeny is complete.
			progeny_id_all_complete = false
			while progeny_id_all_complete == false
				progeny_id.size.times do |p|
					unless progeny_id_this_complete[p]
						branch.each do |bb|
							if bb.id == progeny_id[p]
								progeny_id_this_complete[p] = true
								unless bb.daughter_id1 == "none"
									progeny_id << bb.daughter_id1
									progeny_id_this_complete << false
									progeny_id << bb.daughter_id2
									progeny_id_this_complete << false
								end
							end
						end
					end
				end
				# Test whether all elements of 'progeny_id_this_complete' are true, and if not, keep 'progeny_id_all_complete' as false.
				progeny_id_all_complete = true
				progeny_id_this_complete.each {|x| progeny_id_all_complete = false unless x}
			end # while progeny_id_all_complete == false
		end # unless b.daughter_id1 == "none"
		b.progeny_id = progeny_id

		# Add the first occurrence of the progeny.
		first_occurrence_age_in_progeny = nil
		branch.each do |bb|
			if b.progeny_id.include?(bb.id) and bb.first_occurrence_age != nil
				if first_occurrence_age_in_progeny == nil or bb.first_occurrence_age > first_occurrence_age_in_progeny
					first_occurrence_age_in_progeny = bb.first_occurrence_age
				end
			end
		end			
		b.first_occurrence_age_in_progeny = first_occurrence_age_in_progeny

		# Add the sister of this branch.
		branch.each do |bb|
			if b.parent_id == "treeOrigin"
				if b.id == "b0"
					b.sister_id = "b1"
				elsif b.id == "b1"
					b.sister_id = "b0"
				end
			elsif b.parent_id == bb.id
				if bb.daughter_id1 == b.id
					b.sister_id = bb.daughter_id2
				elsif bb.daughter_id2 == b.id
					b.sister_id = bb.daughter_id1
				end
			end
		end # branch.each do |bb|
	end

	# Analyse all branches to find the number of usable fossils.
	# Criterion 1: Only branches that actually have fossils.
	branch.each {|b| b.fulfilled_criteria << 1 unless b.first_occurrence_age == nil}

	# Criterion 2: Only branches that have fossils that are older than all fossils in the progeny.
	# Scheme C in manuscript. Also used with the FBD process.
	branch.each do |b|
		if b.fulfilled_criteria.include?(1)
			b.fulfilled_criteria << 2 if b.first_occurrence_age_in_progeny == nil or b.first_occurrence_age_in_progeny < b.first_occurrence_age
		end
	end

	# Criterion 11: Among branches that fulfill criterion 2: If multiple branches can be constrained with one and the same fossil, pick all of them.
	# --> Scheme A in manuscript.
	branch.each do |b|
		if b.fulfilled_criteria.include?(2)
			candidates = [b]
			fo = b.first_occurrence_age
			search_on = true
			current_branch = b
			while search_on
				current_id = current_branch.id
				older_id = current_branch.parent_id
				if older_id == "treeOrigin"
					search_on = false
				else
					# Get the parent
					older_branch = nil
					branch.each do |bb|
						if bb.id == older_id
							older_branch = bb
							break
						end
					end
					raise "No parent found of branch #{current_id}!" if older_branch == nil
					if older_branch.first_occurrence_age != nil and older_branch.first_occurrence_age > fo
						search_on = false
					elsif older_branch.first_occurrence_age_in_progeny != nil and older_branch.first_occurrence_age_in_progeny > fo
						search_on = false
					else
						candidates << older_branch
					end
				end
				current_branch = older_branch
			end
			candidates.each {|c| c.fulfilled_criteria << 11}
		end
	end

	# Now, all branches have been characterized.
	# Read the sequence alignment file.
	if random_sampling
		alignment_file_name = "species.fossils.random_sampling.reconstructed.seqs.phy"
	else
		alignment_file_name = "species.fossils.reconstructed.seqs.phy"
	end
	alignment_file = File.open("#{dataset_dir_name}/#{replicate_id}/#{alignment_file_name}")
	alignment_file_lines = alignment_file.readlines
	ntax = alignment_file_lines[0].strip.split(" ")[0].to_i
	if ntax != n_sampled_species
		puts "ERROR: The number of sequences in file #{dataset_dir_name}/#{replicate_id}/#{alignment_file_name} does not match the number of sampled species (#{n_sampled_species}) in the info file!"
		exit 1
	end
	nchar = alignment_file_lines[0].strip.split(" ")[1].to_i
	ids = []
	seqs = []
	alignment_file_lines[1..-1].each do |l|
		alignment_file_ary = l.strip.split(" ")
		ids << alignment_file_ary[0]
		seqs << alignment_file_ary[1]
	end
	unless ids.size == ntax
		puts "ERROR: The number of ids in file #{dataset_dir_name}/#{replicate_id}/#{alignment_file_name} does not match the specified number of taxa!"
		exit 1
	end
	unless seqs.size == ntax
		puts "ERROR: The number of sequences in file #{dataset_dir_name}/#{replicate_id}/#{alignment_file_name} does not match the specified number of taxa!"
		exit 1
	end
	unless seqs[0].length == nchar
		puts "ERROR: The length of sequences in file #{dataset_dir_name}/#{replicate_id}/#{alignment_file_name} does not match the specified number of characters!"
		exit 1
	end
	
	# Read the starting tree.
	if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
		if random_sampling
			starting_tree_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.random_sampling.reconstructed.tre"
		else
			starting_tree_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.reconstructed.tre"
		end
	else
		if scheme_id == "FBD_Range3"
			if random_sampling
				starting_tree_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.random_sampling.reconstructed.fbd_all.tre"
			else
				starting_tree_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.reconstructed.fbd_all.tre"
			end
		else
			if random_sampling
				starting_tree_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.random_sampling.reconstructed.fbd.tre"
			else
				starting_tree_file_name = "#{dataset_dir_name}/#{replicate_id}/species.fossils.reconstructed.fbd.tre"
			end
		end
	end
	starting_tree_file = File.open(starting_tree_file_name)
	starting_tree_lines = starting_tree_file.readlines
	starting_tree_lines[2].match(/(\(.+\))/)
	starting_tree_newick = $1
	starting_tree_newick.gsub!(/\[.*?\]/,"")

	# Prepare the xml string.
	xml = ""
	xml << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
	xml << "<!DOCTYPE SIMULATED [\n"
	xml << "<!ENTITY partitions \"#{analysis_id}\">\n"
	xml << "\<!ENTITY splitpartitions \"#{analysis_id}_1,#{analysis_id}_2,#{analysis_id}_3\">\n"
	xml << "]>\n"
	xml << "<beast beautitemplate=\'Standard\' beautistatus=\'\' namespace=\"beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood\" version=\"2.0\">\n"
	xml << "\n"
	xml << "    <data id=\"#{analysis_id}\" name=\"alignment\">\n"
	if scheme_id == "FBD_Range3"
		ids.size.times {|x| xml << "        <sequence id=\"seq_" + "#{ids[x]}\"".ljust(12) + "taxon=\"" + "#{ids[x]}\"".ljust(12) + "totalcount=\"4\" value=\"#{seqs[x]}\"/>\n"}
	else
		ids.size.times {|x| xml << "        <sequence id=\"seq_" + "#{ids[x]}\"".ljust(8) + "taxon=\"" + "#{ids[x]}\"".ljust(8) + "totalcount=\"4\" value=\"#{seqs[x]}\"/>\n"}
	end
	if ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
		branch.each do |b|
			if b.fulfilled_criteria.include?(2)
				xml << "        <sequence id=\"seq_" + "#{b.id.sub("b","f")}\"".ljust(8) + "taxon=\"" + "#{b.id.sub("b","f")}\"".ljust(8) + "totalcount=\"4\" value=\""
				seqs[0].size.times {xml << "-"}
				xml << "\"/>\n"
			end
		end
	elsif scheme_id == "FBD_Range3"
		branch.each do |b|
			if b.fulfilled_criteria.include?(1)
				b.fossil.size.times do |x|
					xml << "        <sequence id=\"seq_" + "#{b.id.sub("b","f")}_#{x}\"".ljust(12) + "taxon=\"" + "#{b.id.sub("b","f")}_#{x}\"".ljust(12) + "totalcount=\"4\" value=\""
					seqs[0].size.times {xml << "-"}
					xml << "\"/>\n"
				end
			end
		end
	end
	xml << "    </data>\n"
	xml << "\n"
	xml << "    <map name=\"Beta\">beast.math.distributions.Beta</map>\n"
	xml << "    <map name=\"Exponential\">beast.math.distributions.Exponential</map>\n"
	xml << "    <map name=\"InverseGamma\">beast.math.distributions.InverseGamma</map>\n"
	xml << "    <map name=\"LogNormal\">beast.math.distributions.LogNormalDistributionModel</map>\n"
	xml << "    <map name=\"Gamma\">beast.math.distributions.Gamma</map>\n"
	xml << "    <map name=\"Uniform\">beast.math.distributions.Uniform</map>\n"
	xml << "    <map name=\"prior\">beast.math.distributions.Prior</map>\n"
	xml << "    <map name=\"taxon\">beast.evolution.alignment.Taxon</map>\n"
	xml << "    <map name=\"LaplaceDistribution\">beast.math.distributions.LaplaceDistribution</map>\n"
	xml << "    <map name=\"OneOnX\">beast.math.distributions.OneOnX</map>\n"
	xml << "    <map name=\"Normal\">beast.math.distributions.Normal</map>\n"
	xml << "\n"
	xml << "    <taxonset id=\'set_all\' spec=\"TaxonSet\">\n"
	if scheme_id == "FBD_Range3"
		ids.size.times {|x| xml << "        <taxon id=\"" + "#{ids[x]}\"".ljust(12) + "/>\n"}
	else
		ids.size.times {|x| xml << "        <taxon id=\"" + "#{ids[x]}\"".ljust(7) + "/>\n"}
	end
	if ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
		branch.each do |b|
			if b.fulfilled_criteria.include?(2)
				xml << "        <taxon id=\"" + "#{b.id.sub("b","f")}\"".ljust(7) + "/>\n"
			end
		end
	elsif scheme_id == "FBD_Range3"
		branch.each do |b|
			if b.fulfilled_criteria.include?(1)
				b.fossil.size.times do |x|
					xml << "        <taxon id=\"" + "#{b.id.sub("b","f")}_#{x}\"".ljust(12) + "/>\n"
				end
			end
		end
	end
	xml << "    </taxonset>\n"
	xml << "\n"
	xml << "    <LogNormal S=\"@ucldStdev.c:clock\" id=\"LogNormalDistributionModel.c:clock\" meanInRealSpace=\"true\">\n"
	xml << "        <parameter estimate=\"false\" lower=\"0.0\" name=\"M\" upper=\"1.0\" value=\"1.0\"/>\n"
	xml << "    </LogNormal>\n"
	xml << "\n"
	xml << "    <plate var=\'n\' range=\'&partitions;\'>\n"
	xml << "        <plate var=\'k\' range=\'1,2,3\'>\n"
	xml << "            <data data=\"@$(n)\" filter=\"$(k)::3\" id=\"$(n)_$(k)\" spec=\"FilteredAlignment\"/>\n"
	xml << "        </plate>\n"
	xml << "    </plate>\n"
	xml << "\n"
	xml << "    <Gamma id=\"Gamma.0\" name=\"distr\">\n"
	xml << "        <parameter id=\"RealParameter.0\" lower=\"0.0\" name=\"alpha\" upper=\"0.0\" value=\"0.2\"/>\n"
	xml << "        <parameter id=\"RealParameter.01\" lower=\"0.0\" name=\"beta\" upper=\"0.0\" value=\"5.0\"/>\n"
	xml << "    </Gamma>\n"
	xml << "\n"
	xml << "    <OneOnX id=\"OneOnX.0\" name=\"distr\"/>\n"
	xml << "\n"
	xml << "    <Exponential id=\"Exponential.0\" name=\"distr\">\n"
	xml << "        <parameter id=\"RealParameter.02\" lower=\"0.0\" name=\"mean\" upper=\"0.0\" value=\"1\"/>\n"
	xml << "    </Exponential>\n"
	xml << "\n"
	xml << "    <run chainLength=\"100000000\" id=\"mcmc\" spec=\"MCMC\" sampleFromPrior=\'false\' storeEvery=\"100000\">\n"
	xml << "\n"
	xml << "        <state id=\"state\" storeEvery=\"10000\">\n"
	if ["FBX", "FBZ", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "            <parameter id=\"samplingProportionFBD.t:tree1\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">#{true_sampling_proportion}</parameter>\n"
	end
	if ["FBY", "FBZ", "FBD_Range1", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "            <parameter id=\"diversificationRateFBD.t:tree1\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">#{true_net_diversification_rate}</parameter>\n"
	end
	if ["FBD_Range1", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "            <parameter id=\"turnoverFBD.t:tree1\" lower=\"0.0\" name=\"stateNode\" upper=\"1.0\">#{true_turnover}</parameter>\n"
	end
	xml << "            <parameter id=\"ucldStdev.c:clock\" lower=\"0.0\" name=\"stateNode\" upper=\"1\" value=\"0.0001\"/>\n"
	xml << "            <stateNode dimension=\"10\" id=\"rateCategories.c:clock\" spec=\"parameter.IntegerParameter\" value=\"1\"/>\n"
	xml << "            <plate var=\'n\' range=\'&splitpartitions;\'>\n"
	xml << "                <parameter id=\"clockRate.c:$(n)\" name=\"stateNode\" value=\"1.0\"/>\n"
	xml << "                <stateNode id=\"RBcount.s:$(n)\" lower=\"0\" spec=\"parameter.IntegerParameter\" upper=\"5\" value=\"5\"/>\n"
	xml << "                <parameter dimension=\"5\" id=\"RBrates.s:$(n)\" lower=\"0.01\" name=\"stateNode\" upper=\"100.0\" value=\"1\"/>\n"
	xml << "                <parameter id=\"gammaShape.s:$(n)\" name=\"stateNode\" value=\"1.0\"/>\n"
	xml << "            </plate>\n"
	xml << "            <tree id=\"Tree.t:tree1\" name=\"stateNode\">\n"
	if ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
		xml << "                <trait id=\"dateTrait.#{analysis_id}\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"date-backward\" units=\"year\" value=\"\n"
		branch.each do |b|
			if b.fulfilled_criteria.include?(2)
				xml << "                    #{b.id.sub("b","f")} = #{b.first_occurrence_age},\n"
			end
		end
		xml.chomp!(",\n")
		xml << "\">\n"
		xml << "                    <taxa id=\"taxonSet\" spec=\"beast.evolution.alignment.TaxonSet\" alignment=\"@#{analysis_id}\"/>\n"
		xml << "                </trait>\n"
	elsif scheme_id == "FBD_Range3"
		xml << "                <trait id=\"dateTrait.#{analysis_id}\" spec=\"beast.evolution.tree.TraitSet\" traitname=\"date-backward\" units=\"year\" value=\"\n"
		branch.each do |b|
			if b.fulfilled_criteria.include?(1)
				b.fossil.size.times do |x|
					xml << "                    #{b.id.sub("b","f")}_#{x} = #{b.fossil[x].age},\n"
				end
			end
		end
		xml.chomp!(",\n")
		xml << "\">\n"
		xml << "                    <taxa id=\"taxonSet\" spec=\"beast.evolution.alignment.TaxonSet\" alignment=\"@#{analysis_id}\"/>\n"
		xml << "                </trait>\n"		
	end
	xml << "                <taxonset id=\"TaxonSet.#{analysis_id}\" spec=\"TaxonSet\">\n"
	xml << "                    <data idref=\"#{analysis_id}\" name=\"alignment\"/>\n"
	xml << "                </taxonset>\n"
	xml << "            </tree>\n"
	if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
		xml << "            <parameter id=\"birthRate.t:tree\" name=\"stateNode\" lower=\"0.0\" upper=\"100.0\" value=\"1.0\"/>\n"
		xml << "            <parameter id=\"relativeDeathRate.t:tree\" name=\"stateNode\" lower=\"0.0\" upper=\"1.0\" value=\"0.5\"/>\n"
	end
	xml << "        </state>\n"
	xml << "\n"
	xml << "        <distribution id=\"posterior\" spec=\"util.CompoundDistribution\">\n"
	xml << "            <distribution id=\"prior\" spec=\"util.CompoundDistribution\">\n"
	xml << "\n"
	if scheme_id == "FBD" # sampling proportion and diversification rate fixed. # This uses the reparameterization of Gavryushkina et al. 2014, equation 8.
		xml << "                <distribution id=\"FBD.t:tree1\" spec=\"beast.evolution.speciation.SABirthDeathModel\" conditionOnRhoSampling=\"true\" conditionOnRoot=\"true\" tree=\"@Tree.t:tree1\">\n"
		xml << "                    <parameter id=\"turnoverFBD.t:tree1\" lower=\"0.0\" name=\"turnover\" upper=\"1.0\">#{true_turnover}</parameter>\n"
		xml << "                    <parameter id=\"samplingProportionFBD.t:tree1\" lower=\"0.0\" name=\"samplingProportion\" upper=\"1.0\">#{true_sampling_proportion}</parameter>\n"
		xml << "                    <parameter id=\"diversificationRateFBD.t:tree1\" lower=\"0.0\" name=\"diversificationRate\" upper=\"1.0\">#{true_net_diversification_rate}</parameter>\n"
		xml << "                    <parameter id=\"rFBD.t:tree1\" lower=\"0.0\" name=\"removalProbability\" upper=\"1.0\">0.0</parameter>\n"
		xml << "                    <parameter id=\"rhoFBD.t:tree1\" estimate=\"false\" lower=\"0.0\" name=\"rho\" upper=\"1.0\">#{n_sampled_species/n_extant_species.to_f}</parameter>\n"
		xml << "                </distribution>\n"
		xml << "\n"
	elsif scheme_id == "FBX" # sampling proportion estimated and diversification rate fixed.
		xml << "                <distribution id=\"FBD.t:tree1\" spec=\"beast.evolution.speciation.SABirthDeathModel\" conditionOnRhoSampling=\"true\" conditionOnRoot=\"true\" samplingProportion=\"@samplingProportionFBD.t:tree1\" tree=\"@Tree.t:tree1\">\n"
		xml << "                    <parameter id=\"turnoverFBD.t:tree1\" lower=\"0.0\" name=\"turnover\" upper=\"1.0\">#{true_turnover}</parameter>\n"
		xml << "                    <parameter id=\"diversificationRateFBD.t:tree1\" lower=\"0.0\" name=\"diversificationRate\" upper=\"1.0\">#{true_net_diversification_rate}</parameter>\n"
		xml << "                    <parameter id=\"rFBD.t:tree1\" lower=\"0.0\" name=\"removalProbability\" upper=\"1.0\">0.0</parameter>\n"
		xml << "                    <parameter id=\"rhoFBD.t:tree1\" estimate=\"false\" lower=\"0.0\" name=\"rho\" upper=\"1.0\">#{n_sampled_species/n_extant_species.to_f}</parameter>\n"
		xml << "                </distribution>\n"
		xml << "\n"
	elsif scheme_id == "FBY" # sampling proportion fixed and diversification rate estimated.
		xml << "                <distribution id=\"FBD.t:tree1\" spec=\"beast.evolution.speciation.SABirthDeathModel\" conditionOnRhoSampling=\"true\" conditionOnRoot=\"true\" diversificationRate=\"@diversificationRateFBD.t:tree1\" tree=\"@Tree.t:tree1\">\n"
		xml << "                    <parameter id=\"turnoverFBD.t:tree1\" lower=\"0.0\" name=\"turnover\" upper=\"1.0\">#{true_turnover}</parameter>\n"
		xml << "                    <parameter id=\"samplingProportionFBD.t:tree1\" lower=\"0.0\" name=\"samplingProportion\" upper=\"1.0\">#{true_sampling_proportion}</parameter>\n"
		xml << "                    <parameter id=\"rFBD.t:tree1\" lower=\"0.0\" name=\"removalProbability\" upper=\"1.0\">0.0</parameter>\n"
		xml << "                    <parameter id=\"rhoFBD.t:tree1\" estimate=\"false\" lower=\"0.0\" name=\"rho\" upper=\"1.0\">#{n_sampled_species/n_extant_species.to_f}</parameter>\n"
		xml << "                </distribution>\n"
		xml << "\n"
	elsif scheme_id == "FBZ" # sampling proportion estimated and diversification rate estimated from.
		xml << "                <distribution id=\"FBD.t:tree1\" spec=\"beast.evolution.speciation.SABirthDeathModel\" conditionOnRhoSampling=\"true\" conditionOnRoot=\"true\" samplingProportion=\"@samplingProportionFBD.t:tree1\" diversificationRate=\"@diversificationRateFBD.t:tree1\" tree=\"@Tree.t:tree1\">\n"
		xml << "                    <parameter id=\"turnoverFBD.t:tree1\" lower=\"0.0\" name=\"turnover\" upper=\"1.0\">#{true_turnover}</parameter>\n"
		xml << "                    <parameter id=\"rFBD.t:tree1\" lower=\"0.0\" name=\"removalProbability\" upper=\"1.0\">0.0</parameter>\n"
		xml << "                    <parameter id=\"rhoFBD.t:tree1\" estimate=\"false\" lower=\"0.0\" name=\"rho\" upper=\"1.0\">#{n_sampled_species/n_extant_species.to_f}</parameter>\n"
		xml << "                </distribution>\n"
		xml << "\n"
	elsif scheme_id == "FBD_Range1" # sampling proportion fixed, diversification rate and turnover estimated from ranges.
		xml << "                <distribution id=\"FBD.t:tree1\" spec=\"beast.evolution.speciation.SABirthDeathModel\" conditionOnRhoSampling=\"true\" conditionOnRoot=\"true\" diversificationRate=\"@diversificationRateFBD.t:tree1\" turnover=\"@turnoverFBD.t:tree1\" tree=\"@Tree.t:tree1\">\n"
		xml << "                    <parameter id=\"samplingProportionFBD.t:tree1\" lower=\"0.0\" name=\"samplingProportion\" upper=\"1.0\">#{true_sampling_proportion}</parameter>\n"
		xml << "                    <parameter id=\"rFBD.t:tree1\" lower=\"0.0\" name=\"removalProbability\" upper=\"1.0\">0.0</parameter>\n"
		xml << "                    <parameter id=\"rhoFBD.t:tree1\" estimate=\"false\" lower=\"0.0\" name=\"rho\" upper=\"1.0\">#{n_sampled_species/n_extant_species.to_f}</parameter>\n"
		xml << "                </distribution>\n"
		xml << "\n"
	elsif scheme_id == "FBD_Range2" or scheme_id == "FBD_Range3" # sampling proportion, diversification rate, and turnover estimated from ranges.
		xml << "                <distribution id=\"FBD.t:tree1\" spec=\"beast.evolution.speciation.SABirthDeathModel\" conditionOnRhoSampling=\"true\" conditionOnRoot=\"true\" diversificationRate=\"@diversificationRateFBD.t:tree1\" turnover=\"@turnoverFBD.t:tree1\" samplingProportion=\"@samplingProportionFBD.t:tree1\" tree=\"@Tree.t:tree1\">\n"
		xml << "                    <parameter id=\"rFBD.t:tree1\" lower=\"0.0\" name=\"removalProbability\" upper=\"1.0\">0.0</parameter>\n"
		xml << "                    <parameter id=\"rhoFBD.t:tree1\" estimate=\"false\" lower=\"0.0\" name=\"rho\" upper=\"1.0\">#{n_sampled_species/n_extant_species.to_f}</parameter>\n"
		xml << "                </distribution>\n"
		xml << "\n"
	end
	xml << "                <!--                                  Stem clades with fossil constraints                                     -->\n"
	xml << "\n"

	number_of_constrained_branches = 0
	unless scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
		# Prepare an array of all extant species that will be needed multiple times later (only with FBD schemes)
		extant_species = []
		branch.each do |b|
			if b.extant
				if b.species_id.include?("/")
					extant_species << b.species_id.split("/")[-1]
				else
					extant_species << b.species_id
				end
			end
		end
	end
	branch.each do |b|
		unless scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
			# We need to know which extant species descend from this branch, which fossil species are within the clade defined by this branch, and which extant species are not part of this clade.
			extant_species_of_this_branch = []
			fossil_species_of_this_branch = [] # This actually means of this branch or any descendant branch.
			extant_species_not_of_this_branch = []
			if ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
				unless b.first_occurrence_age == nil
					fossil_species_of_this_branch << b.id.sub("b","f") if b.fulfilled_criteria.include?(2)
				end
			elsif scheme_id == "FBD_Range3"
				if b.fulfilled_criteria.include?(1)
					b.fossil.size.times do |x|
						fossil_species_of_this_branch << "#{b.id.sub("b","f")}_#{x}"
					end
				end
			else
				puts "ERROR: Unknown scheme #{scheme_id}!"
				exit 1
			end

			if b.extant
				if b.species_id.include?("/")
					extant_species_of_this_branch << b.species_id.split("/")[-1]
				else
					extant_species_of_this_branch << b.species_id
				end
			end

			if ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
				unless b.extant
					b.progeny_id.each do |p|
						branch.each do |bb|
							if bb.id == p
								if bb.extant
									if bb.species_id.include?("/")
										extant_species_of_this_branch << bb.species_id.split("/")[-1]
									else
										extant_species_of_this_branch << bb.species_id
									end
								end
								unless bb.first_occurrence_age == nil
									fossil_species_of_this_branch << bb.id.sub("b","f") if bb.fulfilled_criteria.include?(2)
								end
							end
						end
					end
				end
			elsif scheme_id == "FBD_Range3"
				unless b.extant
					b.progeny_id.each do |p|
						branch.each do |bb|
							if bb.id == p
								if bb.extant
									if bb.species_id.include?("/")
										extant_species_of_this_branch << bb.species_id.split("/")[-1]
									else
										extant_species_of_this_branch << bb.species_id
									end
								end
								unless bb.first_occurrence_age == nil
									bb.fossil.size.times do |x|
										fossil_species_of_this_branch << "#{bb.id.sub("b","f")}_#{x}"
									end
								end
							end
						end
					end
				end

			else
				puts "ERROR: Unknown scheme #{scheme_id}!"
				exit 1
			end

			extant_species.each do |s|
				unless extant_species_of_this_branch.include?(s)
					extant_species_not_of_this_branch << s
				end
			end
			if extant_species_of_this_branch.size + fossil_species_of_this_branch.size > 1
				xml << "                <distribution id=\"#{b.id}_SA.prior\" stronglyMonophyletic=\"false\" spec=\"beast.evolution.tree.CladeConstraint\" tree=\"@Tree.t:tree1\">\n"
				xml << "                    <taxonsetIn id=\"#{b.id}_SA\" spec=\"TaxonSet\">\n"
				if ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
					extant_species_of_this_branch.each {|s| xml << "                        <taxon idref=\"" + "#{s}\"".ljust(8) + "/>\n"}
					fossil_species_of_this_branch.each {|f| xml << "                        <taxon idref=\"" + "#{f}\"".ljust(8) + "/>\n"}
				elsif scheme_id == "FBD_Range3"
					extant_species_of_this_branch.each {|s| xml << "                        <taxon idref=\"" + "#{s}\"".ljust(12) + "/>\n"}
					fossil_species_of_this_branch.each {|f| xml << "                        <taxon idref=\"" + "#{f}\"".ljust(12) + "/>\n"}
				end
				xml << "                    </taxonsetIn>\n"
				xml << "                    <taxonsetOut spec=\"TaxonSet\">\n"
				if ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
					extant_species_not_of_this_branch.each {|s| xml << "                        <taxon idref=\"" + "#{s}\"".ljust(8) + "/>\n"}
				elsif scheme_id == "FBD_Range3"
					extant_species_not_of_this_branch.each {|s| xml << "                        <taxon idref=\"" + "#{s}\"".ljust(12) + "/>\n"}
				end
				xml << "                    </taxonsetOut>\n"
				xml << "                </distribution>\n"
			end
		end
		constrain_this_branch = false
		constrain_this_branch = true if b.fulfilled_criteria.include?(11) and ["CladeAge", "CladeAge_Range"].include?(scheme_id)
		constrain_this_branch = true if b.fulfilled_criteria.include?(2) and ["FBD", "FBX", "FBY", "FBZ", "FBD_Range1", "FBD_Range2"].include?(scheme_id)
		constrain_this_branch = true if b.fulfilled_criteria.include?(1) and ["FBD_Range3"].include?(scheme_id)
		if constrain_this_branch
			number_of_constrained_branches += 1
			if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
				xml << "                <distribution id=\"#{b.id}_stem.prior\" monophyletic=\"true\" spec=\"beast.math.distributions.FossilPrior\" tree=\"@Tree.t:tree1\">\n"
				xml << "                    <taxonset id=\"#{b.id}_stem\" spec=\"TaxonSet\">\n"
				# Find out which species belong to this branch.
				extant_species_of_this_branch = []
				if b.extant
					if b.species_id.include?("/")
						extant_species_of_this_branch << b.species_id.split("/")[-1]
					else
						extant_species_of_this_branch << b.species_id
					end
				else
					b.progeny_id.each do |p|
						branch.each do |bb|
							if bb.id == p
								if bb.extant
									if bb.species_id.include?("/")
										extant_species_of_this_branch << bb.species_id.split("/")[-1]
									else
										extant_species_of_this_branch << bb.species_id
									end
								end
							end
						end
					end
				end
				if extant_species_of_this_branch.size == 0
					puts "ERROR: No extant species could be found of branch!"
					exit 1
				end

				# Add all species on this branch to a taxon set.
				extant_species_of_this_branch.each {|s| xml << "                        <taxon idref=\"" + "#{s}\"".ljust(8) + "/>\n"}

				# Add a CladeAge constraint for the taxon set.
				if b.first_occurrence_age == nil and b.first_occurrence_age_in_progeny == nil
					puts "ERROR: Both the first occurrence age and the first occurrence age in progeny are nil for branch #{b.id}!"
					exit 1
				elsif b.first_occurrence_age == nil
					constraint_age = b.first_occurrence_age_in_progeny
				elsif b.first_occurrence_age_in_progeny == nil
					constraint_age = b.first_occurrence_age
				else
					if b.first_occurrence_age > b.first_occurrence_age_in_progeny
						constraint_age = b.first_occurrence_age
					else
						constraint_age = b.first_occurrence_age_in_progeny
					end
				end
				xml << "                    </taxonset>\n"
                xml << "                    <fossilDistr\n"
                if scheme_id == "CladeAge"
	                xml << "                    	minOccuranceAge=\"#{constraint_age}\"\n"
					xml << "                    	maxOccuranceAge=\"#{constraint_age}\"\n"
	                xml << "                    	minSamplingRate=\"#{true_sampling_rate}\"\n"
	                xml << "                    	maxSamplingRate=\"#{true_sampling_rate}\"\n"
	                xml << "                    	minDivRate=\"#{true_net_diversification_rate}\"\n"
	                xml << "                    	maxDivRate=\"#{true_net_diversification_rate}\"\n"
	                xml << "                    	minTurnoverRate=\"#{true_turnover}\"\n"
	                xml << "                    	maxTurnoverRate=\"#{true_turnover}\"\n"
	                xml << "                    	minSamplingGap=\"0\"\n"
	                xml << "                    	maxSamplingGap=\"0\"\n"
	                xml << "                    	spec=\"beast.math.distributions.FossilCalibration\"/>\n"
				elsif scheme_id == "CladeAge_Range"
	                xml << "                    	minOccuranceAge=\"#{constraint_age}\"\n"
					xml << "                    	maxOccuranceAge=\"#{constraint_age}\"\n"
	                xml << "                    	minSamplingRate=\"#{0.5 * true_sampling_rate}\"\n"
	                xml << "                    	maxSamplingRate=\"#{1.5 * true_sampling_rate}\"\n"
	                xml << "                    	minDivRate=\"#{0.5 * true_net_diversification_rate}\"\n"
	                xml << "                    	maxDivRate=\"#{1.5 * true_net_diversification_rate}\"\n"
	                xml << "                    	minTurnoverRate=\"#{0.5 * true_turnover}\"\n"
	                xml << "                    	maxTurnoverRate=\"#{1.5 * true_turnover}\"\n"
	                xml << "                    	minSamplingGap=\"0\"\n"
	                xml << "                    	maxSamplingGap=\"0\"\n"
	                xml << "                    	spec=\"beast.math.distributions.FossilCalibration\"/>\n"
				end

                # Finalize this xml part.
                xml << "                </distribution>\n"
			end
			xml << "\n"
		end
	end

	xml << "\n"
	xml << "                <!--                                  Crown clades                                     -->\n"
	xml << "\n"
	crown_group_names = []
	branch.each do |b|
		unless b.extant
			crown_group_names << "#{b.id}_crown.prior"
			if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
				xml << "                <distribution id=\"#{b.id}_crown.prior\" monophyletic=\"true\" spec=\"beast.math.distributions.MRCAPrior\" tree=\"@Tree.t:tree1\" useOriginate=\"false\">\n"
			else
				xml << "                <distribution id=\"#{b.id}_crown.prior\" monophyletic=\"false\" spec=\"beast.math.distributions.MRCAPrior\" tree=\"@Tree.t:tree1\" useOriginate=\"false\">\n"
			end
			xml << "                    <taxonset id=\"#{b.id}_crown\" spec=\"TaxonSet\">\n"
			# Find out which species belong to this branch.
			extant_species_of_this_branch = []
			b.progeny_id.each do |p|
				branch.each do |bb|
					if bb.id == p
						if bb.extant
							if bb.species_id.include?("/")
								extant_species_of_this_branch << bb.species_id.split("/")[-1]
							else
								extant_species_of_this_branch << bb.species_id
							end
						end
					end
				end
			end
			if extant_species_of_this_branch.size == 0
				raise "No extant species could be found of branch #{b.id}!"
			end
			extant_species_of_this_branch.each {|s| xml << "                        <taxon idref=\"" + "#{s}\"".ljust(8) + "/>\n"}
			xml << "                    </taxonset>\n"
			xml << "                </distribution>\n"
			xml << "\n"
		end
	end
	if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
		xml << "                <distribution spec=\"beast.evolution.speciation.BirthDeathGernhard08Model\" id=\"BirthDeath.t:tree\" tree=\"@Tree.t:tree1\"\n"
		xml << "                    birthDiffRate=\"@birthRate.t:tree\"\n"
		xml << "                    relativeDeathRate=\"@relativeDeathRate.t:tree\"\n"
		xml << "                />\n"
		xml << "\n"
	end
	xml << "                <plate var=\'n\' range=\'&splitpartitions;\'>\n"
	xml << "                    <distribution count=\"@RBcount.s:$(n)\" distr=\"@Gamma.0\" id=\"RBprior.s:$(n)\" spec=\"beast.math.distributions.RBPrior\" x=\"@RBrates.s:$(n)\"/>\n"
	xml << "                    <prior id=\"ClockPrior.c:$(n)\" name=\"distribution\" x=\"@clockRate.c:$(n)\" distr=\'@OneOnX.0\'/>\n"
	xml << "                    <prior distr=\"@Exponential.0\" id=\"GammaShapePrior.s:$(n)\" name=\"distribution\" x=\"@gammaShape.s:$(n)\"/>\n"
	xml << "                </plate>\n"
	xml << "\n"
	if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
		xml << "                <prior id=\"BirthRatePrior.t:tree\" name=\"distribution\" x=\"@birthRate.t:tree\">\n"
		xml << "                    <Uniform id=\"Uniform.0\" name=\"distr\" upper=\"1000.0\"/>\n"
		xml << "                </prior>\n"
		xml << "                <prior id=\"DeathRatePrior.t:tree\" name=\"distribution\" x=\"@relativeDeathRate.t:tree\">\n"
		xml << "                    <Uniform id=\"Uniform.01\" name=\"distr\"/>\n"
		xml << "                </prior>\n"
		xml << "\n"
	end
	xml << "                <prior id=\"ucldStdevPrior.c:clock\" name=\"distribution\" x=\"@ucldStdev.c:clock\">\n"
	xml << "                    <Exponential id=\"Exponential.01\" name=\"distr\">\n"
	xml << "                        <parameter estimate=\"false\" id=\"RealParameter.05\" name=\"mean\" value=\"0.3333\"/>\n"
	xml << "                    </Exponential>\n"
	xml << "                </prior>\n"
	xml << "\n"
	if scheme_id == "FBX" or scheme_id == "FBZ"
		xml << "                <prior id=\"samplingProportionFBDPrior.t:tree1\" name=\"distribution\" x=\"@samplingProportionFBD.t:tree1\">\n"
		xml << "                    <Uniform id=\"Uniform.02\" name=\"distr\"/>\n"
		xml << "                </prior>\n"
		xml << "\n"
	end
	if scheme_id == "FBY" or scheme_id == "FBZ"
		xml << "                <prior id=\"diversificationRateFBDPrior.t:tree1\" name=\"distribution\" x=\"@diversificationRateFBD.t:tree1\">\n"
		xml << "                    <Uniform id=\"Uniform.03\" name=\"distr\"/>\n"
		xml << "                </prior>\n"
		xml << "\n"
	end
	if ["FBD_Range1", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "                <prior id=\"diversificationRateFBDPrior.t:tree1\" name=\"distribution\" x=\"@diversificationRateFBD.t:tree1\">\n"
		xml << "                    <Uniform id=\"Uniform.04\" name=\"distr\" lower=\"#{0.5 * true_net_diversification_rate}\" upper=\"#{1.5 * true_net_diversification_rate}\"/>\n"
		xml << "                </prior>\n"
		xml << "\n"
	end
	if ["FBD_Range1", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "                <prior id=\"turnoverFBDPrior.t:tree1\" name=\"distribution\" x=\"@turnoverFBD.t:tree1\">\n"
		xml << "                    <Uniform id=\"Uniform.05\" name=\"distr\" lower=\"#{0.5 * true_turnover}\" upper=\"#{1.5 * true_turnover}\"/>\n"
		xml << "                </prior>\n"
		xml << "\n"
	end
	if ["FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "                <prior id=\"samplingProportionFBDPrior.t:tree1\" name=\"distribution\" x=\"@samplingProportionFBD.t:tree1\">\n"
		xml << "                    <Uniform id=\"Uniform.06\" name=\"distr\" lower=\"#{0.5 * true_sampling_proportion}\" upper=\"#{1.5 * true_sampling_proportion}\"/>\n"
		xml << "                </prior>\n"
		xml << "\n"
	end
	xml << "            </distribution>\n"
	xml << "\n"
	xml << "            <distribution id=\"likelihood\" spec=\"util.CompoundDistribution\" useThreads=\'true\'>\n"
	xml << "                <plate var=\'n\' range=\'&splitpartitions;\'>\n"
	xml << "                    <distribution data=\"@$(n)\" id=\"treeLikelihood.$(n)\" spec=\"TreeLikelihood\" tree=\"@Tree.t:tree1\">\n"
	xml << "                        <siteModel gammaCategoryCount=\"4\" id=\"SiteModel.s:$(n)\" shape=\"@gammaShape.s:$(n)\" spec=\"SiteModel\">\n"
	xml << "                            <parameter estimate=\"false\" id=\"mutationRate.s:$(n)\" name=\"mutationRate\" value=\"1.0\"/>\n"
	xml << "                            <parameter estimate=\"false\" id=\"proportionInvariant.s:$(n)\" lower=\"0.0\" name=\"proportionInvariant\" upper=\"1.0\" value=\"0.0\"/>\n"
	xml << "                            <substModel count=\"@RBcount.s:$(n)\" id=\"RB.s:$(n)\" rates=\"@RBrates.s:$(n)\" spec=\"RB\">\n"
	xml << "                                <frequencies data=\"@$(n)\" id=\"freqs.s:$(n)\" spec=\"Frequencies\"/>\n"
	xml << "                            </substModel>\n"
	xml << "                        </siteModel>\n"
	xml << "                        <branchRateModel clock.rate=\"@clockRate.c:$(n)\" id=\"RelaxedClock.c:$(n)\" rateCategories=\"@rateCategories.c:clock\" spec=\"beast.evolution.branchratemodel.UCRelaxedClockModel\" tree=\"@Tree.t:tree1\">\n"
	xml << "                            <distr idref=\'LogNormalDistributionModel.c:clock\'/>\n"
	xml << "                        </branchRateModel>\n"
	xml << "                    </distribution>\n"
	xml << "                </plate>\n"
	xml << "            </distribution>\n"
	xml << "\n"
	xml << "        </distribution>\n"
	xml << "\n"
	xml << "        <operator id=\"ucldStdevScaler.c:clock\" parameter=\"@ucldStdev.c:clock\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"10.0\"/>\n"
	xml << "        <operator id=\"CategoriesRandomWalk.c:clock\" parameter=\"@rateCategories.c:clock\" spec=\"IntRandomWalkOperator\" weight=\"30.0\" windowSize=\"1\"/>\n"
	xml << "        <operator id=\"CategoriesSwapOperator.c:clock\" intparameter=\"@rateCategories.c:clock\" spec=\"SwapOperator\" weight=\"30.0\"/>\n"
	xml << "        <operator id=\"CategoriesUniform.c:clock\" parameter=\"@rateCategories.c:clock\" spec=\"UniformOperator\" weight=\"30.0\"/>\n"
	xml << "\n"
	xml << "        <plate var=\'n\' range=\'&splitpartitions;\'>\n"
	xml << "            <operator id=\"StrictClockRateScaler.c:$(n)\" parameter=\"@clockRate.c:$(n)\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"3.0\"/>\n"
	xml << "            <operator count=\"@RBcount.s:$(n)\" id=\"RBOperator.s:$(n)\" rates=\"@RBrates.s:$(n)\" spec=\"RBOperator\" weight=\"1.0\"/>\n"
	xml << "            <operator count=\"@RBcount.s:$(n)\" id=\"RBratescaler.s:$(n)\" parameter=\"@RBrates.s:$(n)\" scaleFactor=\"0.5\" spec=\"RBScaleOperator\" weight=\"1.0\"/>\n"
	xml << "            <operator id=\"gammaShapeScaler.s:$(n)\" parameter=\"@gammaShape.s:$(n)\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" weight=\"0.1\"/>\n"
	xml << "        </plate>\n"
	xml << "\n"
	xml << "        <operator id=\"strictClockUpDownOperator.c:$(n)\" scaleFactor=\"0.75\" spec=\"UpDownOperator\" weight=\"100.0\">\n"
	xml << "            <plate var=\'n\' range=\'&splitpartitions;\'>\n"
	xml << "                <up idref=\"clockRate.c:$(n)\"/>\n"
	xml << "            </plate>\n"
	xml << "            <down idref=\"Tree.t:tree1\"/>\n"
	xml << "        </operator>\n"
	xml << "\n"
	if ["FBX", "FBZ", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "        <operator id=\"samplingPScalerFBD.t:tree1\" spec=\"ScaleOperator\" parameter=\"@samplingProportionFBD.t:tree1\" scaleFactor=\"0.75\" weight=\"10.0\"/>\n"
	end
	if ["FBY", "FBZ", "FBD_Range1", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "        <operator id=\"divRateScalerFBD.t:tree1\" spec=\"ScaleOperator\" parameter=\"@diversificationRateFBD.t:tree1\" scaleFactor=\"0.75\" weight=\"10.0\"/>\n"
	end
	if ["FBD_Range1", "FBD_Range2", "FBD_Range3"].include?(scheme_id)
		xml << "        <operator id=\"turnoverScalerFBD.t:tree1\" spec=\"ScaleOperator\" parameter=\"@turnoverFBD.t:tree1\" scaleFactor=\"0.75\" weight=\"10.0\"/>\n"
	end
	if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
		xml << "        <operator id=\"treeScaler.t:tree\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"@Tree.t:tree1\" weight=\"30.0\"/>\n"
		xml << "        <operator id=\"treeRootScaler.t:tree\" rootOnly=\"true\" scaleFactor=\"0.5\" spec=\"ScaleOperator\" tree=\"@Tree.t:tree1\" weight=\"30.0\"/>\n"
		xml << "        <operator id=\"UniformOperator.t:tree\" spec=\"Uniform\" tree=\"@Tree.t:tree1\" weight=\"100.0\"/>\n"
		xml << "        <operator id=\"BirthRateScaler.t:tree\" parameter=\"@birthRate.t:tree\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"10.0\"/>\n"
		xml << "        <operator id=\"DeathRateScaler.t:tree\" parameter=\"@relativeDeathRate.t:tree\" scaleFactor=\"0.75\" spec=\"ScaleOperator\" weight=\"10.0\"/>\n"
	else
		xml << "        <operator id=\"LeafToSAFBD.t:tree1\" spec=\"LeafToSampledAncestorJump\" tree=\"@Tree.t:tree1\" weight=\"10.0\"/>\n"
		xml << "        <operator id=\"SAWilsonBaldingFBD.t:tree1\" spec=\"SAWilsonBalding\" tree=\"@Tree.t:tree1\" weight=\"10.0\"/>\n"
		xml << "        <operator id=\"SAWideFBD.t:tree1\" spec=\"SAExchange\" isNarrow=\"false\" tree=\"@Tree.t:tree1\" weight=\"10.0\"/>\n"
		xml << "        <operator id=\"SANarrowFBD.t:tree1\" spec=\"SAExchange\" tree=\"@Tree.t:tree1\" weight=\"10.0\"/>\n"
		xml << "        <operator id=\"SAUniformOperatorFBD.t:tree1\" spec=\"SAUniform\" tree=\"@Tree.t:tree1\" weight=\"20.0\"/>\n"
		xml << "        <operator id=\"SATreeRootScalerFBD.t:tree1\" spec=\"SAScaleOperator\" rootOnly=\"true\" scaleFactor=\"0.95\" tree=\"@Tree.t:tree1\" weight=\"1.0\"/>\n"
		xml << "        <operator id=\"SATreeScalerFBD.t:tree1\" spec=\"SAScaleOperator\" scaleFactor=\"0.95\" tree=\"@Tree.t:tree1\" weight=\"3.0\"/>\n"
	end
	xml << "\n"
	xml << "        <logger fileName=\"#{full_scheme_id}.log\" id=\"tracelog\" logEvery=\"50000\" model=\"@posterior\" sanitiseHeaders=\"true\" sort=\"smart\">\n"
	xml << "            <log idref=\"posterior\"/>\n"
	xml << "            <log idref=\"likelihood\"/>\n"
	xml << "            <log idref=\"prior\"/>\n"
	xml << "            <plate var=\'n\' range=\'&splitpartitions;\'>\n"
	xml << "                <log idref=\"treeLikelihood.$(n)\"/>\n"
	xml << "                <log idref=\"clockRate.c:$(n)\"/>\n"
	xml << "                <parameter idref=\"RBrates.s:$(n)\" name=\"log\"/>\n"
	xml << "                <parameter idref=\"gammaShape.s:$(n)\" name=\"log\"/>\n"
	xml << "                <log idref=\"RBcount.s:$(n)\"/>\n"
	xml << "            </plate>\n"
	xml << "            <log id=\"TreeHeight.t:tree\" spec=\"beast.evolution.tree.TreeHeightLogger\" tree=\"@Tree.t:tree1\"/>\n"
	if scheme_id == "CladeAge" or scheme_id == "CladeAge_Range"
		xml << "            <log idref=\"BirthDeath.t:tree\"/>\n"
	else
		xml << "            <log idref=\"FBD.t:tree1\"/>\n"
		xml << "            <log idref=\"diversificationRateFBD.t:tree1\"/>\n"
		xml << "            <log idref=\"turnoverFBD.t:tree1\"/>\n"
		xml << "            <log idref=\"samplingProportionFBD.t:tree1\"/>\n"
		xml << "            <log id=\"SACountFBD.t:tree1\" spec=\"beast.evolution.tree.SampledAncestorLogger\" tree=\"@Tree.t:tree1\"/>\n"
	end
	xml << "            <log idref=\"ucldStdev.c:clock\"/>\n"
	xml << "            <log branchratemodel=\"@RelaxedClock.c:#{analysis_id}_1\" id=\"rate.c:clock\" spec=\"beast.evolution.branchratemodel.RateStatistic\" tree=\"@Tree.t:tree1\"/>\n"
	xml << "        </logger>\n"
	xml << "\n"
	xml << "        <logger fileName=\"#{full_scheme_id}.mrca\" logEvery=\"50000\">\n"
	xml << "            <log idref=\"TreeHeight.t:tree\"/>\n"
	crown_group_names.each {|cg| xml << "            <log idref=\"#{cg}\"/>\n"}
	xml << "        </logger>\n"
	xml << "\n"
	xml << "        <logger id=\"screenlog\" logEvery=\"5000\">\n"
	xml << "            <log idref=\"posterior\"/>\n"
	xml << "            <log idref=\"TreeHeight.t:tree\"/>\n"
	xml << "        </logger>\n"
	xml << "\n"
	# xml << "        <logger fileName=\"#{analysis_id}.trees\" id=\"treelog.t:tree\" logEvery=\"50000\" mode=\"tree\">\n"
	# xml << "            <log branchratemodel=\"@RelaxedClock.c:#{analysis_id}_1\" id=\"TreeWithMetaDataLogger.t:tree\" spec=\"beast.evolution.tree.TreeWithMetaDataLogger\" tree=\"@Tree.t:tree1\"/>\n"
	# xml << "        </logger>\n"
	# xml << "\n"

	xml << "        <init spec=\'beast.util.TreeParser\' id=\'tree2\' taxa=\'@#{analysis_id}\' initial=\"@Tree.t:tree1\" IsLabelledNewick=\"true\">\n"
	xml << "            <input name=\'newick\'>\n"
	xml << "                #{starting_tree_newick}\n"
	xml << "            </input>\n"
	xml << "        </init>\n"
	xml << "\n"
	xml << "    </run>\n"
	xml << "\n"
	xml << "</beast>\n"

	xml_file_name = "#{full_scheme_id}.xml"
	Dir.mkdir("#{analysis_dir_name}/#{replicate_id}") unless Dir.exists?("#{analysis_dir_name}/#{replicate_id}")
	xml_file = File.new("#{analysis_dir_name}/#{replicate_id}/#{xml_file_name}","w")
	xml_file.write(xml)
	xml_file.close

	# Feedback.
	puts "Wrote replicate #{replicate_id} with #{number_of_constrained_branches} fossil constraints."

	# Prepare the slurm scripts.
	slurm_string1 = ""
    slurm_string1 << "#!/bin/bash\n"
    slurm_string1 << "\n"
    slurm_string1 << "# Job name:\n"
    slurm_string1 << "#SBATCH --job-name=#{full_scheme_id[0..5]}\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Project:\n"
    slurm_string1 << "#SBATCH --account=nn9244k\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Wall clock limit:\n"
    slurm_string1 << "#\n"
    slurm_string1 << "#SBATCH --time=168:00:00\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Processor and memory usage:\n"
    slurm_string1 << "#SBATCH --ntasks-per-node=4\n"
    slurm_string1 << "#SBATCH --nodes=1\n"
    slurm_string1 << "#SBATCH --mem-per-cpu=5G\n"
    slurm_string1 << "#\n"
    slurm_string1 << "# Outfile:\n"
    slurm_string1 << "#SBATCH --output=#{full_scheme_id}_out.txt\n"
    slurm_string1 << "\n"
    slurm_string1 << "## Set up job environment:\n"
    slurm_string1 << "source /cluster/bin/jobsetup\n"
    slurm_string1 << "module load beagle/4.1\n"
    slurm_string1 << "module load beast2/2.4.2\n"
    slurm_string1 << "\n"
    slurm_string1 << "## Copy input files to the work directory:\n"
    slurm_string1 << "cp #{full_scheme_id}.* $SCRATCH\n"
    slurm_string1 << "\n"
    slurm_string1 << "## Run BEAST:\n"
    slurm_string1 << "cd $SCRATCH\n"
    slurm_string_start = "beast -threads 3 -seed #{rand(100000)} -beagle #{full_scheme_id}.xml\n"
    slurm_string_resume = "beast -threads 3 -seed #{rand(100000)} -beagle -resume #{full_scheme_id}.xml\n"
    slurm_string2 = ""
    slurm_string2 << "\n"
    slurm_string2 << "# Copy files back to the submission directory\n"
    slurm_string2 << "cp #{full_scheme_id}* $SUBMITDIR\n"
    slurm_string2 << "\n"
    slurm_start_file_name = "start.slurm"
    slurm_start_file = File.new("#{analysis_dir_name}/#{replicate_id}/#{slurm_start_file_name}","w")
    slurm_start_file.write(slurm_string1 + slurm_string_start + slurm_string2)
    slurm_start_file.close
    slurm_resume_file_name = "resume.slurm"
    slurm_resume_file = File.new("#{analysis_dir_name}/#{replicate_id}/#{slurm_resume_file_name}","w")
    slurm_resume_file.write(slurm_string1 + slurm_string_resume + slurm_string2)
    slurm_resume_file.close

end
