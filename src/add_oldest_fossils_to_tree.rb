# m_matschiner Wed Jun 20 22:24:05 CEST 2018

# Load the phylsim package.
$libPath = "./phylsim/"
require "./phylsim/main.rb"

# Get the command-line argument.
res_dir = ARGV[0]

# Get the names of the subdirectories in the result directory.
dir_entries = []
Dir.entries(res_dir).each{ |i| dir_entries << i unless i[0] == "." }
dir_entries.sort!

# For each of the subdirectories, read the tree file, read the table with fossil information, and add the fossils to the tree as tips.
dir_entries.each do |i|

	# Read the tree file.
	tree_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.dmp"
	tree = Tree.load(fileName = tree_file_name, verbose = true)

	# Read the file with the branch report.
	branch_report_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.branches.txt"
	branch_report_file = File.open(branch_report_file_name,"r")
	branch_report_lines = branch_report_file.readlines
	branch_report_branch_ids = []
	branch_report_branch_first_occurrences = []
	branch_report_lines.each do |l|
		if l[0..2] == "ID:"
			branch_report_branch_ids << l.split(":")[1].strip
		elsif l[0..20] == "First occurrence age:"
			branch_report_branch_first_occurrences[branch_report_branch_ids.size-1] = l.split(":")[1].strip.to_f
		end
	end
		
	# Get array branch.
	branch = tree.branch

	# Add fossils to their branches.
	branch.each do |b|
		if branch_report_branch_ids.include?(b.id)
			unless branch_report_branch_first_occurrences[branch_report_branch_ids.index(b.id)] == nil
				fossils = [Fossil.new(b.id,branch_report_branch_first_occurrences[branch_report_branch_ids.index(b.id)])]
				b.addFossils(fossils)
			end
		else
			raise "ID of branch with fossil could not be found!"
		end
	end

	# Check whether a branch fulfills criterion 2 (see prepare_BEAST_input.rb).
	branch_fulfills_criterion2 = []
	branch_count = 0
	branch.each do |b|
		branch_fulfills_criterion2[branch_count] = false
		if b.fossil.size > 0
			# Get first occurrence in progeny.
			first_occurrence_age_in_progeny = nil
			branch.each do |bb|
				if b.progenyId.include?(bb.id) and bb.fossil.size > 0
					bb.fossil.each do |f|
						if first_occurrence_age_in_progeny == nil or f.age > first_occurrence_age_in_progeny
							first_occurrence_age_in_progeny = f.age
						end
					end
				end
			end
			# Get first occurrence on this branch.
			first_occurrence_age = 0
			b.fossil.each {|f| first_occurrence_age = f.age if first_occurrence_age < f.age}
			# This branch fulfills criterion 2 if its first occurrence age is greater than that of its progeny
			if first_occurrence_age_in_progeny == nil or first_occurrence_age > first_occurrence_age_in_progeny
				branch_fulfills_criterion2[branch_count] = true
			end
		end
		branch_count += 1
	end

	# Convert fossils into tips.
	new_branch = []

	# Get largest branch id.
	largest_branch_id_number = 0
	branch.each do |bb|
		if bb.id[1..-1].to_i > largest_branch_id_number
			largest_branch_id_number = bb.id[1..-1].to_i
		end
	end

	# Convert fossils to tips.
	branch_count = 0
	new_node_ages = []
	branch.each do |b|
		if branch_fulfills_criterion2[branch_count]

			# Get first occurrence age.
			first_occurrence_age = 0
			b.fossil.each {|f| first_occurrence_age = f.age if first_occurrence_age < f.age}

			new_branch_id_continuation = "b#{largest_branch_id_number+1+new_branch.size}"
			new_branch_id_fossil = "b#{largest_branch_id_number+2+new_branch.size}"

			# Determine the age of the new node.
			new_node_age = b.origin - 0.5 * (b.origin - [b.termination,first_occurrence_age].max)
			while new_node_ages.include?(new_node_age)
				new_node_age += 0.00000000001
			end
			new_node_ages << new_node_age

			# Memorize the termination age of the old branch.
			old_termination = b.termination

			if b.daughterId == ["none","none"]

				# Shorten the existing branch.
				b.updateTermination(new_node_age)
				b.updateEndCause("speciation")
				b.updateDaughterId([new_branch_id_continuation,new_branch_id_fossil])
				b.updateExtant(false)

				# Add an continuation branch.
				continuation_branch = Branch.new(new_branch_id_continuation,new_node_age,old_termination,b.id,["none","none"],"present")
				continuation_branch.addSpeciesId(b.speciesId)
				continuation_branch.updateRate(b.rate)
				continuation_branch.updateExtant(true)
				new_branch << continuation_branch

				# Add a fossil branch.
				fossil_branch = Branch.new(new_branch_id_fossil,new_node_age,first_occurrence_age,b.id,["none","none"],"extinction")
				fossil_branch.addSpeciesId(b.id.sub("b","f"))
				fossil_branch.updateRate(b.rate)
				fossil_branch.updateExtant(false)
				new_branch << fossil_branch
			
			else

				# Get daughters.
				unless b.daughterId == ["none","none"]
					daughter0 = nil
					daughter1 = nil
					branch.each do |d|
						if b.daughterId[0] == d.id
							daughter0 = d
						elsif b.daughterId[1] == d.id
							daughter1 = d
						end
					end
					raise "Daughter was not found!" if daughter0 == nil or daughter1 == nil
				end

				# Shorten the existing branch.
				b.updateTermination(new_node_age)
				b.updateEndCause("speciation")
				b.updateDaughterId([new_branch_id_continuation,new_branch_id_fossil])			

				# Add an continuation branch.
				continuation_branch = Branch.new(new_branch_id_continuation,new_node_age,old_termination,b.id,[daughter0.id,daughter1.id],"speciation")
				continuation_branch.addSpeciesId(b.speciesId)
				continuation_branch.updateRate(b.rate)
				continuation_branch.updateExtant(false)
				new_branch << continuation_branch

				# Update the parent Id of the two daughters.
				daughter0.updateParentId(new_branch_id_continuation)
				daughter1.updateParentId(new_branch_id_continuation)

				# Add a fossil branch.
				fossil_branch = Branch.new(new_branch_id_fossil,new_node_age,first_occurrence_age,b.id,["none","none"],"extinction")
				fossil_branch.addSpeciesId(b.id.sub("b","f"))
				fossil_branch.updateRate(b.rate)
				fossil_branch.updateExtant(false)
				new_branch << fossil_branch

			end

		end
		branch_count += 1
	end

	# Add the new branches to the old branch array.
	new_branch.each do |b|
		branch << b
	end

	# Make sure each branch id is unique.
	branch.each do |b|
		branch.each do |bb|
			unless b == bb
				if b.id == bb.id
					puts "Two branches have the same id!"
					puts "Branch 1:"
					puts b
					puts
					puts "Branch 2:"
					puts bb
					puts
					exit
				end
			end
		end
	end

	# Save the tree with fossils as tips.
	tree.to_newick(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.fbd.tre", branchLengths = "duration", labels = true, plain = false, includeEmpty = true, overwrite = true, verbose = true)

end
