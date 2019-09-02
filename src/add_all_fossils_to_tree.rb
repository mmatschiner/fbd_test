# m_matschiner Fri Jun 22 12:42:35 CEST 2018

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

	["diversified","random"].each do |taxon_sampling_scheme|
		if taxon_sampling_scheme == "diversified"
			tree_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.dmp"
			branch_report_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.branches.txt"
			fossil_tree_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.fbd.tre"
		elsif taxon_sampling_scheme == "random"
			tree_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.dmp"
			branch_report_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.branches.txt"
			fossil_tree_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.fbd_all.tre"
		end

		# Read the tree file.
		tree = Tree.load(fileName = tree_file_name, verbose = true)

		# Read the file with the branch report.
		branch_report_file = File.open(branch_report_file_name,"r")
		branch_report_lines = branch_report_file.readlines
		branch_report_branch_ids = []
		branch_report_lines.each do |l|
			if l[0..2] == "ID:"
				branch_report_branch_ids << l.split(":")[1].strip
			end
		end
		
		# Get array branch.
		branch = tree.branch

		# Check whether a branch fulfills criterion 1 (see prepare_BEAST_input.rb).
		branch_fulfills_criterion1 = []
		branch.size.times do |x|
			branch_fulfills_criterion1[x] = false
			branch_fulfills_criterion1[x] = true if branch[x].fossil.size > 0
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

		branch_count = 0
		new_node_ages = []
		branch.each do |b|
			if branch_fulfills_criterion1[branch_count]

				# Get first occurrence age.
				fossil_indices = []
				fossil_ages = []
				b.fossil.size.times do |x|
					fossil_indices << x
					fossil_ages << b.fossil[x].age
				end
				fossil_ages_sorted = false
				until fossil_ages_sorted
					fossil_ages_sorted = true
					0.upto(fossil_ages.size-2) do |x|
						(x+1).upto(fossil_ages.size-1) do |y|
							if fossil_ages[y] > fossil_ages[x]
								fossil_ages[y],fossil_ages[x] = fossil_ages[x],fossil_ages[y]
								fossil_indices[y],fossil_indices[x] = fossil_indices[x],fossil_indices[y]
								fossil_ages_sorted = false
								break
							end
						end
					end
				end

				# Make sure fossil ages are sorted.
				unless fossil_ages.sort.reverse == fossil_ages
					puts "Fossil ages don't seem to be sorted!"
					puts fossil_ages
					exit
				end

				new_branch_id_continuation = "b#{largest_branch_id_number+1+new_branch.size}"
				new_branch_id_fossil = "b#{largest_branch_id_number+2+new_branch.size}"

				# Determine the age of the new node.
				new_node_age = b.origin - 0.5 * (b.origin - [b.termination,fossil_ages[0]].max)
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

				end

				# Add the first fossil branch.
				fossil_branch = Branch.new(new_branch_id_fossil,new_node_age,fossil_ages[0],b.id,["none","none"],"extinction")
				fossil_branch.addSpeciesId("#{b.id.sub("b","f")}_#{fossil_indices[0]}")
				fossil_branch.updateRate(b.rate)
				fossil_branch.updateExtant(false)
				new_branch << fossil_branch

				# If multiple fossils exist, add all others.
				if fossil_ages.size > 1
					1.upto(fossil_ages.size-1) do |z|
						parent = new_branch.last
						parent.updateEndCause("speciation")
						original_parent_termination = parent.termination
						new_termination_age_for_parent = (parent.origin + parent.termination)/2.0
						parent.updateTermination(new_termination_age_for_parent)
						parent_continuation_branch_id = "b#{largest_branch_id_number+1+new_branch.size}"
						next_fossil_branch_id = "b#{largest_branch_id_number+2+new_branch.size}"
						parent.updateDaughterId([parent_continuation_branch_id,next_fossil_branch_id])
						parent_continuation_branch = Branch.new(parent_continuation_branch_id,new_termination_age_for_parent,original_parent_termination,parent.id,["none","none"],"extinction")
						parent_continuation_branch.addSpeciesId(parent.speciesId)
						parent_continuation_branch.updateRate(b.rate)
						parent_continuation_branch.updateExtant(false)
						new_branch << parent_continuation_branch
						if new_termination_age_for_parent < fossil_ages[z]
							puts "For branch #{next_fossil_branch_id}, the new termination age for parent #{new_termination_age_for_parent} branch is younger than the fossil age #{fossil_ages[z]}!"
							puts "Fossil ages are:"
							puts fossil_ages
							puts "Fossil indices are:"
							puts fossil_indices
							exit
						end
						next_fossil_branch = Branch.new(next_fossil_branch_id,new_termination_age_for_parent,fossil_ages[z],parent.id,["none","none"],"extinction")
						next_fossil_branch.addSpeciesId("#{b.id.sub("b","f")}_#{fossil_indices[z]}")
						next_fossil_branch.updateRate(b.rate)
						next_fossil_branch.updateExtant(false)
						new_branch << next_fossil_branch
					end
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

		# Make sure the branch has daughters if it has speciated.
		branch.each do |b|
			if b.endCause == "speciation" and b.daughterId == ["none","none"]
				puts "Branch has speciated but no daughters!"
				puts b
				exit
			end
		end

		# Save the tree with fossils as tips.
		tree.to_newick(fileName = fossil_tree_file_name, branchLengths = "duration", labels = true, plain = false, includeEmpty = true, overwrite = true, verbose = true)

	end

end
