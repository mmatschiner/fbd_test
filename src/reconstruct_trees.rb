# m_matschiner Wed Jun 20 10:19:03 CEST 2018

# Load the phylsim package.
$libPath = "./phylsim/"
require "./phylsim/main.rb"

# Get the command-line argument.
res_dir = ARGV[0]

# Get the names of the subdirectories in the result directory.
dir_entries = []
Dir.entries(res_dir).each{ |i| dir_entries << i unless i[0] == "." }
dir_entries.sort!

# For each of the subdirectories, read the tree file, reconstruct the tree, reduce the fossil record, and report fossil ages.
dir_entries.each do |i|
	# Read the tree file.
	tree_file_name = "#{res_dir}/#{i}/species.fossils.dmp"
	["diversified","random"].each do |taxon_sampling_scheme|
		tree = Tree.load(fileName = tree_file_name, verbose = true)

		# Reconstruct the tree with the diversified and random sampling schemes.
		if taxon_sampling_scheme == "diversified"
			newick_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.tre"
			dump_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.dmp"
			fossil_report_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.txt"
		elsif taxon_sampling_scheme == "random"
			newick_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.tre"
			dump_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.dmp"
			fossil_report_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.txt"
		end
		tree.reconstruct(samplingScheme = taxon_sampling_scheme, number = 50, focusGroup = nil, verbose = true)
		tree.to_newick(fileName = newick_file_name, branchLengths = "duration", labels = true, plain = false, includeEmpty = true, overwrite = true, verbose = false)
		tree.dump(fileName = dump_file_name, overwrite = true, verbose = false)

		# Remove all fossils except the oldest per branch.
		tree.branch.each do |b|
			fossil = b.fossil
			if fossil.size > 0
				oldest_fossil = fossil[0]
				fossil.each do |f|
					oldest_fossil = f if f.age > oldest_fossil.age
				end
				b.updateFossil([oldest_fossil])
			end
		end

		# Prepare a table with the ages of the oldest fossils per branch.
		fossil_report_string = ""
		tree.branch.each do |b|
			fossil_report_string << "#{b.id}\t#{b.fossil[0].age}\n" if b.fossil.size > 0
		end

		# Write the report to a file.
		fossil_report_file = File.open(fossil_report_file_name, "w")
		fossil_report_file.write(fossil_report_string)
	end
end
