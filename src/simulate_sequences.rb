# m_matschiner Wed Jun 20 12:28:37 CEST 2018

# Load the phylsim package.
$libPath = "./phylsim/"
require "./phylsim/main.rb"

# Get the command-line argument.
res_dir = ARGV[0]

# Get the names of the subdirectories in the result directory.
dir_entries = []
Dir.entries(res_dir).each{ |i| dir_entries << i unless i[0] == "." }
dir_entries.sort!

# For each of the subdirectories, read the reconstructed tree dump files and simulate sequences.
dir_entries.each do |i|

	["diversified","random"].each do |taxon_sampling_scheme|
		if taxon_sampling_scheme == "diversified"
			tree_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.dmp"
			phylip_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.phy"
			info_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.info"
			outgroup_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.outgroup.txt"
			groups_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.monophyleticGroups.txt"
			constraint_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.constraint.txt"
			report_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.branches.txt"
			dump_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.dmp"
		elsif taxon_sampling_scheme == "random"
			tree_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.dmp"
			phylip_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.phy"
			info_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.info"
			outgroup_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.outgroup.txt"
			groups_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.monophyleticGroups.txt"
			constraint_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.constraint.txt"
			report_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.branches.txt"
			dump_file_name = "#{res_dir}/#{i}/species.fossils.random_sampling.reconstructed.seqs.dmp"
		end

		# Read the tree file.
		tree = Tree.load(fileName = tree_file_name, verbose = false)

		# Simulate sequence evolution.
		tree.assignBranchRates(mean = 0.003, standardDeviation = 0.003, autoCorrelation = 100, nStepsPerTimeUnit = 10)
		tree.evolveSequences(3000, substitutionModel = "ecm", threads = 1, verbose = true)

		# Write output.
		tree.writeSequences(fileName = phylip_file_name, format = "phylip", false, false, overwrite = true, verbose = true)
		tree.info(fileName = info_file_name, overwrite = true, verbose = true)
		tree.outgroup(fileName = outgroup_file_name, overwrite = true)
		tree.monophyleticGroups(fileName = groups_file_name, overwrite = true)
		tree.constraint(fileName = constraint_file_name, overwrite = true)
		tree.branchReport(fileName = report_file_name, overwrite = true)
		tree.dump(fileName = dump_file_name, overwrite = true, verbose = true)
	end

end