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

# For each of the subdirectories, read the tree file, reconstruct the tree, reduce the fossil record, and report fossil ages.
dir_entries.each do |i|

	# Read the tree file.
	tree_file_name = "#{res_dir}/#{i}/species.fossils.reconstructed.dmp"
	tree = Tree.load(fileName = tree_file_name, verbose = false)

	# Simulate sequence evolution.
	tree.assignBranchRates(mean = 0.003, standardDeviation = 0.003, autoCorrelation = 100, nStepsPerTimeUnit = 10)
	tree.evolveSequences(3000, substitutionModel = "ecm", threads = 1, verbose = true)

	# Write output.
	tree.writeSequences(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.phy", format = "phylip", false, false, overwrite = true, verbose = true)
	tree.info(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.info", overwrite = true, verbose = true)
	tree.outgroup(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.outgroup.txt", overwrite = true)
	tree.monophyleticGroups(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.monophyleticGroups.txt", overwrite = true)
	tree.constraint(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.constraint.txt", overwrite = true)
	tree.branchReport(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.branches.txt", overwrite = true)
	tree.dump(fileName = "#{res_dir}/#{i}/species.fossils.reconstructed.seqs.dmp", overwrite = true, verbose = true)

end