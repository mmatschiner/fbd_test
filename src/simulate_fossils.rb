# m_matschiner Tue Jun 19 23:19:03 CEST 2018

# Load the phylsim package.
$libPath = "./phylsim/"
require "./phylsim/main.rb"

# Get the command-line argument.
res_dir = ARGV[0]

# Get the names of the subdirectories in the result directory.
dir_entries = []
Dir.entries(res_dir).each{ |i| dir_entries << i unless i[0] == "." }
dir_entries.sort!

# For each of the subdirectories, read the tree file.
dir_entries.each do |i|
	tree_file_name = "#{res_dir}/#{i}/species.dmp"
	tree = Tree.load(fileName = tree_file_name, verbose = true)
	tree.addFossilRecord(0.01, samplingGap = 0, verbose = true)
	tree.dump(fileName = "#{res_dir}/#{i}/species.fossils.dmp", overwrite = true, verbose = false)
	tree.fossilReport(fileName = "#{res_dir}/#{i}/species.fossils.txt", overwrite = true, verbose = false)
end
