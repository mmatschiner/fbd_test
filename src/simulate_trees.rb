# m_matschiner Tue Jun 19 22:59:57 CEST 2018

# Load the phylsim and the fileutils libraries.
$libPath = "./phylsim/"
require "./phylsim/main.rb"
require 'fileutils'

# Get the command-line argument.
res_dir = ARGV[0]

# Set the number of replicate species trees.
number_of_replicates = 20

number_of_replicates.times do |x|
	tree = Tree.generate(lambda = 0.12, mu = 0.06, treeOrigin = 100, present = 0, k = 0, rootSplit = true, np = [4000,5000], npEach = [1,'inf'], checkProbabilities = false, algorithm = "forward", verbose = true, threads = 1)
	dir_name = "#{res_dir}/r#{(x+1).to_s.rjust(3).gsub(" ","0")}"
	FileUtils.mkdir_p("#{dir_name}")
	tree.to_newick(fileName = "#{dir_name}/species.tre", branchLengths = "duration", labels = true, plain = false, includeEmpty = true, overwrite = true, verbose = false)
	tree.dump(fileName = "#{dir_name}/species.dmp", overwrite = false, verbose = false)
end