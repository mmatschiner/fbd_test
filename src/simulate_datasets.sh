# m_matschiner Tue Jun 19 22:57:59 CEST 2018

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/datasets

# Simulate trees with the phylsim package.
ruby simulate_trees.rb ../res/datasets

# Add fossils to trees with the phylsim package.
ruby simulate_fossils.rb ../res/datasets

# Reconstruct the trees with the phylsim package.
ruby reconstruct_trees.rb ../res/datasets

# Simulate sequence evolution with the phylsim package.
ruby simulate_sequences.rb ../res/datasets

# Make a new tree string with fossils as tips.
ruby add_oldest_fossils_to_tree.rb ../res/datasets
ruby add_all_fossils_to_tree.rb ../res/datasets
