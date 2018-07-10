# m_matschiner Wed Jun 20 23:32:28 CEST 2018

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/beast

# Prepare beast input.
for scheme in CladeAge CladeAge_Range FBD FBD_Range1 FBD_Range2 FBD_Range3 # CladeAge FBD FBX FBY FBZ
do
	mkdir -p ../res/beast/${scheme}/replicates
	ruby prepare_beast_input.rb -s ${scheme}
done
