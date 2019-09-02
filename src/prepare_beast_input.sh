# m_matschiner Wed Jun 20 23:32:28 CEST 2018

# Make the output directory if it doesn't exist yet.
mkdir -p ../res/beast

# Prepare beast input.
for analysis_scheme in CladeAge CladeAge_Range FBD FBD_Range1 FBD_Range2 FBD_Range3 CladeAgeRS CladeAge_RangeRS FBDRS FBD_Range1RS FBD_Range2RS
do
	mkdir -p ../res/beast/${analysis_scheme}/replicates
	ruby prepare_beast_input.rb -s ${analysis_scheme}
done
