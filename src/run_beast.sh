# m_matschiner Tue Jul 10 23:13:25 CEST 2018

# Run all beast analyses. This works only on a server with slurm, scripts need to be adapted to run it on other systems.
home=`pwd`
for analysis_dir in ../res/beast/{CladeAge,FBD}*/replicates/r???
do
	cd ${analysis_dir}
	sbatch start.slurm # Replace with resume.slurm if necessary.
	cd ${home}
done
