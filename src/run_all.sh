# m_matschiner Tue Jul 10 23:11:54 CEST 2018

# Simulate datasets.
bash simulate_datasets.sh

# Prepare input for BEAST2.
bash prepare_beast_input.sh

# Run all BEAST2 analyses.
bash run_beast.sh

# Get run statistics.
bash get_run_statistics.sh

# Get age estimates.
bash extract_and_plot_age_estimates.sh

# Get parameter estimates.
bash extract_parameter_estimates.sh