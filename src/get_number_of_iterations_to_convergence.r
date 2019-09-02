# m_matschiner Mon Oct 17 16:40:19 CEST 2016

# Suppress warnings.
options(warn=-1)

# Load the coda library.
library(coda)

# Read the log file specified as the first command line argument.
args <- commandArgs(trailingOnly = TRUE)
log_file_name <- args[1]
burnin <- as.integer(args[2])

# Read the log file as a table.
mcmc <- read.table(log_file_name,header=T)
frequency <- mcmc$Sample[2]

# Set constants.
# burnin_proportion <- 0.5
convergence_ess <- 200

# Get the column names.
col_names <- colnames(mcmc)

# Get the names of those columns that are actually estimated.
uniq_values_per_col <- apply(mcmc, 2, function(x) length(unique(x)))
col_names_not_estimated <- names(which(uniq_values_per_col == 1))
col_names_estimated <- setdiff(col_names, col_names_not_estimated)

# Exclude certain columns.
col_names_excluded <- col_names_estimated[grepl("RB", col_names_estimated)]
col_names_excluded <- append(col_names_excluded, "Sample")
col_names_excluded <- append(col_names_excluded, "rate.clock.variance")
col_names_selected = setdiff(col_names_estimated, col_names_excluded)

# Get the sample at which for the first time all sample sizes are greater than the convergence_ess value.
convergence_sample_size <- 0
for (last_sample in mcmc$Sample[1:length(mcmc$Sample)]){
	if (last_sample > burnin+frequency) {
		lowest_effective_sample_size <- min(effectiveSize(subset(mcmc, Sample > burnin & Sample <= last_sample, select = c(col_names_selected))))
		if (lowest_effective_sample_size > convergence_ess){
			convergence_sample_size <- last_sample
			break
		}
	}
}

# Report the result, or "NA" if no convergence could be found.
if (convergence_sample_size > 0){
	cat(convergence_sample_size, sep="\n")
} else {
	cat("NA", sep="\n")
}