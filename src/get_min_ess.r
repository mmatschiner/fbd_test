# m_matschiner Wed Jul 4 22:44:00 CEST 2018

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

# Find the lowest effective sample size.
cat(min(effectiveSize(subset(mcmc, Sample > burnin, select = c(col_names_selected)))), "\n")