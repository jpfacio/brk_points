library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

##########          ARGUMENT VALIDATION         ##########

# The script expects at least 7 arguments.
# If fewer are provided, execution is interrupted.

if (length(args) < 7) {
  stop("Usage: Rscript brk.R <Control_5> <KO_5> <Control_3> <KO_3> <Control_cov> <KO_cov> <Chromosome> [Output_prefix]")
}

##########          CUSTOM VARIABLES         ##########

# The first four arguments refer to the profile data (5' and 3'),
# which must be in .bedgraph format.
# Arguments 5 and 6 correspond to coverage data (also .bedgraph).
# Argument 7 is the selected chromosome.
# Argument 8 is optional and corresponds to the output prefix.

control_5_profile  <- read.delim(args[1], header = FALSE)
knockout_5_profile <- read.delim(args[2], header = FALSE)
control_3_profile  <- read.delim(args[3], header = FALSE)
knockout_3_profile <- read.delim(args[4], header = FALSE)
coverage_control   <- read.delim(args[5], header = FALSE)
coverage_knockout  <- read.delim(args[6], header = FALSE)

chrom <- args[7]
out_name <- ifelse(length(args) == 8, args[8], "output")

##########          CREATING DIRECTORIES          ##########

# Create 'Raw_Tables' and 'Bedgraph_Version' directories
# if they do not already exist.

if (!dir.exists("Raw_Tables")) dir.create("Raw_Tables")
if (!dir.exists("Bedgraph_Version")) dir.create("Bedgraph_Version")

##########          DATA PROCESSING          ##########

# For normalization, the sum of coverage reads is obtained.

coverage_control_sum  <- sum(coverage_control$V4)
coverage_knockout_sum <- sum(coverage_knockout$V4)

# Separate the specific chromosome from profile data.

control_5_profile  <- control_5_profile[control_5_profile[,1] == chrom, ]
knockout_5_profile <- knockout_5_profile[knockout_5_profile[,1] == chrom, ]
control_3_profile  <- control_3_profile[control_3_profile[,1] == chrom, ]
knockout_3_profile <- knockout_3_profile[knockout_3_profile[,1] == chrom, ]

# Normalize read counts (column 4 of bedgraph) by total coverage.

control_5_profile[,4]  <- control_5_profile[,4]  / coverage_control_sum
knockout_5_profile[,4] <- knockout_5_profile[,4] / coverage_knockout_sum
control_3_profile[,4]  <- control_3_profile[,4]  / coverage_control_sum
knockout_3_profile[,4] <- knockout_3_profile[,4] / coverage_knockout_sum

# Only positions with at least 10 raw reads will be used.
# After normalization, cutoff is scaled accordingly.

cutoff_counts <- 10 / min(coverage_control_sum, coverage_knockout_sum)

control_5_profile  <- control_5_profile[control_5_profile[,4] >= cutoff_counts, ]
knockout_5_profile <- knockout_5_profile[knockout_5_profile[,4] >= cutoff_counts, ]
control_3_profile  <- control_3_profile[control_3_profile[,4] >= cutoff_counts, ]
knockout_3_profile <- knockout_3_profile[knockout_3_profile[,4] >= cutoff_counts, ]

##########          MAIN LOOP          ##########

# w (window) is the maximum tolerated interval to find break points.
# out will store the final result.
# positions corresponds to the number of 5' control positions.

w <- 3
out <- list()
positions <- nrow(control_5_profile)

# Main loop:
# Iterate over each 5' control position and search for
# a matching 3' signal within ±3 bp.

for (i in seq_len(positions)) {
  
  control_positions <- which(
    control_3_profile[,3] - w <= control_5_profile[i,3] &
    control_5_profile[i,3] <= control_3_profile[,3] + w
  )
  
  # If no 3' position is found in the window, skip iteration.
  if (length(control_positions) == 0) next
  
  # Select the 3' position with maximum read count.
  control_positions <- control_positions[which.max(control_3_profile[control_positions,4])]
  
  # Store control 5' and 3' information.
  pos5_c <- control_5_profile[i,3]
  val5_c <- control_5_profile[i,4]
  pos3_c <- control_3_profile[control_positions,3]
  val3_c <- control_3_profile[control_positions,4]
  delta  <- pos5_c - pos3_c
  
  # Find equivalent positions in knockout data.
  ko5_idx <- which(knockout_5_profile[,3] == pos5_c)
  ko3_idx <- which(knockout_3_profile[,3] == pos3_c)
  
  # If the same position does not exist in knockout,
  # use cutoff value for comparison.
  val5_k <- ifelse(length(ko5_idx) == 0, cutoff_counts, knockout_5_profile[ko5_idx,4])
  val3_k <- ifelse(length(ko3_idx) == 0, cutoff_counts, knockout_3_profile[ko3_idx,4])
  
  # Compute log2 fold changes.
  log2_5 <- log2(val5_c / val5_k)
  log2_3 <- log2(val3_c / val3_k)
  
  # Store break point information.
  out[[length(out) + 1]] <- c(
    pos5_c, val5_c,
    pos3_c, val3_c,
    delta,
    pos5_c, val5_k, log2_5,
    pos3_c, val3_k, log2_3
  )
}

# Stop execution if no break points were detected.
if (length(out) == 0) {
  stop("No break points detected.")
}

out <- do.call(rbind, out)

##########          FILTERING          ##########

# Apply log2 fold change filter.
# Delta must be between 1 and 5 bp.
# log2 fold change must be ≥ 1.

log_2_5_filtered <- (1 <= out[,5] & out[,5] <= 5) & (out[,8]  >= 1)
log_2_3_filtered <- (1 <= out[,5] & out[,5] <= 5) & (out[,11] >= 1)

out_log2_5 <- out[log_2_5_filtered, ]
out_log2_3 <- out[log_2_3_filtered, ]

##########          WRITING TABLES          ##########

# Define column names for clarity.

colnames(out) <- c(
  "Pos_5_C","N_5_C",
  "Pos_3_C","N_3_C",
  "Delta",
  "Pos_5_K","N_5_K","log2_5",
  "Pos_3_K","N_3_K","log2_3"
)

path_5 <- paste0("Raw_Tables/", out_name, "_5.csv")
path_3 <- paste0("Raw_Tables/", out_name, "_3.csv")

write.csv(out_log2_5, path_5, row.names = FALSE)
write.csv(out_log2_3, path_3, row.names = FALSE)

##########          BEDGRAPH PARSER          ##########

# Transform filtered CSV files into bedgraph format.

tables <- c(path_5, path_3)
beds   <- c(
  paste0("Bedgraph_Version/", out_name, "_5.bedgraph"),
  paste0("Bedgraph_Version/", out_name, "_3.bedgraph")
)

for (i in seq_along(tables)) {
  
  if (!file.exists(tables[i])) next
  
  df <- read.csv(tables[i])
  log_col <- ifelse(i == 1, "log2_5", "log2_3")
  
  df_bed <- df %>%
    mutate(
      chrom = chrom,
      start = Pos_5_C - 1,
      end   = Pos_5_C,
      score = .data[[log_col]]
    ) %>%
    select(chrom, start, end, score)
  
  write.table(
    df_bed,
    file = beds[i],
    sep = "\t",
    col.names = FALSE,
    row.names = FALSE,
    quote = FALSE
  )
}

########################################

print("Thank you for using Break Points. Don't forget to cite this work.")



