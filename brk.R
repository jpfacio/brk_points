library(dplyr)

args <- commandArgs(trailingOnly=T)

##########          CUSTOM VARIABLES         ##########

# This session refers to the variables that will be passed as input by the user
# The first four arguments refers to the profile data, 5' and 3', which needs
# be in .bedgraph format. The args 5 and 6 expects the coverage data from the
# sequencing, which needs to be in .bedgraph format. The arg 7 and 8 are both
# strings, the first being the selected chromosome and the latter being the 
# prefix of the output (optional).

control_5_profile <- read.delim(args[1], header=F)
knockout_5_profile <- read.delim(args[2], header=F)
control_3_profile <- read.delim(args[3], header=F)
knockout_3_profile <- read.delim(args[4], header=F)
coverage_control <- read.delim(args[5], header=F)
coverage_knockout <- read.delim(args[6], header=F)
chrom <- args[7]
if (length(args) == 8) {out_name <- args[8]} else {out_name <- NULL}

##########          CREATING DIRECTORIES          ##########

# Create 'Raw_Tables' and 'Bedgraph_Version' directories in the working directory
# if they do not exist

if (!dir.exists("Raw_Tables")) {dir.create("Raw_Tables")}

if (!dir.exists("Bedgraph_Version")) {dir.create("Bedgraph_Version")}

##########          DATA PROCESSING          ##########

# For normalization, the sum of the coverage reads are obtained

coverage_control_sum <- sum(coverage_control$V4)
knockout_control_sum <- sum(coverage_knockout$V4)

# The specific chromosome is separated from profile data

control_5_profile = control_5_profile[control_5_profile[,1] == chrom,]
knockout_5_profile = knockout_5_profile[knockout_5_profile[,1] == chrom,]
control_3_profile = control_3_profile[control_3_profile[,1] == chrom,]
knockout_3_profile = knockout_3_profile[knockout_3_profile[,1] == chrom,]

# Based on the normalization variables (coverage_control_sum and knockout_control_sum),
# the fourth column of the .bedgraph (count of reads) are divided by the total sum 
# of coverage data rounded by two decimals

control_5_profile[,4] = signif(control_5_profile[,4]/coverage_control_sum, 2)
knockout_5_profile[,4] = signif(knockout_5_profile[,4]/knockout_control_sum, 2)
control_3_profile[,4] = signif(control_3_profile[,4]/coverage_control_sum, 2)
knockout_3_profile[,4] = signif(knockout_3_profile[,4]/knockout_control_sum, 2)

# Only positions with at least 10 reads will be used

cutoff_counts = signif(10/min(coverage_control_sum,knockout_control_sum), 2)

control_5_profile = control_5_profile[control_5_profile[,4] >= cutoff_counts,]
knockout_5_profile = knockout_5_profile[knockout_5_profile[,4] >= cutoff_counts,]
control_3_profile = control_3_profile[control_3_profile[,4] >= cutoff_counts,]
knockout_3_profile = knockout_3_profile[knockout_3_profile[,4] >= cutoff_counts,]

##########          MAIN LOOP          ##########

# This is the setup for the main loop, the w variable (window) is the maximum 
# interval tolerated to find the break point, out will store the final result
# of the loop to later processing, positions stores the number of rows in the 
# control 5' profile forward, which will be iterated for comparison 

w = 3
out = NULL
positions = nrow(control_5_profile)

# Main loop: Iterate over each end position (column 3 of the bedgraph) of control
# data. For each position, aggregate 3' profile values within a 3-base window 
# (positions-1 to positions+1). Results are stored in the list 'control_positions' 

for (i in 1:positions) {
  control_positions = which(control_3_profile[, 3] - w <= control_5_profile[i, 3] & 
                              control_5_profile[i, 3] <= control_3_profile[, 3] + w)
  
  # If the given position returns nothing, the loop continues
  
  if (length(control_positions) == 0) { 
    next 
  }
  
  # The control_positions variable will store 1 to 3 positions, this step takes the max value of
  # reads (column 4)
  
  control_positions = control_positions[which.max(control_3_profile[control_positions, 4])]
  
  # The brk variable consists of a vector containing the iterated position in the 
  # 5' profile and its associated value, the equivalent 3' position with its value
  # and the difference between those positions (delta).
  
  brk = c(control_5_profile[i, 3], control_5_profile[i, 4], control_3_profile[control_positions, 3], 
          control_3_profile[control_positions, 4], 
          control_5_profile[i, 3] - control_3_profile[control_positions, 3])
  
  # The knockout_positions_5 variable takes exactly the same position spotted in 
  # the control in the knockout data, from the 5' profile. The knockout_positions_3
  # variable do the same, but instead takes the position from the 3' profile (this 
  # will be used for further comparison)
  
  knockout_positions_5 = which(control_5_profile[i, 3] == knockout_5_profile[, 3])
  knockout_positions_3 = which(control_3_profile[control_positions, 3] == 
                                 knockout_3_profile[, 3])
  
  # These variables will store the ratio between the control and knockout, 
  # rounded by two decimals, pivoted in the 5' side. The log5_cutoff uses
  # the previously defined cutoff to compare, in cases where the control 
  # and knockout doesn't match in the identical position
  
  log5_cutoff = round(log2((control_5_profile[i, 4] / 1) / (cutoff_counts / 1)), 2)
  log5 = round(log2((control_5_profile[i, 4] / 1) / 
                      (knockout_5_profile[knockout_positions_5, 4] / 1)), 2)
  
  # Same logic but with the 3' side
  
  log3_cutoff = round(log2((control_3_profile[control_positions, 4] / 1) / (cutoff_counts / 1)), 2)
  log3 = round(log2((control_3_profile[control_positions, 4] / 1) / 
                      (knockout_3_profile[knockout_positions_3, 4] / 1)), 2)
  
  # If the same position isn't present in the knockout but it passed on the
  # script, use the cutoff as comparison
  
  if (length(knockout_positions_5) == 0) { 
    brk = c(brk, control_5_profile[i, 3], cutoff_counts, log5_cutoff)
  } else { 
    brk = c(brk, knockout_5_profile[knockout_positions_5, 3], 
            knockout_5_profile[knockout_positions_5, 4], log5)
  }
  
  # 3' side
  
  if (length(knockout_positions_3) == 0) { 
    brk = c(brk, control_3_profile[control_positions, 3], 
            control_3_profile[control_positions, 4], log3_cutoff)
  } else { 
    brk = c(brk, knockout_3_profile[knockout_positions_3, 3], 
            knockout_3_profile[knockout_positions_3, 4], log3)
  }
  
  # Defines the final vector, containing the break point positions and its 
  # log2 fold change
  
  out = rbind(out, brk)
}

##########          FILTERING          ##########

# Applying the log2 fold change filter, considering the delta value 
# obtained in the main vector. The number 5 stands for 5' side and
# the number 3 stands for 3' side.

log_2_5_filtered = (1 <= out[,5] & out[,5] <= 5) & (out[,8] >= log2(2))
log_2_3_filtered = (1 <= out[,5] & out[,5] <= 5) & (out[,11] >= log2(2))

# Saving a filtered vector for 5' and 3' side

out_log2_5 = out[log_2_5_filtered,]
out_log2_3 = out[log_2_3_filtered,]

##########          WRITING TABLES AND BEDGRAPH          ##########

# This section writes and saves the results in the created
# directories. The 5 and 3 in the variable names refers to
# the 5' and 3' side output, respectively.

if (is.null(out_name)) {
  path_5 <- "Raw_Tables/output_5.csv"
  path_3 <- "Raw_Tables/output_3.csv"
  
  path_5_bed <- "Bedgraph_Version/output_5.bedgraph"
  path_3_bed <- "Bedgraph_Version/output_3.bedgraph"
  
  # The columns names were shortened for cleaning purposes,
  # Pos means Positions, C means control, K for knockout
  # N are the number of reads in the position, 5 and 3
  # refers to 5' and 3' sides and log_2_5 and log_3_5
  # refers to the log2 fold change value in both sides
  
  write.table(out_log2_5, path_5, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_K',
                          'N_5_K', 'log2_5', 'Pos_3_K', 'N_3_K', 'log2_3'), 
              row.names=F, quote=F)
  
  write.table(out_log2_3, path_3, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_M', 
                          'N_5_K', 'log2_5', 'Pos_3_K', 'N_3_K', 'log2_3'),
              row.names=F, quote=F)
} else {
  
  # The paths are saved based on profile side, if the 
  # arg 8 is NULL, the default name is output, with
  # the suffix _5 and _3.
  
  path_5 <- paste("Raw_Tables/", out_name, "_5.csv", sep="")
  path_3 <- paste("Raw_Tables/", out_name, "_3.csv", sep="")
  
  path_5_bed <- paste("Bedgraph_Version/", out_name, "_5.bedgraph", sep="")
  path_3_bed <- paste("Bedgraph_Version/", out_name, "_3.bedgraph", sep="")
  
  write.table(out_log2_5, path_5, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_K', 
                          'N_5_K', 'log2_5', 'Pos_3_K', 'N_3_K', 'log2_3'), 
              row.names=F, quote=F)
  
  write.table(out_log2_3, path_3, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_M', 
                          'N_5_M', 'log2_5', 'Pos_3_M', 'N_3_M', 'log2_3'), 
              row.names=F, quote=F)
}



##########          BEDGRAPH PARSER          ##########

# This session transforms .csv data in .bedgraph

tables <- c(path_5, path_3)

for (i in seq_along(tables)) {
  if (!file.exists(tables[i])) {
    cat("File not found:", tables[i], "\positions")
    next
  }
  
  df <- read.csv(tables[i])
  log2_col <- ifelse(i == 1, "log2_5", "log2_3")
  
  df <- df %>%
    mutate(
      chrom = chrom,
      Initial_Pos = Pos_5_C - 1,
      Final_Pos = Pos_5_C - 1,
      log2 = as.numeric(.data[[log2_col]])
    ) %>%
    select(chrom, Initial_Pos, Final_Pos, log2)
  
  output_file <- ifelse(i == 1, path_5_bed, path_3_bed)
  
  write.table(df, file = output_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}

########################################

print("Thank you for using! Don't forget to cite this work.")
