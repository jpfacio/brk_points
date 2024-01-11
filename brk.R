args <- commandArgs(trailingOnly=T)

##########          CUSTOM VARIABLES         ##########

c3pfwd <- read.delim(args[1], header=F)
k3pfwd <- read.delim(args[2], header=F)
c5pfwd <- read.delim(args[3], header=F)
k5pfwd <- read.delim(args[4], header=F)
coverage_c <- read.delim(args[5], header=F)
coverage_k <- read.delim(args[6], header=F)
chrom <- args[7]
if (length(args) == 8) {out_name <- args[8]} else {out_name <- NULL}
#out_name <- NULL

##########          CREATING DIRECTORIES          ##########

if (system("test -d Raw_Tables", intern=T) == 1) {system("mkdir Raw_Tables")}

if (system("test -d Bedgraph_Version", intern=T) == 1) {system("mkdir Bedgraph_Version")}

##########          MAIN CODE          ##########

cN <- sum(coverage_c$V4)
kN <- sum(coverage_k$V4)

c5pfwd = c5pfwd[c5pfwd[,1] == chrom,]
k5pfwd = k5pfwd[k5pfwd[,1] == chrom,]
c3pfwd = c3pfwd[c3pfwd[,1] == chrom,]
k3pfwd = k3pfwd[k3pfwd[,1] == chrom,]

c5pfwd[,4] = signif(c5pfwd[,4]/cN, 2)
k5pfwd[,4] = signif(k5pfwd[,4]/kN, 2)
c3pfwd[,4] = signif(c3pfwd[,4]/cN, 2)
k3pfwd[,4] = signif(k3pfwd[,4]/kN, 2)

cutoff_counts = signif(10/min(cN,kN), 2)

c5pfwd = c5pfwd[c5pfwd[,4] >= cutoff_counts,]
k5pfwd = k5pfwd[k5pfwd[,4] >= cutoff_counts,]
c3pfwd = c3pfwd[c3pfwd[,4] >= cutoff_counts,]
k3pfwd = k3pfwd[k3pfwd[,4] >= cutoff_counts,]

w = 3

out = NULL
n = nrow(c5pfwd)

for (i in 1:n) {
  k = which(c3pfwd[, 3] - w <= c5pfwd[i, 3] & c5pfwd[i, 3] <= c3pfwd[, 3] + w)
  
  if (length(k) == 0) { 
    next 
  }
  
  k = k[which.max(c3pfwd[k, 4])]
  brk = c(c5pfwd[i, 3], c5pfwd[i, 4], c3pfwd[k, 3], c3pfwd[k, 4], c5pfwd[i, 3] - c3pfwd[k, 3])
  
  l = which(c5pfwd[i, 3] == k5pfwd[, 3])
  m = which(c3pfwd[k, 3] == k3pfwd[, 3])
  
  log5_cutoff = round(log2((c5pfwd[i, 4] / 1) / (cutoff_counts / 1)), 2)
  log5 = round(log2((c5pfwd[i, 4] / 1) / (k5pfwd[l, 4] / 1)), 2)
  
  log3_cutoff = round(log2((c3pfwd[k, 4] / 1) / (cutoff_counts / 1)), 2)
  log3 = round(log2((c3pfwd[k, 4] / 1) / (k3pfwd[m, 4] / 1)), 2)
  
  if (length(l) == 0) { 
    brk = c(brk, c5pfwd[i, 3], cutoff_counts, log5_cutoff)
  } else { 
    brk = c(brk, k5pfwd[l, 3], k5pfwd[l, 4], log5)
  }
  
  if (length(m) == 0) { 
    brk = c(brk, c3pfwd[k, 3], c3pfwd[k, 4], log3_cutoff)
  } else { 
    brk = c(brk, k3pfwd[m, 3], k3pfwd[m, 4], log3)
  }
  
  out = rbind(out, brk)
}


log2_5_f = (1 <= out[,5] & out[,5] <= 5) & (out[,8] >= log2(2))
log2_3_f = (1 <= out[,5] & out[,5] <= 5) & (out[,11] >= log2(2))

out_log2_5 = out[log2_5_f,]
out_log2_3 = out[log2_3_f,]

if (out_name == NULL) {
  path_5 <- "Raw_Tables/output_5.csv"
  path_3 <- "Raw_Tables/output_3.csv"
  
  path_5_bed <- "Bedgraph_Version/output_5.bedgraph"
  path_3_bed <- "Bedgraph_Version/output_3.bedgraph"
  
  write.table(out_log2_5, path_5, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_M',
                          'N_5_M', 'log2_5', 'Pos_3_M', 'N_3_M', 'log2_3'), row.names=F, quote=F)
  
  write.table(out_log2_3, path_3, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_M', 
                          'N_5_M', 'log2_5', 'Pos_3_M', 'N_3_M', 'log2_3'), row.names=F, quote=F)
} else {
  path_5 <- paste("Raw_Tables/", out_name, "_5.csv", sep="")
  path_3 <- paste("Raw_Tables/", out_name, "_3.csv", sep="")
  
  path_5_bed <- paste("Bedgraph_Version/", out_name, "_5.bedgraph", sep="")
  path_3_bed <- paste("Bedgraph_Version/", out_name, "_3.bedgraph", sep="")
  
  write.table(out_log2_5, path_5, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_M', 
                          'N_5_M', 'log2_5', 'Pos_3_M', 'N_3_M', 'log2_3'), row.names=F, quote=F)
  
  write.table(out_log2_3, path_3, sep=",", 
              col.names=c('Pos_5_C', 'N_5_C', 'Pos_3_C', 'N_3_C', 'Delta', 'Pos_5_M', 
                          'N_5_M', 'log2_5', 'Pos_3_M', 'N_3_M', 'log2_3'), row.names=F, quote=F)
}



##########          BEDGRAPH PARSER          ##########

library(dplyr)

tables <- c(path_5, path_3)

for (i in seq_along(tables)) {
  if (!file.exists(tables[i])) {
    cat("File not found:", tables[i], "\n")
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
