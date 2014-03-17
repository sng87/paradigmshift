#!/usr/bin/env Rscript
    args = commandArgs(TRUE)
    focus_gene = args[1]
    real_file = args[2]
    null_file = args[3]
    
    real_data = read.table(real_file)$V1
    null_data = read.table(null_file)$V1
    null_data = null_data[!is.nan(null_data)]
    min_val = min(real_data, null_data)
    if (min_val > 0) {
        min_val = 0
    }
    max_val = max(real_data, null_data)
    if (max_val < 0) {
        max_val = 0
    }
    min_val = 1.4*min_val
    max_val = 1.4*max_val

    pdf(paste(focus_gene, ".significance.pdf", sep = ""), height = 5, width = 5)
    plot(density(null_data), main = focus_gene, xlim = c(min_val, max_val))
    abline(v = real_data[1], col = "red")
    dev.off()
    
