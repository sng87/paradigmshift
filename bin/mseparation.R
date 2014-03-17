#!/usr/bin/env Rscript
    args = commandArgs(TRUE)
    focus_gene = args[1]
    positive_file = args[2]
    negative_file = args[3]
    
    positive_data = read.table(positive_file)
    negative_data = read.table(negative_file)
    x_min = min(positive_data$V1, negative_data$V1)
    if (x_min > 0) {
        x_min = 0
    }
    x_max = max(positive_data$V1, negative_data$V1)
    if (x_max < 0) {
        x_max = 0
    }
    x_min = 1.4*x_min
    x_max = 1.4*x_max
    
    y_max = 1.1*max(density(positive_data$V1)$y, density(negative_data$V1)$y)
    
    pdf(paste(focus_gene, ".msepartion.pdf", sep = ""), height = 5, width = 5)
    plot(density(positive_data$V1), col = "red", xlim = c(x_min, x_max), ylim = c(0, y_max), main = focus_gene, xlab = "")
    par(new = T)
    plot(density(negative_data$V1), xlim = c(x_min, x_max), ylim = c(0, y_max), main = "", xlab = "", ylab = "")
    dev.off()
    
