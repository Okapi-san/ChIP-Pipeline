#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args[1]
args[2]
b <- as.numeric(args[2])
x <- scan(args[1])
mat <- matrix(x, ncol = b, byrow = TRUE)
jpeg(file="filename.jpg")
image(mat)
dev.off()
q()
