#!/usr/bin/env Rscript
# type1error.R
# Kamil Slowikowski
# November 26, 2013

library(data.table)
library(reshape2)
library(ggplot2)
library(gap)

args <- commandArgs(TRUE)
base <- args[1]

# Read the null pvalues.
header <- c("name", "pvalue", "nulls_observed", "nulls_tested", "replicate")
null_pvalues <- read.delim(file.path(base, "null_pvalues.txt"),
                           header=F,
                           col.names=header)
null_pvalues <- data.table(null_pvalues)

# Find the number of p-values below various thresholds.
fun <- function(x, alpha) { sum(x < alpha) / length(x) }
df <- null_pvalues[, list(fun(pvalue, 0.5),
                        fun(pvalue, 0.1),
                        fun(pvalue, 0.05),
                        fun(pvalue, 0.01),
                        fun(pvalue, 0.005)), by=name]

# Prepare the data for plotting.
mf <- melt(df, id.vars=c("name"))
mf$variable <- as.character(mf$variable)
mf$variable[mf$variable == "V1"] <- "P < 0.5"
mf$variable[mf$variable == "V2"] <- "P < 0.1"
mf$variable[mf$variable == "V3"] <- "P < 0.05"
mf$variable[mf$variable == "V4"] <- "P < 0.01"
mf$variable[mf$variable == "V5"] <- "P < 0.005"
mf$variable <- factor(mf$variable, levels=unique(mf$variable))

N <- max(null_pvalues$replicate) + 1

ggplot(mf) +
    geom_hline(yintercept=c(0.5, 0.1, 0.05, 0.01, 0.005), alpha=0.5) +
    geom_point(aes(x=name, y=value, color=variable)) +
    scale_y_log10(breaks=c(0.5, 0.1, 0.05, 0.01, 0.005)) +
    scale_color_discrete(name="Threshold") +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    labs(title=paste(sep="",
                     "Proportion of p-values below each threshold for\n",
                     N, " random sets of SNPs"),
         y="Proportions",
         x="Conditions")

ggsave(file.path(base, "type1error.pdf"), width=6, height=5)
