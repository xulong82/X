rm(list = ls())
setwd("~/Dropbox/edX")

# loess(): local polynomial regression fitting
# matrixStats::rowMads()
# genefilter::rowttests(), qvalue::qvalue()

# boxplot with formula: boxplot(x ~ idx) | boxplot(split(x, idx))

# SVD: U*D is distance between genes and D*V distance between samples
# Hereby, use svd$d * t(svd$v), instead of svd$v, to represent sample distance

# loess() do local fitting and smoothing
# predict() do prediction from results of various model fitting functions

# prop.table()
d = read.csv("assoctest.csv")
chisq.test(table(d))
fisher.test(table(d))

pop.var = var(bwt.nonsmoke)
vars = replicate(1000, var(sample(bwt.nonsmoke, 50)))
mean(vars > 1.5 * pop.var)
