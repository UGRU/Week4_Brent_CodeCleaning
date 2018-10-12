load("subsampled.RM.out.RData")
head(subsampled.RM.out)

dat <- subsampled.RM.out
rm(subsampled.RM.out)

comm_by_chromosome <- lapply(split(dat, dat$Chromosome), function(y) t(sapply(split(y, list(cut(y[,"End"], as.integer(max(y$End)/as.integer(Window.Size))))), function(x) table(x$Family))))