load("subsampled.RM.out.RData")
head(subsampled.RM.out)

dat = subsampled.RM.out
rm(subsampled.RM.out)

head(dat)

#no its a lil window thingy
Window.Size = 150000

comm_by_chromosome = lapply(split(dat, dat$Chromosome), function(y) t(sapply(split(y, list(cut(y[,"End"], as.integer(max(y$End)/as.integer(Window.Size))))), function(x) table(x$Family))))

									     
# flatter and more readable... but much much slower solution.... also the df is sideways relative to brent's
chr_counts = function(df){
	fam_count_chr = data.frame(Family = factor(unique(df$Family)))
	fam_count_chr
	
	windows = seq(100000, max(df$End), 100000)
	
	for(seg in windows){
		fam_count_chr[,paste(seg)] = 0 

		subset_a = df[df$End < seg ,]

		subset = subset_a[df$End > (seg-100000), ]
		for(trans in fam_count_chr$Family){

			count = length(subset$Family[subset$Family == trans])
			fam_count_chr[fam_count_chr$Family == trans, paste(seg)]	= count

		}
	}

	return(fam_count_chr)
}
	

list_of_df = split(dat, dat$Chromosome)

chr_counts_df = lapply(list_of_df, chr_counts)

chr_counts_df

#########
#########
# other example count the number of TE from each family on each chromosome.
#########
#########

all_chr = unique(dat$Chromosome)
length(all_chr) #humans!

all_fams = factor(unique(dat$Family))

length(all_fams)
levels(all_fams)

comm_by_chromosome = data.frame(Family = all_fams)
comm_by_chromosome

for(chr in all_chr){
	comm_by_chromosome[,chr] = 0 

	subset = dat[dat$Chromosome == chr ,]

	for(trans in comm_by_chromosome$Family){
		count = length(subset$Family[subset$Family == trans])
		comm_by_chromosome[comm_by_chromosome$Family == trans, chr]	= count
	}

}

head(comm_by_chromosome) #dataframe with chr as columns and family counts in rows


#challenge - one line the above with a groupby
