library(data.table)
library(magrittr)

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
									     
# This lil pipe does the job a simpler
counts_chr = dat %>%
  group_by(Chromosome, Family) %>%
  count()

# version from Karl --------

library(tidyverse)

karl = dat %>% group_by(Chromosome) %>% 
  mutate(position = cut(End,as.integer(max(End)/as.integer(Window.Size)))) %>% 
  ungroup() %>% 
  group_by(Chromosome, position, Family) %>% 
  summarise(abundance = n()) %>% 
  spread(Family, abundance, fill = 0) %>% 
  ungroup() %>% 
  arrange(position) %>% 
  split(f = .$Chromosome)

# Jack's solution ---------

windows <- dat %>%
  mutate(window = round(End/100000)*100000) %>%
  group_by(Chromosome, window, Family) %>%
  summarise(obs = n()) %>%
  spread(Family, obs, fill = 0) %>% 
  ungroup() # I needed to add this line to work with the results afterwards

# compare the next line with the last line of the karl solution, to see how you can use base functions into a piping approach by using the "." notation as a placeholder for the end result of the previous pipe.
chromosomes <- split(windows, windows$Chromosome)

# comparing the 3 solutions ---------

# comparing column totals ----
comm_by_chromosome[[1]] %>% 
  colSums()

karl[[1]] %>% dplyr::select(-Chromosome, -position) %>% 
  colSums()

chromosomes[[1]] %>% dplyr::select(-Chromosome, -window) %>% 
  colSums()

# these are identical to Brent's original solution
# e.g. total SINE/MIR abundance is 64 in all solutions.
# except that certain columns without any elements are missing in karl and Jack's solution

# comparing row totals

comm_by_chromosome[[1]] %>% 
  rowSums()

karl[[1]] %>% dplyr::select(-Chromosome, -position) %>% 
  rowSums()

chromosomes[[1]] %>% dplyr::select(-Chromosome, -window) %>% 
  rowSums()

# here there are some discrepancies
# Brent's original solution maintains all windows, even if there are no elements found in them
# karl and Jack's solution are very similar, but due to small differences in how the intervals are computed, there some small differences in the row totals

# print Brent's solution w/o the empty rows
comm_by_chromosome[[1]] %>% 
  rowSums() %>% 
  unname() %>% # this line was necessary to help printing the sums without the window names
  .[. > 0] # another example of how the "." notation helps you with using base R functions

# looks very similar to the other solutions, so you will need to decide if the missing rows (and columns) are important or not
# in a community ecology analysis they would removed anyway, so in this case not.

# comparing speeds of the 3 solutions --------

library(microbenchmark)

# the code below will take some time. Repeat only if you want to replicate it yourself
# I just copy-pasted the above code inside the benchmark function


jack = dat %>%
  mutate(window = round(End/100000)*100000) %>%
  group_by(Chromosome, window, Family) %>%
  summarise(obs = n()) %>%
  spread(Family, obs, fill = 0) %>% 
  ungroup()

dat2 <- as.data.table(dat)
kennel <- dat2[, window := round(End, -5)] %>%
  dcast(Chromosome+window~Family, 
        value.var = 'window', 
        fun.aggregate = length)

# equivalent results
all.equal(as.data.frame(jack), as.data.frame(kennel))


speedcomp = microbenchmark(
  brent = lapply(split(dat, dat$Chromosome), function(y) t(sapply(split(y, list(cut(y[,"End"], as.integer(max(y$End)/as.integer(Window.Size))))), function(x) table(x$Family)))),
  karl = dat %>% group_by(Chromosome) %>% 
    mutate(position = cut(End,as.integer(max(End)/as.integer(Window.Size)))) %>% 
    ungroup() %>% 
    group_by(Chromosome, position, Family) %>% 
    summarise(abundance = n()) %>% 
    spread(Family, abundance, fill = 0) %>% 
    ungroup() %>% 
    arrange(position) %>% 
    split(f = .$Chromosome),
  jack = dat %>%
    mutate(window = round(End/100000)*100000) %>%
    group_by(Chromosome, window, Family) %>%
    summarise(obs = n()) %>%
    spread(Family, obs, fill = 0) %>% 
    ungroup(),
  kennel = dat2[, window := round(End, -5)] %>%
              dcast(Chromosome+window~Family, 
                    value.var = 'window', 
                    fun.aggregate = length),
  times = 10
)


speedcomp
# Unit: milliseconds
# expr        min         lq       mean     median         uq        max neval
# brent 2112.20703 2134.67042 2229.21085 2210.63347 2229.02289 2617.13170    10
# karl   88.08862   89.90596  125.59915   95.26605  131.61082  301.28707    10
# jack   40.07870   40.25981   49.24898   40.69559   42.94650  119.14077    10
# kennel   16.13447   16.49930   21.66882   18.35953   20.09653   54.88795    10

# jack's code is the clear winner, one order of magnitude faster compared to karl's, and again compared to brent's!
