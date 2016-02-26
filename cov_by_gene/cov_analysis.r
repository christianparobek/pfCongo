## An R script to determine whether coverage over HRP2 is less than expected
## Given the coverage over the rest of the genes in the genome
## Christian P
## 26 February 2015
##
## USAGE: Rscript cov_analysis.r <in.res>


#############################
######### PARSE CMD #########
#############################

args = commandArgs(trailingOnly=TRUE)


#############################
######### READ DATA #########
#############################

data <- read.table(args[1])
  # a table of bedtools coverage -d output

names(data) <- c("chr", "start", "end", "count", "depth")
  # name the columns appropriately

data <- within(data, id <- paste(chr, start, end, sep='.'))
  # create a unique ID column, which is the chr start and end for each gene


#############################
######### AGGREGATE #########
#############################

agg <- aggregate(depth ~ id, data=data, FUN=sum)
  # sum depth column as aggregate by id column

geneloc <- matrix(unlist(strsplit(agg$id, split = "[.]")), ncol = 3, byrow = TRUE)
  # during aggregate(), we lose the useful data from chr, start, and end columns
  # so get split those out into a matrix here

new <- cbind(agg, geneloc)
  # make a new data table with the aggregated info
  # has chr, start, and end info

names(new) <- c("id", "total_ct", "chr", "start", "end")
  # give this new table some appropriate names


#############################
######### NORMALIZE #########
#############################

new$lengths <- as.numeric(as.character(new$end)) - as.numeric(as.character(new$start))
  # get lengths of each gene

new$norm <- new$total_ct / new$lengths
  # normalize coverage by gene length


#############################
######### PLOT COV ##########
#############################

hist(new$norm, breaks = 100000, xlim = c(0,2))
plot(density(new$norm), xlim = c(0,2))
  # plot coverage


#############################
######### PLOT COV ##########
#############################

hrp2 <- new[new$id == "Pf3D7_08_v3.1374236.1375299",]$norm # hrp2
hrp3 <- new[new$id == "Pf3D7_13_v3.2840727.2841703",]$norm # hrp3

print("hrp2")
sum(hrp2 < new$norm)/length(new$norm)
sum(hrp2 <= new$norm)/length(new$norm)


print("hrp3")
sum(hrp3 < new$norm)/length(new$norm)
sum(hrp3 <= new$norm)/length(new$norm)


