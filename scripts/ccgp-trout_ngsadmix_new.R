### CCGP O. mykiss
## ngsAdmix outputs --> barplots of ancestry proportions
## working with NGSadmix outputs from mega-post-bcf-exploratory-snakeflows pipeline
## output files for each K and rep are in folders called: K_#_rep_#
## cyp viii-2024


## first find out which K best describes the structure in your sample set

# specify directory with ngsAdmix output folders
ngsadmix_dir <- "ngsadmix" # thin_10_1, 1675913 SNPs

N_K <- 10    # set number of K runs you have
total_reps <- 4  # set number of reps you ran

## the following code chooses a K value based on the highest delta K

# pull all log files from the ngsAdmix output folders
log_files <- list.files(ngsadmix_dir, pattern = ".log", full.names = T, recursive=T)

# read all logs
all_logs <- lapply(1:length(log_files), FUN = function(i) readLines(log_files[i]))

# for each log file, pull out the line that starts with "best like=" from all logs, just target 'b'
# make a list of these line values
library(stringr)
bestlikes_str_list <- sapply(1:length(log_files), FUN= function(x) all_logs[[x]][which(str_sub(all_logs[[x]], 1, 1) == 'b')])

# make a dataframe with 2:N_K and N_reps to add the likelihood values
loglikes <- data.frame(K = rep(2:N_K, each=total_reps))

# then add the log likelihood value (first number in the string)
loglikes$loglike<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", bestlikes_str_list) ))

# calculate delta K and probability
# choose the K with the highest value
# example:
#   1         2         3         4         5
# Inf       Inf  41.50002  29.32132 394.18113
# in this case, it's 5 
# since this is the last K that was run, you would need to run more K
# in case the best K is greater than 5

tapply(loglikes$loglike, loglikes$K, FUN= function(x) mean(abs(x))/sd(abs(x)))
# thin_10_1 (filtered by  MAF>0.05 & F_MISSING<0.1 & FMT/DP>4 & FMT/DP<15)
#            2            3            4            5            6            7            8            9           10 
# 1.827007e+02          Inf 2.004670e+01 8.865771e+07 7.460984e+03 1.322731e+02 4.086051e+02 1.238569e+02 2.326289e+02 
chosen_k <- 5 

## next prepare your metadata
# read the sample_list.txt file from your ngsAdmix run
# the qopt file is ordered by the sample order in this file
# this code will look for sample_list.txt in the ngsadmix directory
sample_list <- read.table(file.path(ngsadmix_dir,"ccgp_samples.txt"))
colnames(sample_list) <- c("sample_ID")
# make an order column so that the sample list order can be maintained
sample_list$order_number  <- 1:nrow(sample_list)

# read in your metadata
# IMPORTANT: make sure your metadata contains ALL samples in your sample list
# if these do not match, your ngsadmix plot will not be accurate (not good!)
metadata_file <- "../metadata/CCGP-trout_sample-coordinates_correctly-geo-ranked.csv"
metadata <- read.csv(metadata_file,header=T)

# merge the sample list with metadata, keeping the order of the sample list
sample_info <- merge(sample_list,metadata,by="sample_ID")
# now reorder samples in correct order after merge
sample_info <- sample_info[order(sample_info$order_number), ]
# sanity check time: 
# the number of rows in your merged sample list should be the same as 
# the number of rows in the sample_list.txt
# if not, you didn't have ALL your sample list samples in the metadata bonk!
nrow(sample_info) # merged sample list
nrow(metadata) # metadata
# there are 4 missing sequence, check that they don't appear in sample_info
# M052862,M052864,M052867,M057325
nrow(sample_list) # sample list
# 586 samples

# set your barplot color palette here
# note that colors will be assigned to clusters randomly
trout_palette<-c("steelhead"="#BFBDD1", "redband"="#BBD3D1",
                 "little_kern"="#95BF96", "kern_river"="#467047","golden"="#6D7046", 
                 "eagle_lake"="#2f4040", "mccloud"="#546A7C", "hatchery_SilverKing"="#997a8d",
                 "hatchery_Coleman"="#B66878", "hatchery_CrystalLake"="#D3928C", 
                 "hatchery_HotCreek"="#FF9999","hatchery_MountWhitney"="#8E7068","hatchery_MountShasta"="#D991A2")

#trout_palette<-c("chinook"="salmon","steelhead_1"="#7794A1","steelhead_2"="#7794A1", "redband"="#ACD6E8",
#                 "coastal_cutthroat"="#BFBDD1","golden"="#6D7046","mccloud"="#BBD3D1",
#                 "lahontan"="#D3928C","paiute"="#8E7068")

# order the merged sample list by some variable and save that order
# for example, first order by subspecies, then by watershed, then by population name
#ord <- order(sample_info$subspecies_name,sample_info$WATERSHED,sample_info$WATER_NAME)
#ord <- order(sample_info$order_group, sample_info$rank)
ord <- order(sample_info$overall_rank)


labels <- sample_info$species_plotting_locality[ord]

# create a vector of labels to use in plot 
# the middle sample will be labelled by population name, the rest will be left blank (NA)
water_name <- sample_info$species_plotting_locality
sum_pops <- sapply(unique(water_name[ord]), function(x){sum(water_name[ord]==x)})
sum_pops.df <- as.data.frame(sum_pops)
sum_pops.df$pop_name <- row.names(sum_pops.df) 
labels <- c()
for(pop in sum_pops.df$pop_name) {
  # first label
  # labels <- c(labels, c(pop,rep(NA,sum_pops.df[sum_pops.df$pop_name == pop, 'sum_pops']-1)))
  # middle label
  if(sum_pops.df[sum_pops.df$pop_name == pop, 'sum_pops'] == 1) { labels <- c(labels,c(pop)) }
  else { labels <- c(labels, c(rep(NA,ceiling(sum_pops.df[sum_pops.df$pop_name == pop, 'sum_pops']/2)-1),pop,rep(NA,floor(sum_pops.df[sum_pops.df$pop_name == pop, 'sum_pops']/2)))) }
}
labels
cumsum_pops <- cumsum(sum_pops)
cumsum_pops

## now plot all of your ngsAdmix results
# for each K that you ran, this loop will output a pdf figure with all replicates plotted
# change the K values in this list to match the Ks you ran
for(k in c(2,3,4,5,6,7,8,9,10)) {
  
  
  # make pdf
  output_file_name <- paste0("Figures/ccgp-trout_ngsadmix_K",k,".thin10-1.pdf")
  pdf(output_file_name,width=14,height=10)
  par(mfrow=c(total_reps+2,1), mar=c(1,4,1,1))
  
  # loop through replicates
  for(rep in c(1:total_reps)) {
    # read in the qopt ancestry proportions
    q <- read.table(file.path(ngsadmix_dir,paste0("K_",k,"_rep_",rep),"output.qopt"))
    q$sample_ID <- sample_list$sample_ID
    q$replicate <- rep
    #if(rep==1) { q_all_reps = q }
    #else { q_all_reps <- rbind(q_all_reps,q)
    
    # for visual purposes: if on the last replicate, plot the population name
    if(rep==total_reps) {
      barplot(t(q)[,ord], # order the q values by the saved order
              col=trout_palette,
              #names=sample_info$plot_water_name[ord], # order the population names by the saved order
              names=labels, # apply the labels defined before
              las=2,
              space=0,
              border=NA,
              cex.names=0.8,
              xlab="",
              ylab=paste0("admixture proportions for K=",k))
      # add line between populations
      abline(v=cumsum(sapply(unique(water_name[ord]), function(x){sum(water_name[ord]==x)})), col = "white", lwd = 0.5)
    }
    # for visual purposes: if on any but the last replicate, plot without a name
    else {
      barplot(t(q)[,ord],
              col=trout_palette,
              las=2,
              space=0,
              border=NA,
              xlab="",
              ylab=paste0("admixture proportions for K=",k))
      # add line between populations
      abline(v=cumsum(sapply(unique(water_name[ord]), function(x){sum(water_name[ord]==x)})), col = "white", lwd = 0.5)
    }
  }
  dev.off()
}



### Final figure
## individual K plots so that cluster colors match the group
## rep chosen based on clustering with most replicates
pdf(paste0("Figures/FigX_ngsadmix.pdf"),width=14,height=8)
par(mfrow=c(6,1), mar=c(1,6,1,1))

# K=4 (rep 4)
trout_palette_4<-c( "extra2" = "#FFB26B",
                    "hatchery"="#997a8d", 
                    "steelhead"="#FFD56F",
                    "extra1" = "#FF7B54")
q <- read.table(file.path(ngsadmix_dir,"K_4_rep_4","output.qopt"))
barplot(t(q)[,ord],
        col=trout_palette_4,
        las=2,
        space=0,
        border=NA,
        xlab="",
        ylab=paste0("admixture proportions (K=4)"))
# add line between populations
abline(v=cumsum(sapply(unique(water_name[ord]), function(x){sum(water_name[ord]==x)})), col = "white", lwd = 0.5)

# K=3 (rep 4)
trout_palette_3<-c( "extra1" = "#FF7B54",
                    "steelhead"="#FFD56F",
                    "hatchery"="#997a8d")
q <- read.table(file.path(ngsadmix_dir,"K_3_rep_4","output.qopt"))
barplot(t(q)[,ord],
        col=trout_palette_3,
        las=2,
        space=0,
        border=NA,
        xlab="",
        ylab=paste0("admixture proportions (K=3)"))
# add line between populations
abline(v=cumsum(sapply(unique(water_name[ord]), function(x){sum(water_name[ord]==x)})), col = "white", lwd = 0.5)


# K=2 (rep 1)
trout_palette_2<-c("steelhead"="#FFD56F",
                   "hatchery"="#997a8d")
q <- read.table(file.path(ngsadmix_dir,"K_2_rep_1","output.qopt"))
barplot(t(q)[,ord], 
        col=trout_palette_2,
        names=labels, 
        las=2,
        space=0,
        border=NA,
        cex.names=1,
        xlab="",
        ylab=paste0("admixture proportions (K=2)"))


# add line between populations
abline(v=cumsum(sapply(unique(water_name[ord]), function(x){sum(water_name[ord]==x)})), col = "white", lwd = 0.5)

dev.off()
