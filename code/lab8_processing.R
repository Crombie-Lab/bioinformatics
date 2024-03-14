library(tidyverse)

# set the working directory to allow relative paths to project directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#------------------------------------------------------------------------------#
# RUN ONCE ONLY: Setup the required packages using the BiocManager package manager
#------------------------------------------------------------------------------#
# install the BiocManager packge
install.packages("BiocManager")

# use BiocManager to install ggtree
# If you get a prompt (Update all/some/none? [a/s/n]: ), enter "a"
BiocManager::install("ggtree")

# use BiocManager to install ggmsa
# If you get a prompt (Update all/some/none? [a/s/n]: ), enter "a"
BiocManager::install("ggmsa")

# use BiocManager to install treeio
# If you get a prompt (Update all/some/none? [a/s/n]: ), enter "a"
BiocManager::install("treeio")

# use BiocManager to install adephylo
# If you get a prompt (Update all/some/none? [a/s/n]: ), enter "a"
BiocManager::install("adephylo")


#------------------------------------------------------------------------------#
# Part 1: Examine the ebola MSA file ebov.mafft.fasta 
#------------------------------------------------------------------------------#
# load the ggmsa package
# The package documentation is here: https://yulab-smu.top/ggmsa/index.html
library(ggmsa)

# load the .fasta file containing the nucleotide msa - this can take 5 minutes b/c
# the MSA is enormous. To help save time, lets just look at the 100 bp region
# in the middle of the nucleoprotein gene which spans from 470 bp and runs to 2689 bp
# in the ebola genome. Let's assume that the MSA positions are close to the
# physical positions of the genome: Middle of the nucleoprotein (1480bp - 1580bp).
msa <- ggmsa::ggmsa(msa = "data/raw/ebov.mafft.fasta", 1480, 1580,
                    color = "Chemistry_NT", font = "DroidSansMono",
                    char_width = 0.5, seq_name = TRUE)

# now add some fancy features with a seqLogo plot and a consensus bar to look purdy
msa_add <- msa + geom_seqlogo() + geom_msaBar()

# save the msa views to plots so we can see what ggmsa does for us. This can take a 
# few minutes too.
ggsave(msa, filename = "plots/middle_NP.png", height = 14, width = 14)
ggsave(msa_add, filename = "plots/middle_NP_fancy.png", height = 14, width = 14)

#------------------------------------------------------------------------------#
# Part 2: Read in the unrooted tree file and explore ggtree plotting syntax
#------------------------------------------------------------------------------#
# load the required packges
library(ggtree)
library(tidytree)

# use the read.tree function from treeio package
unrooted_tree <- treeio::read.tree("data/raw/unrooted.raxml.bestTree")

# plot the tree using ggtree
# The ggtree book is here: https://yulab-smu.top/treedata-book/
unrooted_plot <- ggtree::ggtree(unrooted_tree, layout="equal_angle")
# View it as is - This should look a lot like Figure 2A but rotated differently
unrooted_plot

# lets add a scale to show the branch length divergence
unrooted_plot2 <- ggtree::ggtree(unrooted_tree, layout="equal_angle") +
                      geom_treescale(x = 0.005, y = 0.02, width = 0.001, offset = -0.0015) # add a scale bar with specific position
# view this plot
unrooted_plot2

# now we can add in some custom clade labels so we can see the years like in Fig 2 
# of the paper. For this we need to see the internal node names within the tree object.
# Get the internal node numbers for the most recent common ancestor of each viral year.
tidy_unrooted_tree <- tidytree::as_tibble(unrooted_tree) %>%
  tidyr::separate(label, into = c("virus", "year", "id"), remove = F) %>%
  dplyr::filter(!is.na(virus)) %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(mrca_node = min(parent) + 1)

# We can use the mrca node from above to label clades. The 2002 isolate is single so
# I can't label it the same way, but you get the idea
unrooted_plot3 <- unrooted_plot2 +
  geom_cladelab(node = 114, label = "76'77'", offset = 0.0001, fontsize = 3) +
  geom_cladelab(node = 120, label = "95'", offset = 0.0002, fontsize = 3) +
  geom_cladelab(node = 117, label = "94'\n96'", offset = -0.0002, fontsize = 3) +
  geom_cladelab(node = 123, label = "07'08'", offset = 0.0001, fontsize = 3) +
  geom_cladelab(node = 103, label = "14'", offset = 0.0003, fontsize = 3) 
unrooted_plot3

# save the tree plot
ggsave(plot = unrooted_plot3, filename = "plots/fig2A.png", width = 5, height = 5)

#------------------------------------------------------------------------------#
# Part 3: Read in the 2014 and 1976 rooted trees and explore the data
#------------------------------------------------------------------------------#
library(ggtree)
library(tidytree)

# use the read.tree function from treeio package
rooted_2014_tree <- treeio::read.tree("data/raw/root2014.raxml.bestTree")

# plot the tree using ggtree
# The ggtree book is here: https://yulab-smu.top/treedata-book/
# here we'll use the default 
rooted_2014_plot <- ggtree::ggtree(rooted_2014_tree, layout = "rectangular")
# view it
rooted_2014_plot

# OK, let's read in the 1976 rooted tree, which is the rooting they choose.
# use the read.tree function from treeio package
rooted_1976_tree <- treeio::read.tree("data/raw/root1976.raxml.bestTree")
# plot the tree using ggtree
rooted_1976_plot <- ggtree::ggtree(rooted_1976_tree, layout = "rectangular")
  
# lets try collapsing some of the 2014 sequences based on the internal node
rooted_1976_plot2 <- ggtree::ggtree(rooted_1976_tree, layout = "rectangular") %>%
  ggtree::collapse(node = 119) + # collapse some of the 2014 samples
  #geom_tiplab() +
  ggplot2::xlim(0, 0.06) + # get some extra space on the x-axis of the plot to not crop labels
  geom_treescale(x = 0.005, y = 27.5, width = 0.001, offset = 0.5) +  # add a scale bar with specific position
  geom_strip("EBOV_2007_HQ613403", "EBOV_2007_KC242789", label = "DRC, 07'-08'",
             offset.text= 0.001, fontsize = 3, color = "salmon") +
  geom_strip("EBOV_2014_G3687", "EBOV_2014_KJ660347", label = "West Africa 14'",
             offset.text= 0.001, fontsize = 3, color = "darkgreen") +
  geom_strip("EBOV_2002_KC242800", "EBOV_2002_KC242800", label = "Gabon, 02'",
             offset.text= 0.001, fontsize = 3, color = "orange") +
  geom_strip("EBOV_1995_AY354458", "EBOV_1995_KC242796", label = "Zaire (DRC), 95'",
             offset.text= 0.001, fontsize = 3, color = "powderblue") +
  geom_strip("EBOV_1996_KC242793", "EBOV_1996_KC242794", label = "Gabon, 94'-96'",
             offset.text= 0.001, fontsize = 3, color = "powderblue") +
  geom_strip("EBOV_1977_KC242791", "EBOV_1976_KC242801", label = "Zaire (DRC), 76'-77'",
             offset.text= 0.001, fontsize = 3, color = "royalblue")
  
# View it  
rooted_1976_plot2

# save the nice tree
ggsave(plot = rooted_1976_plot2, filename = "plots/fig2C.png", width = 5, height = 5) 

#------------------------------------------------------------------------------#
# Part 4: See if we can plot the correlations from root to tip for the trees
#------------------------------------------------------------------------------#
# load the adephylo package 
library(adephylo)

# calculate the distances to output a named vector with distances from the root
# this outputs a named vector
root14_dist <- adephylo::distRoot(x = rooted_2014_tree, tips = "all", method = "patristic")
root14_dist
# we can make a dataframe out of this using the tibble function
dist14_df <- tibble::tibble(sample = names(root14_dist),
                            root_tip_dist = as.numeric(root14_dist),
                            root = "2014 root") %>%
  tidyr::separate(sample, into = c("virus", "year", "id"), remove = F)

# do the same for the 1976 rooted tree
root76_dist <- adephylo::distRoot(x = rooted_1976_tree, tips = "all", method = "patristic")
root76_dist
# we can make a dataframe out of this using the tibble function
dist76_df <- tibble::tibble(sample = names(root76_dist),
                            root_tip_dist = as.numeric(root76_dist),
                            root = "1976 root") %>%
  tidyr::separate(sample, into = c("virus", "year", "id"), remove = F)

# now we need to combine the two dataframes to get it all in one
dist_df <- bind_rows(dist14_df, dist76_df)

# Plot the results! No R^2 here, but it's good enough
dist_plot <- ggplot(dist_df) +
  aes(x = as.numeric(year), y = root_tip_dist) +
  geom_point() +
  facet_wrap(~root, ncol = 1, scales = "free") +
  geom_smooth(method = "lm") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = "Year", y = "Root-to-tip divergence")

# look at it
dist_plot

# save it
ggsave(plot = dist_plot, filename = "plots/fig2B.png", width = 5, height = 5)

#------------------------------------------------------------------------------#
# Part 5: Combine the plots using cowplot
#------------------------------------------------------------------------------#
# try putting it all together.
fig2 <- cowplot::plot_grid(unrooted_plot3, dist_plot, rooted_1976_plot2, ncol = 3, align = "bt", labels = c("A", "B", "C"))
fig2

# save it
ggsave(plot = fig2, filename = "plots/fig2.png", width = 7.5, height = 3.75)

