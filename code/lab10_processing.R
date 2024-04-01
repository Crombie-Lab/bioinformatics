library(tidyverse)

# set the working directory to allow relative paths to project directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

#------------------------------------------------------------------------------#
# RUN ONCE ONLY: Install vcfR package to handle vcf
#------------------------------------------------------------------------------#
install.packages("vcfR")
library(vcfR)

#------------------------------------------------------------------------------#
# Step 1: read in the vcf and extract the genotypes (gt) at the variant positions
#------------------------------------------------------------------------------#
# make sure your path to variants.vcf is correct!
vcf <- read.vcfR(file = "data/processed/variants.vcf", verbose = T)
# pull out the genotypes (gt) at variant positions from the vcf
gt <- extract.gt(vcf, element = c('GT'), as.numeric = TRUE)


#------------------------------------------------------------------------------#
# Step 2: Clean up the gt data so we can visualize the findings
#------------------------------------------------------------------------------#
# lets get the physical positions of the variants from the gt object
row.names(gt) # see the row names have what we need?
pos <- as.numeric(stringr::str_replace(row.names(gt), pattern = "KJ660346.2_", replacement = ""))

proc_gt <- as_tibble(gt) %>% # set gt to a tibble
  mutate(pos = pos, .before = 1) %>% # add the variant positions
  tidyr::pivot_longer(cols = -pos, values_to = "genotype") %>% # reshape
  dplyr::mutate(sample = case_when(name == "bam/SRR1553430.bam" ~ "EM112",
                                   name == "bam/SRR1553449.bam" ~ "G3670.1",
                                   name == "bam/SRR1553451.bam" ~ "G3676.1",
                                   name == "bam/SRR1553455.bam" ~ "G3677.1",
                                   name == "bam/SRR1553463.bam" ~ "G3682.1",
                                   name == "bam/SRR1553605.bam" ~ "NM042.1"),
                .before = name) %>% # add sample names 
  dplyr::mutate(cluster = case_when(name == "bam/SRR1553430.bam" ~ "SL3",
                                    name == "bam/SRR1553449.bam" ~ "SL1",
                                    name == "bam/SRR1553451.bam" ~ "SL1",
                                    name == "bam/SRR1553455.bam" ~ "SL2",
                                    name == "bam/SRR1553463.bam" ~ "SL2",
                                    name == "bam/SRR1553605.bam" ~ "SL3"),
                .before = name) %>% # add cluster names  from paper
  dplyr::mutate(samp_cluster = paste(sample, cluster, sep = " | "), .before = name) %>%
  dplyr::arrange(cluster, sample, pos) # arrange by cluster

# look at data
View(proc_gt)
#------------------------------------------------------------------------------#
# Step 3: Plot the variants similar to Fig4A
#------------------------------------------------------------------------------#
# lets plot the data
ggplot(proc_gt) +
  aes(x = pos, y = samp_cluster, fill = genotype) +
  geom_tile()
# hmm, why are the clusters out of order? Notice they are plotted
# alphabetically starting at the bottom of the y axis and moving up?
# we need to set the order of this variable to overide ggplot's default behavior

# we can order them properly by setting the order of samp_cluster as a factor
ggplot(proc_gt) +
  aes(x = pos, y = factor(samp_cluster, levels = rev(unique(proc_gt$samp_cluster))), fill = genotype) +
  geom_tile()

# now, notice how there are more variants than we exected given Fig4 in the paper?
# we can remove the variants that are only present in one sample, other than 10218, which is the variant
# we are def interested in. variants 1 (72bp), 4 (3388bp), 7 (10005bp) 
proc_gt2 <- proc_gt %>%
  dplyr::filter(!(pos %in% c(72, 3388, 10005)))

ggplot(proc_gt2) +
  aes(x = pos, y = factor(samp_cluster, levels = rev(unique(proc_gt$samp_cluster))), fill = genotype) +
  geom_tile()
# ok, we see only 6 variants shared across SL1, SL2, and SL3. This is not expected compared to the paper,
# Also, sample EM112 does not have a variant called at 10,218bp. We could do a bunch of invstigative work to understand
# were our analysis differed from the papers, but we see in priciple how they build figure 4A.





