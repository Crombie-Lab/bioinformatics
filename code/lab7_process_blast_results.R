library(tidyverse)

# set the working directory to allow relative paths to project directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load blast results
blast <- read_tsv("data/raw/rhab_blast.tsv2")

# lets clean up the subject title a bit using the stringr package from tidyverse
# first, lets get rid of the "UNVERIFIED: " section for some of these subjects
blast2 <- blast %>%
  mutate(sub = str_replace(string = `subject title`, pattern = "UNVERIFIED: ", replacement = ""))

# now lets just use the first two words of the sub variable to pull the taxa name of the subject
blast3 <- blast2 %>%
  mutate(tax = word(string = sub, 1, 2))
  
# now lets fix all the C. for Caenorhabditis
blast4 <- blast3 %>%
  mutate(tax_fixed = str_replace(string = tax, pattern = "C\\.", replacement = "Caenorhabditis"))

# ok, good. Now we need to fix the `query_id`, I just want the sample name "S-XXXXX", not all the extra stuff
# also, I don't want the ____control samples
blast5 <- blast4 %>%
  filter(!str_detect(string = `query id`, pattern = "control")) %>% # filter out queries with "control" in name
  mutate(query_fixed = str_extract(string = `query id`, pattern = "[^_]*"), .before = 1) # match everything before the first _

# lastly, get just the best hit using distinct, since the top hit is on top this works without arrange
blast_proc <- blast5 %>%
  group_by(query_fixed) %>%
  distinct(query_fixed, .keep_all = T) %>%
  select(query = query_fixed,
         query_len = `query length`,
         sub_len = `subject length`,
         align_len = `alignment length`,
         perc_ident = `% identity`,
         bit_score = `bit score`,
         tax_id = tax_fixed)

#------------------------------------------------------------------------------#
# lets plot the results
#------------------------------------------------------------------------------#
# basic
bar <- ggplot(blast_proc) +
  aes(x = tax_id) +
  geom_bar()
bar

# lets see that x-axis rotated and get rid of the x-axis title
bar2 <- bar +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks.x = element_blank())
bar2

#------------------------------------------------------------------------------#
# Can you make a scatter plot of perc_ident vs align_len?
# Can you color the points by bit_score
# Does this mean anything about our confidence in the tax_id
# HINT - remember geom_point()?
#------------------------------------------------------------------------------#
ggplot(blast_proc) +
  geom_point(<YOU FILL IN>)

#------------------------------------------------------------------------------#
# Can you save your plots, including the scatter plot?
# HINT - remember ggsave()?
#------------------------------------------------------------------------------#
ggsave(<YOU FILL IN>)