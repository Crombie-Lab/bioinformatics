library(tidyverse)

# set the working directory to allow relative paths to project directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load blast results -YOU NEED TO HAVE THE EXACT PROJECT DIR STRUCTURE TO HAVE THIS COMMAND WORK!
blast <- read_tsv("data/raw/rhab_blast_results.tsv")

# let's glimpse the results in R
glimpse(blast)

# OK, now let's view the data in R, scroll around to see what you got.
view(blast)

#------------------------------------------------------------------------------#
# Processing the data with R
#------------------------------------------------------------------------------#
# lets clean up the subject title a bit using the stringr package from tidyverse
# first, lets get rid of the "UNVERIFIED: " section for some of these subjects

# see help for stringr::str_replace
?stringr::str_replace()

# use str_replace to clean up "UNVERIFIED"
blast2 <- blast %>%
  mutate(sub = str_replace(string = `subject title`, pattern = "UNVERIFIED: ", replacement = ""))
# view it to confirm the "UNVERIFIED: " patterns are gone. 
view(blast2)

#------------------------------------------------------------------------------#
# now lets extract the first two words of the new sub variable get the species of the subject.
# here we'll use the stringr::word function

# see help for stringr::word
?stringr::word()

# Extract the species name of the subject and assign to tax variable
blast3 <- blast2 %>%
  mutate(tax = word(string = sub, 1, 2))

# view it and see how we did? Hmm, can you see all the C. elegans, we want the full genus
view(blast3)

#------------------------------------------------------------------------------#
# now lets replace all the "C. " patterns with "Caenorhabditis" using str_replace again
blast4 <- blast3 %>%
  mutate(tax_fixed = str_replace(string = tax, pattern = "C\\.", replacement = "Caenorhabditis"))

# view it - how'd we do?
view(blast4)

#------------------------------------------------------------------------------#
# ok, good. Now we need to fix the `query_id`, I just want the sample name "S-XXXXX", not all the extra stuff
# also, I don't want the ____control samples.

# Look at the help for str_extract, which is the only new function here.
?stringr::str_extract()

# use it
blast5 <- blast4 %>%
  filter(!str_detect(string = `query id`, pattern = "control")) %>% # filter out queries with "control" in name
  mutate(query_fixed = str_extract(string = `query id`, pattern = "[^_]*"), .before = 1) # match everything before the first _

# view it
view(blast5)

#------------------------------------------------------------------------------#
# lastly, get just the best hit using distinct, since the best scoring hit is on top, this works without arrange.
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
# view it
view(blast_proc)

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
# Can you color the points by bit_score?
# Does this mean anything about our confidence in the tax_id?
# HINT - remember geom_point()?
#------------------------------------------------------------------------------#
ggplot(blast_proc) +
  geom_point(<YOU FILL IN>)

#------------------------------------------------------------------------------#
# Can you save your plots, including the scatter plot?
# HINT - remember ggsave()?
#------------------------------------------------------------------------------#
ggsave(<YOU FILL IN>)