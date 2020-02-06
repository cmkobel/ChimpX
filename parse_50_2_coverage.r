#### This file parses the coverage files so that
#### the normalized median for each gene can be 
#### tabulated and visualized

library(tidyverse)
library(ape)
setwd('~/bioinformatics/chimpx/50_2_33/')

filenames =  c("Akwaya-Jean_cov.txt",
               "Andromeda_cov.txt",
               "Bosco_cov.txt",
               "Bwambale_cov.txt",
               "Carl_cov.txt",
               "Carolina_cov.txt",
               "Clara_cov.txt",
               "Clint_cov.txt",
               "Damian_cov.txt",
               "Dennis_cov.txt",
               "Dirk1_cov.txt",
               "Donald_cov.txt",
               "Doris_cov.txt",
               "Dylan_cov.txt",
               "Frits_cov.txt",
               "Harriet_cov.txt",
               "Jimmie_cov.txt",
               "Julie-A959_cov.txt",
               "Julie-LWC21_cov.txt",
               "Kidongo_cov.txt",
               "Koby_cov.txt",
               "Koto_cov.txt",
               "Marco_cov.txt",
               "Marlies_cov.txt",
               "Marlon_cov.txt",
               "Nakuu_cov.txt",
               "Pat_cov.txt",
               "Pearl_cov.txt",
               "Ruud_cov.txt",
               "Simliki_cov.txt",
               "Taweh_cov.txt",
               "Vaillant_cov.txt",
               "Vincent_cov.txt")
# Read the first item to get the structure
coverage_data = read_delim(paste0("data/", filenames[1]),
                  "\t", escape_double = FALSE, col_names = c('chrom', 'position', 'Akwaya-Jean'), 
                  trim_ws = TRUE) %>% select(position, `Akwaya-Jean`)

# Read the rest of the items
for (i in filenames[2:length(filenames)]) {
    pretty_name = str_split(i, "_")[[1]][1]
    i = paste0("data/", i)
    
    print(i)
    
    
    new_coverage_data = read_delim(i, 
                          "\t", escape_double = FALSE, col_names = c('chrom', 'position', pretty_name), 
                          trim_ws = TRUE)[c(2,3)]
    
    coverage_data = full_join(coverage_data, new_coverage_data, by = c("position"))
    
}
rm(new_coverage_data)

#saveRDS(coverage_data, file = "coverage_data_checkpoint_1.rds")
#coverage_data = readRDS("coverage_data_checkpoint_1.rds")


# Now we will load in the gff file, read the genes out, and merge it to the main coverage_data
annotation = read.gff("../pantro3.gff3")

# Pull out the gene names.
genes = annotation %>%
    select(-score, -phase) %>%
    filter(type == "gene") %>%
    mutate(short_name = substring(str_extract(attributes, regex("Name=\\w+")),6),
           long_name = str_sub(attributes, 9, 26)) %>%
    mutate(description = substring(str_extract(attributes, regex("description=[^\\[]+")), 13)) %>% 
    
    
    mutate(name = if_else(is.na(short_name), long_name, short_name)) %>% 
    mutate(length = abs(start - end))
    select(-short_name)



# mark the present genes at each position
position_and_names = coverage_data[,] %>% select(position) %>%
    group_by(position) %>% 
    mutate(present_gene = paste0(which(position > genes$start & position < genes$end), collapse = " ")) %>% 
    mutate(present_gene_name = 
               paste0(genes$name[which(position > genes$start & position < genes$end)], collapse = " ")) %>% 
    ungroup() %>% 
    filter(present_gene != "")

#saveRDS(position_and_names, "checkpoint_2_position_and_names.rds")
#position_and_names = readRDS("checkpoint_2_position_and_names.rds")




# Her samler jeg annotering og coverage data
full = left_join(position_and_names, coverage_data) %>% select(-position, -present_gene) %>%  group_by(present_gene_name) %>% 
    summarise_all(list(med = median), na.rm = T)
# Hvornår skal der egentlig normaliseres?


# Nu skal full transponeres, således at hver kolonne er et gen
# transpose
full_t = full %>%
    gather(key = subject, value = value, 2:ncol(full)) %>% 
    spread(key = names(full)[1], value = 'value')

# normalize ift. DMD som i pantro3 kaldes ENSPTRG...44033
full_t_norm = full_t %>%  mutate_at(
    vars(-subject), funs(. / full_t$ENSPTRG00000044033)
)

#saveRDS(full_t_norm, "checkpoint_3_full_t_norm.rds")
#full_t_norm = readRDS("checkpoint_3_full_t_norm.rds")

#remove "_med"
full_t_norm$subject = str_sub(full_t_norm$subject, 1, -5)


# load subject metadata
metadata = read_delim("~/bioinformatics/chimpx/all_38.tsv",
                      "\t", escape_double = FALSE, trim_ws = TRUE)

# change the julie names
metadata$subject = metadata$subject %>% recode(Julie_A959 = "Julie-A959",
                                               Julie_LWC21 = "Julie-LWC21")

# add metadata to full_t_norm
complete = left_join(full_t_norm, metadata) # Has genes as columns


# tidy version of complete
complete_gathered = complete %>%
    gather(key = "gene_name", value = "CN", -subject, -project, -species, -sex, -path, -PE_pairs) %>% 
    filter(subject != "Simliki") %>% 
    left_join(select(genes, gene_name = name, description, length)) # description and length will be missing from the overlapping regions


# This is the final data, which should be sufficient for all plotting.
#saveRDS(complete_gathered, "checkpoint_4_complete_gathered.rds")
#complete_gathered = readRDS("checkpoint_4_complete_gathered.rds")


# this sorted list is used for arranging the genes in the plot in a descending manner.
sorted = complete_gathered %>% group_by(gene_name, description, length) %>%
    summarise(sd = sd(CN), median = median(CN)) %>% arrange(desc(median)) %>% 
    mutate(ratio = sd/median) %>% 
    select(gene_name, median, length, sd, ratio, description)


# to do: add gene lenghs hre
write_tsv(sorted, "sorted.tsv")

sorted %>% ggplot(aes(sd, median, color = length)) + geom_point() # mere varians på mere kopierede gener
sorted %>% ggplot(aes(length, 1/sd)) + geom_point() # mere data på længere gener
sorted %>% ggplot(aes(median, length)) + geom_point() # mere varians på mere kopierede gener




# histogram of sd
plot_sd = complete_gathered %>% group_by(gene_name, sex) %>% 
    summarise(sd = sd(CN))  %>% ggplot(aes(sd, color = sex)) + geom_histogram(binwidth = 0.02)
plot_sd


# histogram of all subjects together
plot_hist = complete_gathered %>% ggplot(aes(CN)) + geom_histogram(binwidth = .25) + 
    scale_x_continuous(breaks = 0:10)
plot_hist




# main plot sorted by sd
complete_gathered %>% group_by(gene_name) %>% 
    ggplot(aes(x = gene_name, y = CN, color = sex)) +
    geom_hline(yintercept = 1, size = 0.2, color = "grey50") +
    geom_point(alpha = 0.5) +
    scale_x_discrete(limits = c(arrange(sorted, desc(median))$gene_name[1:20], "ENSPTRG00000044033")) + 
    scale_y_continuous(breaks = c(0:10)) +
    #theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) + 
    labs(x = "gene", y = "copy number", caption = "genes are arranged in the order of descending median of CN")
ggsave("phase2_main_sd.pdf", plot = last_plot(), height = 10, width = 18)



# main plot sorted by median
complete_gathered %>% group_by(gene_name) %>% 
    ggplot(aes(x = gene_name, y = CN, color = sex)) +
    geom_hline(yintercept = 1, size = 0.2, color = "grey50") +
    geom_point(alpha = 0.7) +
    scale_x_discrete(limits = c(arrange(sorted, desc(sd))$gene_name[1:20], "ENSPTRG00000044033")) + 
    scale_y_continuous(breaks = c(0:10)) +
    #theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) + 
    labs(x = "gene", y = "copy number", caption = "genes are arranged in the order of descending standard deviation of CN")
ggsave("phase2_main_median.pdf", plot = last_plot(), height = 10, width = 18)







# boxplot for each subject
complete_gathered %>% ggplot(aes(subject, CN, color = sex)) +
    geom_hline(yintercept = 1, color = "grey50") +
    geom_boxplot() + 
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 0)) + 
    scale_y_continuous(breaks = 0:10)

# main plot but with boxplot
complete_gathered %>%  ggplot(aes(gene_name, CN, color = species)) + geom_boxplot() + 
    scale_x_discrete(limits = c(sorted$gene_name[1:10], "ENSPTRG00000044033")) +
    scale_y_continuous(breaks = c(0:10)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) 
    

# number of individuals in each species
number_of_genes_present = 919 # a bit hard to calculate manually.
complete_gathered$species %>% table/number_of_genes_present



# ellioti

























