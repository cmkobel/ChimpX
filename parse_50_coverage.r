#### This file parses the coverage files so that
#### the normalized median for each gene can be 
#### tabulated and visualized

library(tidyverse)
library(ape)

filenames =  c("50th_Andromeda_cov.txt",
               "50th_Bosco_cov.txt",
               "50th_Bwambale_cov.txt",
               "50th_Carl_cov.txt",
               "50th_Carolina_cov.txt",
               "50th_Clara_cov.txt",
               "50th_Clint_cov.txt",
               "50th_Dennis_cov.txt",
               "50th_Dirk1_cov.txt",
               "50th_Donald_cov.txt",
               "50th_Dylan_cov.txt",
               "50th_Frits_cov.txt",
               "50th_Harriet_cov.txt",
               "50th_Jimmie_cov.txt",
               "50th_Julie_A959_cov.txt",
               "50th_Koto_cov.txt",
               "50th_Marco_cov.txt",
               "50th_Marlies_cov.txt",
               "50th_Marlon_cov.txt",
               "50th_Nakuu_cov.txt",
               "50th_Pat_cov.txt",
               "50th_Pearl_cov.txt",
               "50th_Ruud_cov.txt",
               "50th_Simliki_cov.txt",
               "50th_Taweh_cov.txt",
               "50th_Vaillant_cov.txt",
               "50th_Vincent_cov.txt")

coverage_data = read_delim("50th_Andromeda_cov.txt", 
                  "\t", escape_double = FALSE, col_names = c('chrom', 'position', 'Andromeda'), 
                  trim_ws = TRUE) %>% select(position, Andromeda)

for (i in filenames[2:length(filenames)]) {
    
    print(i)
    pretty_name = str_split(i, "_")[[1]][2]
    
    new_coverage_data = read_delim(i, 
                          "\t", escape_double = FALSE, col_names = c('chrom', 'position', pretty_name), 
                          trim_ws = TRUE)[c(2,3)]
    
    coverage_data = full_join(coverage_data, new_coverage_data, by = c("position"))
}

# saveRDS(coverage_data, file = "coverage_data_checkpoint_1.rds")
coverage_data = readRDS("data_checkpoint_1.rds")


# Now we will load in the gff file, read the genes out, and merge it to the main coverage_data
annotation = read.gff("../pantro3.gff3")

# Pull out the gene names.
genes = annotation %>%
    select(-score, -phase) %>%
    filter(type == "gene") %>%
    mutate(short_name = substring(str_extract(attributes, regex("Name=\\w+")),6),
           long_name = str_sub(attributes, 9, 26)) %>% 
    mutate(name = if_else(is.na(short_name), long_name, short_name)) %>% 
    select(-short_name, -long_name)



# mark the present genes at each position
position_and_names = coverage_data[,] %>% select(position) %>%
    group_by(position) %>% 
    mutate(present_gene = paste0(which(position > genes$start & position < genes$end), collapse = " ")) %>% 
    mutate(present_gene_name = 
               paste0(genes$name[which(position > genes$start & position < genes$end)], collapse = " ")) %>% 
    ungroup() %>% 
    filter(present_gene != "")
saveRDS(position_and_names, "checkpoint_2_position_and_names.rds")

#position_and_names %>% mutate(nchar = nchar(present_gene_name)) %>% View



full = left_join(position_and_names, coverage_data) %>% select(-position, -present_gene) %>%  group_by(present_gene_name) %>% 
    summarise_all(list(med = median), na.rm = T)
# Hvornår skal der egentlig normaliseres?


# Nu skal full transponeres, således at hver kolonne er et gen

# transpose
full_t = full %>%
    gather(key = subject, value = value, 2:ncol(full)) %>% 
    spread(key = names(full)[1], value = 'value')

# normalize
full_t_norm = full_t %>%  mutate_at(
    vars(-subject), funs(. / full_t$ENSPTRG00000044033)
)

saveRDS(full_t_norm, "checkpoint_3_full_t_norm.rds")

full_t_norm$subject = str_sub(full_t_norm$subject, 1, -5)

# Fix wrong names. More will come later
full_t_norm$subject[which(full_t_norm$subject == "Julie")] = "Julie_A959"

# load subject metadata
# load subject metadata?
metadata = read_delim("~/bioinformatics/chimpx/all_38.tsv",
                      "\t", escape_double = FALSE, trim_ws = TRUE)
complete = left_join(full_t_norm, metadata)

complete_gathered = complete %>% gather(key = "gene_name", value = "CN", -subject, -project, -species, -sex, -path, -PE_pairs)

sorted = complete_gathered %>% group_by(gene_name) %>% 
    summarise(sd = sd(CN)) %>% arrange(desc(sd))

plot_sd = complete_gathered %>% group_by(gene_name) %>% 
    summarise(sd = sd(CN))  %>% ggplot(aes(sd)) + geom_histogram()
plot_sd


complete %>% gather(key = "gene_name", value = "CN", -subject, -project, -species, -sex, -path, -PE_pairs) %>% group_by(gene_name) %>% 
    ggplot(aes(x = gene_name, y = CN, color = sex)) +
    geom_point(alpha = 0.7) +
    scale_x_discrete(limits = c(sorted$gene_name[1:100], "ENSPTRG00000044033")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))

ggsave("phase2_main.png", height = 10, width = 18)




complete %>% gather(key = "gene_name", value = "CN", -subject, -project, -species, -sex, -path, -PE_pairs) %>% group_by(gene_name) %>% 
    ggplot(aes(x = gene_name, y = sd(CN, na.rm = T), color = sex)) +
    geom_point(alpha = 0.7) +
    scale_x_discrete(limits = c(sorted$gene_name[1:100], "ENSPTRG00000044033")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))

















%>% gather(key = "subject", value = "median_coverage", ends_with("_med"))

top100 = done %>% group_by(present_gene_name) %>% summarise(sd = sd(median_coverage)) %>% arrange(desc(sd)) 

%>%
    ggplot(aes(x = present_gene_name, y = median_coverage))
