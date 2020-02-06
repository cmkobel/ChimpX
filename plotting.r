library(tidyverse)
setwd('~/bioinformatics/chimpx/50_2_33/')

complete_gathered = readRDS("checkpoint_4_complete_gathered.rds")
#complete_gathered$gene_name = complete_gathered$gene_name %>% factor %>%  recode("control" = "ENSPTRG00000044033")
complete_gathered$gene_name = str_replace(complete_gathered$gene_name, " ", " ∩\n")
complete_gathered$gene_name = str_replace(complete_gathered$gene_name, "ENSPTRG00000044033", "control (DMD)")






# this sorted list is used for arranging the genes in the plot in a descending manner.
sorted = complete_gathered %>% group_by(gene_name, description, length) %>%
    summarise(sd = sd(CN), median = median(CN)) %>% arrange(desc(median)) %>% 
    mutate(ratio = sd/median) %>% 
    select(gene_name, median, length, sd, ratio, description)



sorted %>% ggplot(aes(sd, median, color = length)) + geom_point() # mere varians på mere kopierede gener
sorted %>% ggplot(aes(length, 1/sd, color = median)) + geom_point() # mere data på længere gener
sorted %>% ggplot(aes(median, length, color = sd)) + geom_point() # mere varians på mere kopierede gener




# distribution CN
complete_gathered %>% ggplot(aes(CN)) + geom_histogram(binwidth = .08) + 
    scale_x_continuous(breaks = 0:10)
# opdelt i køn
complete_gathered %>% ggplot(aes(CN, fill = sex)) +
    geom_density(alpha = 0.5) + # forskellige antal af hvert køn
    #geom_histogram(binwidth = 0.1, alpha = 1) +
    facet_grid(sex~.) +
    #geom_vline(xintercept = 1, alpha = 0.9, linetype = "dotted") + 
    scale_x_continuous(breaks=1:10)



# histogram of sd
# Not interesting
plot_sd = complete_gathered %>% group_by(gene_name, sex) %>% 
    summarise(sd = sd(CN))  %>% ggplot(aes(sd, color = sex)) + geom_histogram(binwidth = 0.02)
plot_sd




# main plot sorted by descending median
complete_gathered %>% #group_by(gene_name) %>% 
    ggplot(aes(x = gene_name, y = CN, color = sex)) +
    geom_hline(yintercept = 1, size = 0.2, color = "grey50") +
    geom_point(alpha = 0.5) +
    #scale_x_discrete(limits = c(arrange(sorted, desc(median))$gene_name[1:26], "control (DMD)")) + 
    scale_x_discrete(limits = c( arrange(sorted, desc(median))%>% filter(median >= 1.5) %>% pull(gene_name), "control (DMD)" )) + 
    
    scale_y_continuous(breaks = c(0:10)) +
    #theme_linedraw() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) + 
    labs(x = "gene", y = "copy number", caption = "genes are arranged in the order of descending median of CN")
ggsave("phase2_main_sd.pdf", plot = last_plot(), height = 10, width = 18)



# main plot sorted by descending sd median
complete_gathered %>% #group_by(gene_name) %>% 
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










