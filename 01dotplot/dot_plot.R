library(readr)
library("ape")

# Bugs:
# * genes spanning more than one window are not shown.

# Load annotation and location data from windowing.
annotation = read.gff("C:/Users/admin/Sync/Bioinformatik/Project Chimp X/00reference/pan_tro_3.0/Pan_troglodytes.Pan_tro_3.0.97.chromosome.X.gff3")
locations <- read_csv("C:/Users/admin/Sync/Bioinformatik/Project Chimp X/00reference/pan_tro_3.0/windows/locations.csv")

# Pull out the gene names.
genes = annotation %>% select(-score, -phase) %>% filter(type == "gene") %>% mutate(name = str_extract(attributes, regex("Name=\\w+"))) %>% mutate(name = substring(name, 6)) 

options(scipen=10000)

number_of_regions = 500

# Loading and plotting loop.
for (i in 1:500) {
  dots <- read_delim(paste0('0dotplots/window_', i, '.dotplot'), "\t", escape_double = FALSE, trim_ws = TRUE)
  plot_start = locations$start[which(locations$i == i)]
  plot_end = locations$end[which(locations$i == i)]
  
  # 
  # if (toggle_just == "top") {
  #   toggle_just = "bottom"
  #   toggle_nudge = 2000
  # } else if(toggle_just == "bottom") {
  #   toggle_just = "top"
  #   toggle_nudge = -2000
  # }
  # 
  
  within_range = genes %>%
    na.omit() %>% 
    filter(start >= plot_start | end <= plot_end | (start < plot_start & end > plot_end)) %>%
    mutate(start_ = ifelse(start > plot_start, start, plot_start)) %>%          # Within window start
    mutate(end_ = ifelse(end > plot_end, plot_end, end)) %>%                    # Within window end
    mutate(toggle_just = if_else(row_number() %% 2 == 0, "top", "bottom")) %>%  # alternated text adjustment
    mutate(toggle_nudge = if_else(row_number() %% 2 == 0, -2000, 2000)) %>%     #
    mutate(physical_start = if_else(strand == "+", start, end)) %>%             # physical start
    mutate(physical_end = if_else(strand == "+", end, start)) %>%               #
    mutate(toggle_nudge_x = if_else(strand == "+", 2000, -2000))                # Backwards direction of end-hook
    
  
    
  print(paste(i, plot_start, nrow(within_range)))
    

  
  # load physical positions. Necessary for gene annotation.
  if (length(dots) == 2) {
    dots[,1] = dots[, 1]+plot_start
    dots[,2] = dots[, 2]+plot_start
    
    
    # variables used for label positioning
    x_min = min(dots[,1], na.rm=T)
    x_max = max(dots[,1], na.rm=T)
    y_min = min(dots[,2], na.rm=T)
    y_max = max(dots[,2], na.rm=T)
    
    
    # plot with gene annotation
    p = dots %>% ggplot(aes(x = Pan_tro_3, y = Pan_tro_3_1)) +
      geom_path() +
      labs(title = paste0('Pan_tro_3.0 region ', i, '/', number_of_regions),
           subtitle = paste0('(', plot_start, ':', plot_end,')'))
      
      # Used for moving segments (cancelled).
    y_offset = 0
    
    # If any genes exist in the annotation for this window.
    if (nrow(within_range) > 0) {
      
      # Gene line
      p = p + geom_segment(data = within_range,
                           mapping = aes(x = start_,
                                         y = plot_start+y_offset,
                                         xend = end_,
                                         yend = plot_start+y_offset),
                           size = 1,
                           color = "cyan4",
                           alpha = 0.8)
      
      # Gene name text
      p = p + geom_text(data = within_range, mapping = aes(x = start, y = plot_start + y_offset,
                        label = name),
                        hjust = "left", vjust = within_range$toggle_just,
                        color = "cyan4",
                        nudge_y = y_offset+within_range$toggle_nudge)
                        #nudge_y = unit(1, "npc"))
      
      # Gene physical start of genes segment
      p = p + geom_point(data = within_range, mapping = aes(x = physical_start, y = plot_start + y_offset, xend = physical_start, yend = plot_start-toggle_nudge),
                         color = "cyan4",
                         size = 1,
                         shape = 15)
                         #size = 1)
      
      # Gene physical end of genes hook
      p = p + geom_segment(data = within_range, mapping = aes(x = physical_end, y = plot_start + y_offset, xend = physical_end-toggle_nudge_x, yend = plot_start-toggle_nudge),
                         color = "cyan4",
                         size = 0.5)
                         #size = 1)
      
        
                          
      }
      #p = p + xlim(x_min, x_max)
      #p = p + ylim(y_min-2000, y_max)
      p = p + xlim(plot_start, plot_end)
      p = p + ylim(plot_start-2000, plot_end)
          
      p = p  + coord_fixed(ratio = 1)
      p = p + theme_minimal()
      p = p + theme(legend.position = "none")
    
    ggsave(file=paste0('1plots\\window_', i, '.png'), plot = p,
           units = "mm", height = 300, width = 300, dpi = 150)#last_plot())
       
  }
}



# 
# 
# dots <- read_delim(paste0('0dotplots/window_1.dotplot'), "\t", escape_double = FALSE, trim_ws = TRUE)
# 
# library(tidyverse)
# dots %>% ggplot(aes(x = Pan_tro_3, y = Pan_tro_3_1)) + geom_path() + theme_bw() +
#   geom_segment(data = within_range, mapping = aes(x = start, y = 0, xend = start, yend = 20000),
#              color = 'cyan4') +
#   geom_label(data = within_range, mapping = aes(x = start, y = 0,
#                                                 label = name), alpha = 0.5, hjust = "left", color = "cyan 3")
# 
# #length(dots)
# 
# 
# dots %>% mutate(var = if_else(row_number() %% 2 == 0, 'a', 'b'))



