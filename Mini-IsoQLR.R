# Author: Gonzalo Nunez Moreno
# Date: 05/07/2021

library(ggplot2)
#library(ggrepel)
library(plyr)
library(optparse)
library(cowplot)



#############
# Arguments # 
#############
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input gff3 file (output from GMAP)", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-l", "--logfile"), type="character", default=NULL, 
              help="error file from GMAP", metavar="character"),
  make_option(c("-r", "--runame"), type="character", default=NULL, 
              help="run name that will be used as a file prefix and in figure titles", metavar="character"),
  
  make_option(c("-k", "--known_sites"), type="character", default=NULL, 
              help="4 column, tab delimiter file with known splice sites (no header). Columns: type of splece site (\"start\"/\"end\" of the exon), splice site position, minor position, mayor position. Minor and mayor positions are used to asign all the splice sites to the consensus splice site", metavar="character"),
  make_option(c("-s", "--filter_sites"), type="character", default=NULL, 
              help="2 column, tab delimiter file with known splice sites (no header) that reads must have to be included in the analysis. This option can be used to filter PCR artifacts by taking into account only reads with both ends (where primers hibridate).  Columns: type of splece site (\"start\"/\"end\" of the exon) and splice site position", metavar="character"),
  
  make_option(c("-b", "--beginning"), type="integer", default=0, 
              help="beginning position of the segment of study (trimming)", metavar="integer"),
  make_option(c("-f", "--final"), type="integer", default=NULL, 
              help="final position of the segment of study (trimming)", metavar="integer"),
  
  make_option(c("-t", "--threshold"), type="double", default=5, 
              help="threshold (0-100%) used to filter the breakpoints present in mode than x % of the reads [default= %default]", metavar="integer"),
  make_option(c("-p", "--padding"), type="integer", default=5, 
              help="number of bases to each side from a defined break point to consider a read as part of that group [default= %default]", metavar="integer"),
  make_option(c("-a", "--abundance"), type="double", default=5, 
              help="only isoforms with a percentage equal or higher will be displayed on the combined plot [default= %default]", metavar="integer")
  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


gff3_path = opt$input
plot_path = opt$outputdir
log_map_path = opt$logfile
run_name = opt$runame
known_sites_file = opt$known_sites
filter_sites_file = opt$filter_sites
inicio = opt$beginning
final = opt$final
breakpoint_freq_threshold = opt$threshold
breakpoint_padding = opt$padding
abundance = opt$abundance

# gff3_path = "/home/gonzalo/tblab/home/gonzalo/pax6/post_trimm_multiplex1/mapped/cluster_cons_BARCODE02.gff3"
# plot_path = "/home/gonzalo/tblab/home/gonzalo/pax6/post_trimm_multiplex1/plots/plot_BARCODE02/"
# log_map_path = "/home/gonzalo/tblab/home/gonzalo/pax6/post_trimm_multiplex1/mapped/log_BARCODE02.err"
# run_name = "Multiplex1_BARCODE02"
# known_sites_file = NULL
# filter_sites_file = NULL
# inicio = 1200
# final = 3500
# breakpoint_freq_threshold = 5
# breakpoint_padding = 5
# abundance = 5

# gff3_path = "/home/gonzalo/tblab/home/gonzalo/pax6/post_trimm_multiplex2_2/mapped/cluster_cons_BARCODE08.gff3"
# plot_path = "/home/gonzalo/tblab/home/gonzalo/pax6/post_trimm_multiplex2_2/plots/plot_BARCODE08/"
# log_map_path = "/home/gonzalo/tblab/home/gonzalo/pax6/post_trimm_multiplex2_2/mapped/log_BARCODE08.err"
# run_name = "Multiplex1_BARCODE02"
# known_sites_file = "/home/gonzalo/tblab/home/gonzalo/pax6/211115_minigenes_multiplex3_2/knownsites.txt"
# filter_sites_file = "/home/gonzalo/tblab/home/gonzalo/pax6/211115_minigenes_multiplex3_2/filtersites.txt"
# inicio = 0
# final = NULL
# breakpoint_freq_threshold = 5
# breakpoint_padding = 5
# abundance = 5


################
# Data loading #
################

gff3 = read.delim(gff3_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
gff3 = gff3[gff3$V3 == "exon",c("V4","V5","V9", "V6")]
colnames(gff3) = c("start", "end", "id", "score")
gff3$score = as.numeric(gff3$score)

gff3$id=gsub("^.*Name=", "", gff3$id, perl = TRUE)
gff3$id=gsub(";.*$", "", gff3$id, perl = TRUE)


num_reads_initial = length(unique(gff3$id))
# print(paste("Reads iniciales: ", num_reads_initial, sep = ""))


############
# Trimming #
############
# Trimming of unwanted exons
if (is.null(final)) {
  final = max(gff3$end) + 1
}

gff3_filtered = gff3[gff3$end > inicio & gff3$start < final, ]

num_reads_post_trimming = length(unique(gff3_filtered$id))
# print(paste("Reads post-trimming: ", num_reads_post_trimming, sep = ""))


##########################
# Break point definition #
##########################

# Break points frequency
start_count = table(gff3_filtered$start)
end_count = table(gff3_filtered$end)


# Most frequent breakpoints
start_count = start_count[start_count>(num_reads_post_trimming*breakpoint_freq_threshold/100)]
end_count = end_count[end_count>num_reads_post_trimming*breakpoint_freq_threshold/100]

start_position = data.frame(breakp = as.numeric(names(start_count)), freq = as.numeric(start_count), stringsAsFactors = F)
end_position = data.frame(breakp = as.numeric(names(end_count)), freq = as.numeric(end_count), stringsAsFactors = F)



# Break point boundaries
start_position$left = start_position$breakp - breakpoint_padding
start_position$right = start_position$breakp + breakpoint_padding
end_position$left = end_position$breakp - breakpoint_padding
end_position$right = end_position$breakp + breakpoint_padding



for (i in 1:(nrow(start_position)-1)){
  if (start_position[i, "right"] > start_position[i+1, "left"]){
    media_pos = (start_position[i, "right"] + start_position[i+1, "left"])/2
    if (media_pos%%1==0) { # Check if the number is decimal or interger
      start_position[i, "right"] = media_pos - 1
      start_position[i+1, "left"] = media_pos + 1
    }else{
      start_position[i, "right"] = floor(media_pos)
      start_position[i+1, "left"] = ceiling(media_pos)
    }
  }
}



for (i in 1:(nrow(end_position)-1)){
  if (end_position[i, "right"] > end_position[i+1, "left"]){
    media_pos = (end_position[i, "right"] + end_position[i+1, "left"])/2
    if (media_pos%%1==0) { # Check if the number is decimal or interger
      end_position[i, "right"] = media_pos - 1
      end_position[i+1, "left"] = media_pos + 1
    }else{
      end_position[i, "right"] = floor(media_pos)
      end_position[i+1, "left"] = ceiling(media_pos)
    }
  }
}





# Known break points
if (!is.null(known_sites_file)){
  known_sites = read.delim(known_sites_file, header = F, stringsAsFactors = F)
  colnames(known_sites) = c("tipo", "site", "lower", "upper")
  known_sites$tipo = tolower(known_sites$tipo)
  
  for (i in 1:nrow(known_sites)){
    if (known_sites[i, "tipo"] == "start"){
      start_position = start_position[!(start_position$breakp >= known_sites[i, "lower"] & start_position$breakp <= known_sites[i, "upper"]), ]
      start_position[start_position$breakp < known_sites[i, "lower"] & start_position$right >= known_sites[i, "lower"], "right"] = known_sites[i, "lower"] - 1
      start_position[start_position$breakp > known_sites[i, "lower"] & start_position$left <= known_sites[i, "upper"], "left"] = known_sites[i, "upper"] + 1
      start_position = rbind(start_position, c(known_sites[i, "site"],
                                               table(gff3_filtered$start)[as.character(known_sites[i, "site"])],
                                               known_sites[i, "lower"],
                                               known_sites[i, "upper"]))
    } else if (known_sites[i, "tipo"] %in% c("end", "stop")) {
      end_position = end_position[!(end_position$breakp >= known_sites[i, "lower"] & end_position$breakp <= known_sites[i, "upper"]), ]
      end_position[end_position$breakp < known_sites[i, "lower"] & end_position$right >= known_sites[i, "lower"], "right"] = known_sites[i, "lower"] - 1
      end_position[end_position$breakp > known_sites[i, "lower"] & end_position$left <= known_sites[i, "upper"], "left"] = known_sites[i, "upper"] + 1
      end_position = rbind(end_position, c(known_sites[i, "site"],
                                           table(gff3_filtered$end)[as.character(known_sites[i, "site"])],
                                           known_sites[i, "lower"],
                                           known_sites[i, "upper"]))
    }
  }
}
  

  





# Break point asignment
gff3_filtered$start_tag = NA
for (i in 1:nrow(start_position)){
  gff3_filtered[gff3_filtered$start >= start_position[i, "left"]  & gff3_filtered$start <= start_position[i, "right"], "start_tag"] = start_position[i, "breakp"]
}


gff3_filtered$end_tag = NA
for (i in 1:nrow(end_position)){
  gff3_filtered[gff3_filtered$end >= end_position[i, "left"]  & gff3_filtered$end <= end_position[i, "right"], "end_tag"] = end_position[i, "breakp"]
}





###################
# Exon difinition #
###################


# Read filtering if start/end do not match with with the estimated break points
exon_id_filter = gff3_filtered[rowSums(is.na(gff3_filtered)) > 0, "id"]
gff3_filtered_by_exon = gff3_filtered[!gff3_filtered$id %in% exon_id_filter,]

num_reads_post_exon_filtering = length(unique(gff3_filtered_by_exon$id))
# print(paste("Reads post exon filtering: ", num_reads_post_exon_filtering, sep = ""))


# Read filtering selected by argument
if (!is.null(filter_sites_file)){
  filter_sites = read.delim(filter_sites_file, header = F, stringsAsFactors = F)
  colnames(filter_sites) = c("tipo", "site")
  filter_sites$tipo = tolower(filter_sites$tipo)
  
  #filter_start = filter_sites[filter_sites$tipo == "start", "site"]
  #filter_end = filter_sites[filter_sites$tipo %in% c("end", "stop"), "site"]
  
  read_id_filter = c()
  for (i in 1:nrow(filter_sites)){
    if (filter_sites[i, "tipo"] == "start"){
      read_id_filter = c(read_id_filter, gff3_filtered_by_exon[gff3_filtered_by_exon$start_tag == filter_sites[i, "site"],"id"])
    } else if (filter_sites[i, "tipo"] %in% c("end", "stop")) {
      read_id_filter = c(read_id_filter, gff3_filtered_by_exon[gff3_filtered_by_exon$end_tag == filter_sites[i, "site"],"id"])
    }
  }
  
  final_id_filter = names(table(read_id_filter))[table(read_id_filter) == nrow(filter_sites)]
  final_id_filter_out = unique(gff3_filtered_by_exon[!gff3_filtered_by_exon$id %in% final_id_filter,"id"])
  gff3_filtered_by_exon = gff3_filtered_by_exon[gff3_filtered_by_exon$id %in% final_id_filter,]
  num_reads_post_site_filtering = length(unique(gff3_filtered_by_exon$id))
}



# Exon difinition
gff3_filtered_by_exon$coordinates = paste(gff3_filtered_by_exon$start_tag, gff3_filtered_by_exon$end_tag, sep = "_")



######################
# Isoform definition #
######################
# Isoform definition
read_list = split(gff3_filtered_by_exon, gff3_filtered_by_exon$id)
read_isoform = unlist(lapply(read_list, function(x) paste(sort(x$coordinates), collapse = "-")))

# frequency calculation
isoform_frequencies = data.frame(sort(table(read_isoform), decreasing = TRUE))
isoform_frequencies$read_isoform = as.character(isoform_frequencies$read_isoform)
isoform_frequencies$perc = isoform_frequencies$Freq * 100 / sum(isoform_frequencies$Freq)
isoform_frequencies$perc = unlist(lapply(isoform_frequencies$perc, function(x) round(x, 1)))









#########
# Plots #
#########

##################################
# Alignment plot of the isoforms #
##################################
isoforms = as.character(isoform_frequencies$read_isoform)
isoform_num = 1
x_pos = c()
tipo = c()
exon = c()
y_pos = c()
isoform_id_all = c()
perc = c()
star_stop  =c()

for (i in 1:nrow(isoform_frequencies)) {
  
  #isoform_id = paste("Iso", isoform_num, ": ", isoform_frequencies[i,2] ," reads (", isoform_frequencies[i,3], "%)", sep = "")
  isoform_id = paste("Iso", isoform_num, " (", isoform_frequencies[i,3], "%)", sep = "")
  isoform_id_all = c(isoform_id_all, isoform_id)
  isoform_num = isoform_num + 1
  
  positions = strsplit(isoform_frequencies[i,1], split = "-")[[1]]
  exon_num = 1
  
  for (j in positions){

    start = strsplit(j, split = "_")[[1]][1]
    end = strsplit(j, split = "_")[[1]][2]
    
    x_pos = c(x_pos, start, end)
    tipo = c(tipo, "start", "end")
    exon = c(exon, exon_num, exon_num)
    y_pos = c(y_pos, isoform_id, isoform_id)
    perc = c(perc, isoform_frequencies[i,3], isoform_frequencies[i,3])
    star_stop = c(star_stop, paste(start, end, sep = "-"), paste(start, end, sep = "-"))
    
    exon_num = exon_num + 1
  }
}

x_pos = as.integer(x_pos)
df_pos = data.frame(x_pos, tipo, exon, y_pos, perc, star_stop)

df_pos$paste_pair = paste(df_pos$exon, df_pos$y_pos, sep = "_")
df_pos$y_pos = ordered(df_pos$y_pos, c(rev(isoform_id_all),"Break points"))

# Plot of all isoforms
num_isoforms = max(c(length(unique(df_pos$y_pos)), length(unique(df_pos$star_stop)) + 1), 4)

p0 = ggplot(df_pos, aes(x = x_pos, y = y_pos, group = paste_pair, color = star_stop)) + 
  geom_line(size = 2) +
  theme_bw() +
  ylab("Isoform ID\n(% of reads)") +
  xlab("Position") +
  labs(color = "Exon coordinates") +
  scale_x_continuous(position = "top") +
  ggtitle(run_name)
#p0

ggsave(path = plot_path, filename = paste(run_name, ".all_isoform_information.pdf", sep = ""), p0, device="pdf", width = 10,
       height = 1 + num_isoforms/4, limitsize = F)
ggsave(path = plot_path, filename = paste(run_name, ".all_isoform_information.jpeg", sep = ""), p0, device="jpeg", width = 10,
       height = min(1 + num_isoforms/4, 100), limitsize = F)





reference = data.frame(x_pos = unique(df_pos[df_pos$perc >= abundance, "x_pos"]),
                       perc = 100,
                       y_pos = factor("Break points", levels = c(rev(isoform_id_all), "Break points")),
                       paste_pair = seq(1,length(unique(df_pos[df_pos$perc >= abundance, "x_pos"])),1))

df_plot = rbind.fill(df_pos, reference)
df_plot$y_pos = ordered(df_plot$y_pos, c(rev(isoform_id_all),"Break points"))





#####################
# BREAK POINT Plots #
#####################

# Break point superplot
start_df = data.frame(table(gff3[, "start"]))
start_df$tipo = "start"
end_df = data.frame(table(gff3[, "end"]))
end_df$tipo = "end"

df_superplot = rbind(start_df, end_df)

df_superplot$Var1 = as.numeric(as.character(df_superplot$Var1))
df_superplot$Freq = as.numeric(as.character(df_superplot$Freq))
df_superplot$perc = df_superplot$Freq/num_reads_initial*100

p3 = ggplot(df_superplot, aes(x=Var1, y=perc, fill = tipo)) +
  #geom_density(adjust = 1/5, alpha=0.5, position="identity") +
  #geom_histogram(aes(y=..count../num_reads_initial*100),alpha=0.7, position="identity", bins = round((max(gff3$start)-min(gff3$start)))) +
  geom_col(position="identity", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, max(df_superplot$Var1), 10)) +
  #scale_y_continuous(trans = "log10") +
  ylab("Percentage of reads") +
  xlab("Break point position") +
  labs(fill = "Break point\ntype") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, size = 2))
#p3

ggsave(path = plot_path, filename = paste(run_name, ".breakpoint_superplot.pdf", sep = ""), p3, device="pdf", width = 100,
       height = 5, limitsize = F)
ggsave(path = plot_path, filename = paste(run_name, ".breakpoint_superplot.jpeg", sep = ""), p3, device="jpeg", width = 100,
       height = 5, limitsize = F)





# Break point plot
df_plot = data.frame(posicion = c(gff3$start, gff3$end),
                     tipo = rep(c("start", "end"), each = nrow(gff3)),
                     read_id = gff3$id,
                     stringsAsFactors = FALSE)

p4 = ggplot(df_plot, aes(x=posicion, fill = tipo)) +
  geom_histogram(aes(y=..count../num_reads_post_trimming*100),alpha=0.5, 
                 position="identity", bins = 500) +  
  scale_x_continuous(breaks = seq(inicio,final,round((final-inicio)/12, 0)), limits = c(inicio-10,final+10)) +
  ylab("Percentage of reads (after trimming)") +
  xlab("Break point position") +
  labs(fill = "Break point\ntype") +
  theme_bw()
#p4

ggsave(path = plot_path, filename = paste(run_name, ".breakpoints.pdf", sep = ""), p4, device="pdf", width = 10,
       height = 5, limitsize = F)
ggsave(path = plot_path, filename = paste(run_name, ".breakpoints.jpeg", sep = ""), p4, device="jpeg", width = 10,
       height = 5, limitsize = F)







#################
# Combined plot #
#################

isoform_frequencies_plot = isoform_frequencies[isoform_frequencies$perc>=1,]
isoform_frequencies_plot$iso = paste("Iso", row.names(isoform_frequencies_plot), sep = "")
isoform_frequencies_plot  = rbind(isoform_frequencies_plot, c("-", 0, 100-sum(isoform_frequencies_plot$perc), "Other"))
isoform_frequencies_plot$iso = factor(isoform_frequencies_plot$iso, levels = c("Other", rev(isoform_id_all)))
isoform_frequencies_plot$perc = as.numeric(isoform_frequencies_plot$perc)


p11 = ggplot(df_pos[df_pos$perc >= abundance,], aes(x = x_pos, y = y_pos, color = star_stop)) + 
  geom_line(aes(group = paste_pair), size = 2) +
  ylab("Isoform ID\n(% of reads)") +
  xlab("Position") +
  labs(color = "Exon coordinates") +
  labs(color = "Exon coordinates") +
  geom_hline(yintercept="Break points", alpha = 0.10, 
             color = "red", size=20) +
  #scale_y_discrete(position = "right") +  
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_x_continuous(breaks = seq(inicio,final,round((final-inicio)/12, 0)), limits = c(inicio-10,final+10)) +
  coord_cartesian(c(inicio,final)) +
  theme(plot.margin = unit(c(0.7, 0.7, 0, 0.7), "cm")) +
  theme(plot.background = element_rect(fill = "white", colour = "white")) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(run_name)
#p11



p44 = ggplot(df_plot, aes(x=posicion, fill = tipo)) +
  geom_histogram(aes(y=..count../num_reads_post_trimming*100),alpha=0.5, 
                 position="identity", bins = 500) +
  scale_x_continuous(breaks = seq(inicio,final,round((final-inicio)/12, 0)), limits = c(inicio-10,final+10)) +
  ylab("% of reads") +
  xlab("Position") +
  labs(fill = "Break point\ntype") +
  #scale_y_continuous(position = "right") + 
  coord_cartesian(c(inicio,final)) +
  theme(plot.margin = unit(c(0, 0.7, 0.7, 0.7), "cm")) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "white", colour = "white"),
        legend.title = element_text(size = 5)) #+
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#p44




num_iso = max(c(length(unique(df_pos[df_pos$perc >= abundance,"y_pos"])), length(unique(df_pos[df_pos$perc >= abundance,"star_stop"])) + 1), 4)

p_vertical_2 = cowplot::plot_grid(p11, p44, ncol=1, align='v', rel_heights = c(num_iso/3,2), axis = "lr")
#p_vertical_2

ggsave(path = plot_path, filename = paste(run_name, ".combined_plot_vertical.pdf", sep = ""), p_vertical_2, device="pdf", width = 10,
       height = 2 + num_iso/3, limitsize = F)
ggsave(path = plot_path, filename = paste(run_name, ".combined_plot_vertical.jpeg", sep = ""), p_vertical_2, device="jpeg", width = 10,
       height = 2 + num_iso/3, limitsize = F)













##########
# Tables #
##########
# num_reads_initial
# num_reads_post_trimming
# num_reads_post_exon_filtering
# 
# length(setdiff(gff3$id, gff3_filtered$id))
# length(setdiff(gff3_filtered$id, gff3_filtered_by_exon$id))



# All reads - class
log_map = read.delim(log_map_path, header = FALSE, stringsAsFactors = FALSE)
no_map_df = log_map[grepl("No paths found for", log_map$V1),]
no_map_df = data.frame(read_ids = gsub("No paths found for ", "", no_map_df), group = "Unmapped", stringsAsFactors = FALSE)

if (identical(setdiff(gff3$id, gff3_filtered$id), character(0))){
  vector_reads = data.frame()
}else{
  vector_reads = data.frame(read_ids = setdiff(gff3$id, gff3_filtered$id), group = "Only_vector_reads", stringsAsFactors = FALSE)
}

no_consensous_breakpoint_reads = data.frame(read_ids = unique(exon_id_filter), group = "No_consensous_breakpoint_reads", stringsAsFactors = FALSE)
read_isoform_df = data.frame(read_ids = names(read_isoform), group = as.character(read_isoform), stringsAsFactors = FALSE)

if (!is.null(filter_sites_file)){
  final_id_filter_out_reads = data.frame(read_ids = unique(final_id_filter_out), group = "Filter_out_reads", stringsAsFactors = FALSE)
  all_read_group = rbind(no_map_df, vector_reads, no_consensous_breakpoint_reads, final_id_filter_out_reads, read_isoform_df)
}else{
  all_read_group = rbind(no_map_df, vector_reads, no_consensous_breakpoint_reads, read_isoform_df)
}


write.table(all_read_group, file = paste(plot_path, "/", run_name, ".read_clasification.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




# Isoform freq
if (!is.null(filter_sites_file)){
  other_groups_df = data.frame(read_isoform = c("Unmapped", "Only_vector_reads", "No_consensous_breakpoint_reads", "Filter_out_reads"), 
                               Freq = c(nrow(no_map_df), nrow(vector_reads), nrow(no_consensous_breakpoint_reads), nrow(final_id_filter_out_reads)),
                               perc = c("-", "-", "-","-"), stringsAsFactors = FALSE)
}else{
  other_groups_df = data.frame(read_isoform = c("Unmapped", "Only_vector_reads", "No_consensous_breakpoint_reads"), 
                               Freq = c(nrow(no_map_df), nrow(vector_reads), nrow(no_consensous_breakpoint_reads)),
                               perc = c("-", "-", "-"), stringsAsFactors = FALSE)
}
all_groups = rbind(other_groups_df, isoform_frequencies)
all_groups$perc_total = round(all_groups$Freq*100/sum(all_groups$Freq), 2)
colnames(all_groups) =  c("Group", "N_of_reads", "Prerc_partial", "Perc_total")

write.table(all_groups, file = paste(plot_path, "/", run_name, ".isoform_freq.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




# Break points
start_position$type = "start"
end_position$type = "end"

start_position = merge(start_position, data.frame(table(gff3_filtered$start_tag)), by.x = "breakp", by.y = "Var1")
end_position = merge(end_position, data.frame(table(gff3_filtered$end_tag)), by.x = "breakp", by.y = "Var1")

start_end_position = rbind(start_position, end_position)
start_end_position$relative_freq = round(start_end_position$freq/num_reads_post_trimming*100, 2)
start_end_position$increase = round((start_end_position$Freq - start_end_position$freq)/start_end_position$freq*100, 2)
start_end_position = start_end_position[,c("breakp", "type", "freq", "relative_freq", "left", "right", "Freq", "increase")]
colnames(start_end_position) = c("breakpoint", "type", "number_of_reads", "percent_of_reads", "left_point", "right_point", "final_number_of_reads", "increase_percentage")

write.table(start_end_position, file = paste(plot_path, "/", run_name, ".breakpoints_info.tsv", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



