#loading necessary packages
library(ggplot2)
library(ggpubr)
library(eoffice)
#creat two lists to store data and plot
data_list <- list()
plot_list <- list()
#get work directory as a character string
dir <- getwd()
#list all files in the directory
files <- list.files(dir, pattern = "^D3.*xpo\\.txt$")
#remove the suffix of the file names
names <- unique(sub("_r[123]","",files))
names <- sub("\\.txt","",names)
#read data from each file and store them in data_list
for (i in seq_along(files)) {
  file <- files[i]
  name <- names[i]
  data <- read.table(file, header = TRUE,sep = "\t",skip=3,row.names = 1)
  data_list[[file]] <- data
}
#merge data from different data frame
Reduce(function(x,y) merge(x,y,by="nm",all=T),lapply(data_list, function(x) { x$nm <- rownames(x); x })) -> merged_data
#rename the columns and rows
colnames(merged_data) <- c("nm", names(data_list))
rownames(merged_data) <- merged_data$nm
#remove the nm column
merged_data[,-1] -> merged_data
#calculate the mean of each row
rowMeans(merged_data[,1:ncol(merged_data)]) -> merged_data$mean
#add a column to store the rownames
merged_data$rowname <- rownames(merged_data)
#plot the data and save them as pdf files
for (i in names) {
  plot_list[[i]] <- ggplot(merged_data, aes(x = rowname, y = mean, group = 1)) +
   geom_point(size = 0.1) +
   geom_line()+ 
   theme_bw() + 
   xlab("Wavelength (nm)") + 
   ylab("Relative intensity") + 
   theme(legend.position = "none")+
   ggtitle(name)+
   scale_x_discrete(breaks =  seq(380, 710, by = 70))
  topptx(plot_list[[i]],filename = paste0(names,".pptx"), width = 6, height = 4)
  }
#set the file name
filename <- paste0(names,("_avg"))
#save the data as a txt file
write.table(merged_data,paste0(filename),sep="\t",row.names = F,quote = F)
#merge the multiple plots into one use ggpubr
combined_plot <- ggarrange(plotlist = plot_list, ncol = 4, nrow = 5)
#get the file which end by _avg
file_pairs <- combn(list.files(pattern = "\\_avg$"),2,simplify = F)
#creat two lists to store data and plot
pair_list <- list()
plot_list <- list()
#read data from each file and store them in data_list(pairwise)
for (i in seq_along(file_pairs)){
  pair <- file_pairs[[i]]
  pair_data <- list()
  #load data from each file from the "pair"
  pair_data[[pair[1]]] <- read.table(pair[1], header = TRUE, sep = "\t")
  #add a column to store the group name(filename)
  pair_data[[pair[1]]]$group <- pair[1]
  #remove the first three columns
  pair_data[[pair[[1]]]] <- pair_data[[pair[[1]]]][,-(1:3)]
  #same as above,but for the second file of the "pair"
  pair_data[[pair[2]]] <- read.table(pair[2], header = TRUE, sep = "\t")
  pair_data[[pair[2]]]$group <- pair[2]
  pair_data[[pair[[2]]]] <- pair_data[[pair[[2]]]][,-(1:3)]
  #merge the data from two files as the same header
  rbind(pair_data[[pair[[1]]]],pair_data[[pair[[2]]]]) -> pair_data
  #plot the data and save them into the plot_list
  plot_list[[paste(pair[1], pair[2], sep = "_vs_")]] <- ggplot(pair_data, aes(x = rowname, y = mean, group = group, color = group)) +
  geom_point(size = 0.1) +
  geom_line() + 
  theme_bw() + 
  xlab("Wavelength (nm)") + 
  ylab("Relative intensity") + 
  theme(legend.position = "right") + 
  ggtitle(paste(pair[1], pair[2], sep = " vs ")) + 
  scale_x_continuous(breaks = seq(380, 710, by = 70))
  #use the ggsave function to save the plot as pdf files
  ggsave(filename = paste(pair[1],"_vs_",pair[2],".pdf"),plot_list[[paste(pair[1], pair[2], sep = "_vs_")]], width = 6, height = 4)
  #export the pair wise data to the pair_data to check
  pair_list[[paste(pair[1], pair[2], sep = "_vs_")]] <- pair_data
}
#get the file which end by _avg
mean_files <- list.files(pattern = "\\_avg$")
#creat two lists to store data and plot
filename <- list()
plotname <- list()
#read data from each file and store them in data_list(pairwise)
for (i in seq_along(mean_files)){
  files <- mean_files[i]
  read.table(files,header = TRUE,sep = "\t") -> filename[[files]]
  p <- ggplot(filename[[files[1]]], aes(x = rowname, y = mean, group = 1)) + geom_point(size = 0.1) + geom_line() + theme_bw() + xlab("Wavelength (nm)") + ylab("Relative intensity") + theme(legend.position = "right") + ggtitle(files[1]) + scale_x_continuous(breaks = seq(380, 710, by = 70))
  ggsave(filename = paste(files,".pdf"),p, width = 6, height = 4)
}
