#AUC
library(DescTools)
files <- list.files(pattern = "\\.txt$")
names <- sub("\\.txt","",files)
gsub("^([^_]*)_([^_]*)_", "", names) -> groupname
data_list <- list()
cal_list <- list()
for (i in seq_along(files)) {
  file <- files[[i]]
  name <- names[i]
  data <- read.table(file, header = TRUE,sep = "\t",skip=3)
  data_list[[file]] <- data
  AUC(data_list[[file]]$X,data_list[[file]]$Y,method = "step") -> cal_list[[name]]
}
as.data.frame(cal_list) -> cal_list
write.table(t(cal_list),paste(groupname,"_AUC.txt",sep = ""),sep = "\t",row.names = F,quote = F)
