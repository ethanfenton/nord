#This should analyze qPCR data and give a Fold Change via the DDCt method
# Specify the following 4 variables:
  # file_path, housekeeping_gene, experimental_groups, control_groups
library(tidyverse)

#file_path <- "G:\\My Drive\\Labs\\Nord\\CHD8\\qPCR\\14-Jul-2022_CHD8_Cortex_Validation_plate2_fixed.csv"
#file_path <- "G:\\My Drive\\Labs\\Nord\\CHD8\\qPCR\\14-Jul-2022_CHD8_Cortex_Validation_plate1.csv"
#file_path <- "G:\\My Drive\\Labs\\Nord\\CHD8\\qPCR\\15-Jul-2022_CHD8_Cerebellum_Validation_plate1__new_Ct.csv"
#file_path <- "G:\\My Drive\\Labs\\Nord\\CHD8\\qPCR\\15-Jul-2022_CHD8_Cerebellum_Validation_plate1_new_Ct.csv"

file_path <- "G:\\My Drive\\Labs\\Nord\\CHD8\\qPCR\\18-Jul-2022_CHD8_Cortex_downreg_plate1.csv"
housekeeping_gene <- "Gapdh"

# the order of the control and experimental groups should align with one another
# since the ddct comparison will use the order of these lists
experimental_groups <- c('1036_cx', '1039_cx')
control_groups <- c('1037_cx', '1038_cx')

df <- read.csv2(file_path, row.names=NULL)
skips <- which(grepl("Well",df[,1]))
df <- read.csv2(file_path, skip=skips, sep=",")

#check for No Ct values
if (nrow(df[which(df$Ct == "No Ct"),])){
  print(paste("Warning there is a No Ct value in your dataset, consider removing this"))
  print(df[which(df$Ct == "No Ct"),])
  answer <- readline("Would you like to remove the No Ct well(s)? y/n \nDefaults to y")
  if (answer == "y" | answer == ""){
    df <- df[which(df$Ct != "No Ct"),]
  }
}

mean_cts <- unique(select(df, Sample.name, Mean.Ct, Gene, Std.Dev..Ct))
mean_cts <- mean_cts[!apply(mean_cts == "", 1, all),]

for (row in 1:nrow(mean_cts)){
  if (mean_cts$Std.Dev..Ct[row] > 0.5) {
    print(paste("Warning", mean_cts$Sample.name[row], mean_cts$Gene[row],"has a high StDev"))
    print(df[df$Sample.name == mean_cts$Sample.name[row] & df$Gene==mean_cts$Gene[row],])
  }
}

deltaCt <- c()
for (row in rownames(mean_cts)){
  for (row2 in rownames(mean_cts)){
    if(mean_cts[row,"Sample.name"] == mean_cts[row2,'Sample.name'] & mean_cts[row2,"Gene"]==housekeeping_gene){
      deltaCt[row] <- as.numeric(mean_cts[row,"Mean.Ct"]) - as.numeric(mean_cts[row2,"Mean.Ct"])
    }
  }
}

mean_cts$deltaCt <- deltaCt
ddCt_list <- c()
group_list <- c()
gene_list <- c()
for(gene in unique(mean_cts$Gene)){
  for(group in 1:length(experimental_groups)){
    dCt_value_1 <- mean_cts[mean_cts$Gene==gene & mean_cts$Sample.name==experimental_groups[group], "deltaCt"]
    dCt_value_2 <- mean_cts[mean_cts$Gene==gene &  mean_cts$Sample.name==control_groups[group], "deltaCt"]
    ddCt_value <- dCt_value_1 - dCt_value_2
    if (length(ddCt_value) == 0){
      ddCt_value <- 0
    }
    ddCt_list <- append(ddCt_list, ddCt_value)
    group_list <- append(group_list, group)
    gene_list <- append(gene_list, gene)
  }
}
fold_change_list <- 2^-(ddCt_list) 
ddCt <- data.frame(group=group_list, gene=gene_list, ddCt=ddCt_list, fold_change=fold_change_list)
ddCt$group[which(ddCt$group==1)] <- "male"
ddCt$group[which(ddCt$group==2)] <- "female"
print(ddCt)
                   