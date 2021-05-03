# This is a script to make the Yeo 17 parcellation
library(ciftiTools)
ciftiTools.setOption('wb_path','/Applications/workbench')
yeo <- read_cifti("Yeo2011_17Networks_ContiguousParcels.32k_fs_LR.dlabel.nii")
plot(yeo, fname = "~/Desktop/506_yeo_parcellation.png", borders = FALSE,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")
yeo$meta$cifti
# Need to grab the labels and rework them into the 17-network parcellation
new_labels <- rownames(yeo$meta$cifti$labels$`#1`)
new_labels <- new_labels[-1] # Remove the "???" label
new_labels <- substring(new_labels,15,18) # Remove the "[rl]h.17Networks_" string from the labels
new_labels <- gsub("-.*","",new_labels) # Gets rid of the subnetwork label
new_labels <- sub("Defi","Medial Wall", new_labels)
# Make a color palette for the 17 classes (+ medial wall)
better_labels <- data.frame(
  mine = 1:17,
  theirs = c(
    "Visual A", "Visual B", "Somatomotor A", "Somatomotor B", "Dorsal Attention B",
    "Dorsal Attention A", "Salience/VenAttn A", "Salience/VenAttn B", "Limbic A", "Limbic B",
    "Control C", "Control A", "Control B", "Temporal Parietal", "Default C",
    "Default A", "Default B"
  ),
  colors = c("purple","red","cornflowerblue","cyan","forestgreen",
             "green","mediumpurple2","pink","lightgreen","darkolivegreen4",
             "grey50","orange","darkred","blue","darkblue",
             "yellow","indianred2")
)
rgb_colors_pal <- col2rgb(better_labels$colors) / 255
old_labels <- yeo$meta$cifti$labels$`#1`
label_df <- data.frame(
  Key = c(1:19),
  Red = c(1,rgb_colors_pal[1,],1),
  Green = c(1,rgb_colors_pal[2,],1),
  Blue = c(1,rgb_colors_pal[3,],1),
  Alpha = c(0,rep(1,17),0)
)
rownames(label_df) <- c("???",better_labels$theirs, "Medial.Wall")
num_key <- as.numeric(new_labels)
num_key[is.na(num_key)] <- 18
num_key <- c(1, 1 + num_key) # Added 1 for the first key value back in for the "???" label
yeo_17 <- yeo
yeo_17$data$cortex_left <- as.matrix(num_key[1 + yeo_17$data$cortex_left[,1]])
yeo_17$data$cortex_right <- as.matrix(num_key[1 + yeo_17$data$cortex_right[,1]])
yeo_17$meta$cifti$labels$`#1` <- label_df

plot(yeo_17, fname = "~/Desktop/506_yeo_17_parcellation.png", borders = F,
     surfL = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.L.inflated.32k_fs_LR.surf.gii",
     surfR = "/Volumes/GoogleDrive/My Drive/MEJIA_LAB/data/Q1-Q6_R440.R.inflated.32k_fs_LR.surf.gii")

write_cifti(yeo_17, "~/Desktop/yeo_17_parcellation.dlabel.nii")
