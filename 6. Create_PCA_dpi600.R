# Edit this script to run in R
# Change the working directory, the prefix which should be the same as your plink output file, the metadata file.
# Change the labels and colours in the plot
# Change the name of the plot

library(showtext)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()
library(viridis)
library(scales)

workdir <- "/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/pca" # Working directory with plink files
prefix <- "mito_only_melas_plusglobal" # Prefix for plink files
metadata <- "/mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/genomics_database_melas2019plusglobal/genomics_database_melas2019plusglobal_vcf/melas_2019_plusglobal_filtering/pca/metadata_melasplusglobal_clusters_v2.csv" # File path to metadata

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "sample"))

dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

# PCA #
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE) # Multidimensional Scaling - might take a while
# saveRDS(cmd, paste0(prefix, ".dist.rds") # save to RDS format
#cmd <- readRDS(file.path(workdir, paste0(prefix, ".dist.rds"))
vars <- calc_variance_explained(cmd) # Calculations of variance explained

# Overlay region, country info
# Convert cmdscale output to a dataframe and set row names as a 'sample' column
df <- as.data.frame(cmd$points, stringsAsFactors = FALSE)
df$sample <- rownames(df)  # Set rownames as the 'sample' column

# Add metadata columns
df$country <- gsub("_", " ", desc$country)
df$island <- gsub("_", " ", desc$island)
df$species <- gsub("_", " ", desc$species)
df$Population <- gsub("_", " ", desc$Population)

# Correctly reorder the dataframe to make 'sample' the first column
df <- df[, c("sample", setdiff(names(df), "sample"))]

# Rename PCA columns to start with 'PC', if they start with 'V' 
# (assuming all PCA columns currently start with 'V', adjust if different)
colnames(df)[-c(1, ncol(df)-1, ncol(df))] <- gsub("^V", "PC", colnames(df)[-c(1, ncol(df)-1, ncol(df))])

color_by <- "Population" # specify if colored by region or country

## changing colour scheme to be with viridis, with color_by being a discrete variable

my_colours <- c("An.melas Cameroon" = "#7678ed", "An.melas The Gambia" = "darkorange2", 
                "An.melas Bijagos Group A" = "deeppink", "An.melas Bijagos Group B" = "#f0e442")

# plot
png("mito_only_PCA_coloured_by_cluster_dpi_600_larger_text.png", width = 7000, height = 7000, res = 600) 
ggplot(data = df, aes(x = PC1, y = PC2, color = !!sym(color_by))) +
    geom_point(size = 5) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"), y = paste0("PC2", " (", vars["PC2"], "%)"), title = "Mitochondria") +
    scale_color_manual(values = my_colours) +
    scale_x_continuous(labels = label_number()) +
    scale_y_continuous(labels = label_number()) +
    theme_classic(base_size = 120) +
    theme(legend.position = "bottomleft", 
    legend.direction = "vertical", 
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 10, r = 40, b = 30, l = 10, unit = "pt"),
    legend.text = element_text(size = 120),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5))
dev.off()
#
