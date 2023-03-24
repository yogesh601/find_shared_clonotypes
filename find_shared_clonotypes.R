## Function for looking for clonal/ clonotype expension in seurat object
# With this function you can directly visualize the clonotypes from your classical seurat object

```{r}
library(Seurat)
library(VennDiagram)

find_shared_clonotypes <- function(seurat_obj, samples, min_cells = 1) {
  
  # Extract the clonotype information from the Seurat object
  clonotypes <- FetchData(seurat_obj, vars = c("cdr3s_aa", "Tissue")) # You can change the Tissue/condition/celltype anythingamong which you want to compare the clonal expension also change the same at every place where it appear
  
  # Filter out clonotypes that appear in fewer than min_cells samples
  clonotypes_counts <- table(clonotypes$cdr3s_aa)
  clonotypes_filtered <- names(clonotypes_counts)[clonotypes_counts >= min_cells]
  clonotypes <- clonotypes[clonotypes$cdr3s_aa %in% clonotypes_filtered, ]
  
  # Filter out clonotypes that are not present in the specified samples
  clonotypes <- clonotypes[clonotypes$Tissue %in% samples, ] ## Change name
  
 # Find shared clonotypes for each sample
  shared_clonotypes <- list()
  total_clonotypes <- list()
  for (sample in samples) {
    clonotypes_sample <- clonotypes[clonotypes$Tissue == sample, "cdr3s_aa"] ## Change name
    shared_clonotypes[[sample]] <- unique(clonotypes_sample[duplicated(clonotypes_sample)])
    total_clonotypes[[sample]] <- unique(clonotypes_sample)
  }
  
  # Print summary of shared clonotypes and total clonotypes
  cat("Summary of shared clonotypes:/n")
  for (sample in samples) {
    cat(paste0(sample, ": ", length(shared_clonotypes[[sample]]), "/n"))
  }
  cat("/n")
  cat("Summary of total clonotypes:/n")
  for (sample in samples) {
    cat(paste0(sample, ": ", length(total_clonotypes[[sample]]), "/n"))
  }
  
  # Plot venn diagram of shared clonotypes
  venn_list <- list()
  for (sample in samples) {
    venn_list[[sample]] <- shared_clonotypes[[sample]]
  }
  if (length(venn_list) > 1) {
    venn.plot <- venn.diagram(venn_list, filename=NULL, fill=rainbow(length(venn_list)), alpha=c(0.5), label.col=c("black"), cat.cex = 3, fontface=3, cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "orange", "green"), cat.fontface=1)
    grid.draw(venn.plot)
  } else {
    cat("Cannot plot venn diagram with only one sample./n")
  }
  
  # Return list of shared clonotypes
  return(list(shared_clonotypes = shared_clonotypes, total_clonotypes = total_clonotypes))

}

# Here now just run the function
cells <- c("Inflamed", "Tumor", "Blood", "Normal", "Lymphnode") # assign the Tissue/celltype/condition
shared_clonotypes <-find_shared_clonotypes(Seurat_object, cells)

## Thats all now you can visualize the output of this shared_clonotypes in venn diagram or other form. Although this function will also plot the venn diagram for you
# Note dont exceed the limitatio of venn diagram I think It will take maximum 6 type of tissues.
# Write me for any error.
# Thanks

