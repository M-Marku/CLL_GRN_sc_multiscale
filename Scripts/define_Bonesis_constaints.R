input_files <- c("genes_M0_cond.txt", "genes_M3_cond.txt", "genes_M27_cond.txt")

# Genes that appear anywhere in the filtered GRN
all_GRN_genes <- unique(c(CLL_GRN$Source, CLL_GRN$Target))

# Loop over the files
for (file in input_files) {
  
  # Read gene list
  gene_list <- read.table(file, stringsAsFactors = FALSE)[,1]
  
  # Keep only genes present in filtered GRN
  matched_genes <- gene_list[gene_list %in% all_GRN_genes]
  
  # Output filename (e.g., genes_list1_filtered.txt)
  output_file <- sub(".txt$", "_filtered.txt", file)
  
  # Write filtered gene list
  write.table(
    matched_genes,
    file = output_file,
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  message("Processed: ", file, " → ", output_file)
}
