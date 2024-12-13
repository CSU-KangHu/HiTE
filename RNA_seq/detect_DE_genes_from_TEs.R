# 检查并安装缺失的包
load_or_install <- function(package, from_bioconductor=FALSE) {
  if (!require(package, character.only = TRUE)) {
    if (from_bioconductor) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package)
    } else {
      install.packages(package, dependencies = TRUE)
    }
    library(package, character.only = TRUE)
  }
}

# 使用函数加载所需的包
load_or_install("argparser")
load_or_install("readr")
load_or_install("dplyr")
load_or_install("stringr")
load_or_install("tidyr")


# Create a parser for command-line arguments
p <- arg_parser("Gene expression and TE association analysis")

# Add command-line arguments for the gene expression file and gene-TE associations file
p <- add_argument(p, "gene_express_file", help = "Path to gene expression file (e.g., gene_express.table)", type = "character")
p <- add_argument(p, "gene_te_associations_file", help = "Path to gene-TE associations file (e.g., gene_te_associations.tsv)", type = "character")

# Parse the command-line arguments
args <- parse_args(p)

# Print the input files (optional)
cat("Gene expression file:", args$gene_express_file, "\n")
cat("Gene-TE associations file:", args$gene_te_associations_file, "\n")

# Read the gene expression file
gene_express <- read_delim(args$gene_express_file, 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) %>%
  mutate(across(-1, ~ sapply(strsplit(as.character(.), ","), function(x) {
    last_val <- tail(x, 1)
    # If the last value is "NA", replace it with actual NA
    if (last_val == "NA") NA else as.numeric(last_val)
  })))

# Read the gene-TE associations file
gene_te_associations <- read_delim(args$gene_te_associations_file, 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)

# Process gene-TE associations to add species and calculate distance
gene_te_associations <- gene_te_associations %>%
  mutate(
    Species = Genome_name,
    Distance = if_else(
      TE_end >= Gene_start & TE_start <= Gene_end,  # If TE is inside the gene
      0,  # Set distance to 0
      pmin(abs(TE_end - Gene_start), abs(TE_start - Gene_end))  # Otherwise, calculate the minimum distance
    )
  )

# Reshape gene expression data to long format
gene_express <- gene_express %>%
  pivot_longer(
    cols = -gene_id,               # Keep gene_id column, and reshape the others (species columns)
    names_to = "Species",           # Create a new column for species names
    values_to = "expression"        # Create a new column for expression values
  )

# Merge gene expression data with gene-TE associations
gene_position_express <- gene_te_associations %>%
  right_join(gene_express, by = c("Gene_name" = "gene_id", "Species"))

# Calculate species present and species with NA for each gene
species_info <- gene_position_express %>%
  group_by(Gene_name) %>%
  summarise(
    species_na = paste(unique(Species[is.na(Position)]), collapse = ", "),
    species_present = paste(unique(Species[!is.na(Position)]), collapse = ", ")
  )

# Group gene expression by gene and position
group_express <- gene_position_express %>%
  group_by(Gene_name, Position)

# Calculate the minimum, maximum, and absolute difference for each gene by position
diff_express_detail <- group_express %>%
  summarise(
    min_expression = min(expression, na.rm = TRUE),
    max_expression = max(expression, na.rm = TRUE),
    abs_diff_expression = abs(max(expression, na.rm = TRUE) - min(expression, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = Position,
              values_from = c(min_expression, max_expression, abs_diff_expression),
              values_fill = list(min_expression = NA, max_expression = NA, abs_diff_expression = NA)) %>%
  rename_with(~ gsub("min_expression", "min", .), starts_with("min_expression")) %>%
  rename_with(~ gsub("max_expression", "max", .), starts_with("max_expression")) %>%
  rename_with(~ gsub("abs_diff_expression", "abs_diff", .), starts_with("abs_diff_expression"))


# 这里我们认为两组表达量[A, ..., B]和[C, ..., D] （假设B>A, D>C, C>B, 数组有序）之间有显著差异，需要满足C>B, (C-B)>(B-A) and (C-B)>(D-C)
significant_direct <- diff_express_detail %>%
  mutate(
    Upstream_direct = if_else(
      max_Upstream < min_NA & (min_NA - max_Upstream > abs_diff_NA) & (min_NA - max_Upstream > abs_diff_Upstream),
      'down',
      if_else(
        min_Upstream > max_NA & (min_Upstream - max_NA > abs_diff_NA) & (min_Upstream - max_NA > abs_diff_Upstream),
        'up',
        'ns'
      )
    ),
    Inside_direct = if_else(
      max_Inside < min_NA & (min_NA - max_Inside > abs_diff_NA) & (min_NA - max_Inside > abs_diff_Inside),
      'down',
      if_else(
        min_Inside > max_NA & (min_Inside - max_NA > abs_diff_NA) & (min_Inside - max_NA > abs_diff_Inside),
        'up',
        'ns'
      )
    ),
    Downstream_direct = if_else(
      max_Downstream < min_NA & (min_NA - max_Downstream > abs_diff_NA) & (min_NA - max_Downstream > abs_diff_Downstream),
      'down',
      if_else(
        min_Downstream > max_NA & (min_Downstream - max_NA > abs_diff_NA) & (min_Downstream - max_NA > abs_diff_Downstream),
        'up',
        'ns'
      )
    )
  ) %>%
  select(Gene_name, Upstream_direct, Inside_direct, Downstream_direct)

# Calculate fold changes for different positions compared to NoInsertion
fold_changes <- group_express %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = Position, values_from = mean_expression, values_fill = list(mean_expression = NA)) %>%
  rename(NoInsertion = 'NA') %>%
  mutate(
    # Calculate fold change between Upstream, Inside, Downstream and NoInsertion
    Upstream_vs_NoInsertion_logFoldChange = log2(Upstream + 1) - log2(NoInsertion + 1),
    Inside_vs_NoInsertion_logFoldChange = log2(Inside + 1) - log2(NoInsertion + 1),
    Downstream_vs_NoInsertion_logFoldChange = log2(Downstream + 1) - log2(NoInsertion + 1),
  )

# Identify significant genes based on fold change and position relationships
significant_genes <- fold_changes %>%
  left_join(significant_direct, by='Gene_name') %>%
  mutate(
    Upstream_vs_NoInsertion_logFoldChange = if_else(Upstream_direct == 'ns', NA_real_, Upstream_vs_NoInsertion_logFoldChange),
    Inside_vs_NoInsertion_logFoldChange = if_else(Inside_direct == 'ns', NA_real_, Inside_vs_NoInsertion_logFoldChange),
    Downstream_vs_NoInsertion_logFoldChange = if_else(Downstream_direct == 'ns', NA_real_, Downstream_vs_NoInsertion_logFoldChange)
  ) %>%
  filter(
    (Upstream_vs_NoInsertion_logFoldChange > 1 & Upstream_direct == 'up') | 
      (Upstream_vs_NoInsertion_logFoldChange < -1 & Upstream_direct == 'down') | 
      (Inside_vs_NoInsertion_logFoldChange > 1 & Inside_direct == 'up') | 
      (Inside_vs_NoInsertion_logFoldChange < -1 & Upstream_direct == 'down') | 
      (Downstream_vs_NoInsertion_logFoldChange > 1 & Downstream_direct == 'up') | 
      (Downstream_vs_NoInsertion_logFoldChange < -1 & Downstream_direct == 'down')
  ) %>%
  mutate(
    max_abs_value = pmax(abs(Upstream_vs_NoInsertion_logFoldChange), 
                         abs(Inside_vs_NoInsertion_logFoldChange), 
                         abs(Downstream_vs_NoInsertion_logFoldChange), na.rm = TRUE)
  ) %>%
  arrange(desc(max_abs_value)) %>%
  select(-max_abs_value) %>%
  rename(
    'NoInsertion (mean_TPM)' = NoInsertion,
    'Upstream (mean_TPM)' = Upstream,
    'Downstream (mean_TPM)' = Downstream,
    'Inside (mean_TPM)' = Inside
  ) %>%
  select(-c(Upstream_direct, Inside_direct, Downstream_direct))

write.table(gene_position_express, "all_gene_TEs_details.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(significant_genes, "DE_genes_from_TEs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
