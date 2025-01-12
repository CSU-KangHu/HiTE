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
load_or_install("ggplot2")

# Create a parser for command-line arguments
p <- arg_parser("Gene expression and TE association analysis")

# Add command-line arguments for the gene expression file and gene-TE associations file
p <- add_argument(p, "gene_express_file", help = "Path to gene expression file (e.g., gene_express.table)", type = "character")
p <- add_argument(p, "gene_te_associations_file", help = "Path to gene-TE associations file (e.g., gene_te_associations.tsv)", type = "character")
# 添加 min_sample_threshold 和 min_diff_threshold 参数
p <- add_argument(p, "--min_sample_threshold", help = "Minimum sample threshold (default: 5)", type = "numeric", default = 5)

p <- add_argument(p, "--min_diff_threshold", help = "Minimum diff threshold (default: 10)", type = "numeric", default = 10)

# Parse the command-line arguments
args <- parse_args(p)

# Print the input files (optional)
cat("Gene expression file:", args$gene_express_file, "\n")
cat("Gene-TE associations file:", args$gene_te_associations_file, "\n")
gene_express_file <- args$gene_express_file
gene_te_associations_file <- args$gene_te_associations_file
min_sample_threshold <- args$min_sample_threshold
min_diff_threshold <- args$min_diff_threshold

# Read the gene expression file
gene_express <- read_delim(gene_express_file,
                           delim = "\t", escape_double = FALSE,
                           trim_ws = TRUE) %>%
  mutate(across(-1, ~ sapply(strsplit(as.character(.), ","), function(x) {
    last_val <- tail(x, 1)
    # If the last value is "NA", replace it with actual NA
    if (last_val == "NA") NA else as.numeric(last_val)
  })))

# Read the gene-TE associations file
gene_te_associations <- read_delim(gene_te_associations_file,
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
  ) %>%
  filter(!is.na(expression))

# Merge gene expression data with gene-TE associations
gene_position_express <- gene_te_associations %>%
  right_join(gene_express, by = c("Gene_name" = "gene_id", "Species")) %>%
  mutate(Position = replace_na(Position, 'No_Insertion'))

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
    abs_diff_expression = abs(max(expression, na.rm = TRUE) - min(expression, na.rm = TRUE)),
    Group_Size = n(), .groups = 'drop'
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = Position, 
              values_from = c(min_expression, max_expression, abs_diff_expression, Group_Size), 
              values_fill = list(min_expression = NA, max_expression = NA, abs_diff_expression = NA, Group_Size = NA)) %>%
  rename_with(~ gsub("min_expression", "min", .), starts_with("min_expression")) %>%
  rename_with(~ gsub("max_expression", "max", .), starts_with("max_expression")) %>%
  rename_with(~ gsub("abs_diff_expression", "abs_diff", .), starts_with("abs_diff_expression"))

# 这里我们认为两组表达量[A, ..., B]和[C, ..., D] （假设B>A, D>C, C>B, 数组有序）之间有显著差异，需要满足abs(C-B) > min_diff_threshold
# 定义计算方向的函数
assign_direction <- function(min_no_insert, max_insert, min_insert, max_no_insert, min_diff_threshold, size_no_insert, size_insert, min_sample_threshold = 5) {
  # 检查 size_no_insert 和 size_insert 是否为 NA
  is_na <- is.na(size_no_insert) | is.na(size_insert)
  
  # 检查样本数是否满足阈值
  below_threshold <- max(size_no_insert, size_insert) < min_sample_threshold
  
  # 判断表达变化方向
  result <- case_when(
    is_na | below_threshold ~ 'ns',  # 如果样本数为 NA 或样本数不足，直接返回 'ns'
    (min_no_insert - max_insert) > min_diff_threshold ~ 'down',
    (min_insert - max_no_insert) > min_diff_threshold ~ 'up',
    TRUE ~ 'ns'
  )
  
  return(result)
}

# 定义计算差异绝对值的函数
calculate_diff <- function(min_no_insert, max_insert, min_insert, max_no_insert, size_no_insert, size_insert, min_sample_threshold = 5) {
  # 检查 size_no_insert 和 size_insert 是否为 NA
  is_na <- is.na(size_no_insert) | is.na(size_insert)
  
  # 检查样本数是否满足阈值
  below_threshold <- max(size_no_insert, size_insert) < min_sample_threshold
  
  # 计算表达差异
  case_when(
    is_na ~ 0,  # 如果样本数为 NA，直接返回 0
    below_threshold ~ 0,  # 如果样本数不足，直接返回 0
    (min_no_insert - max_insert) > 0 ~ abs(min_no_insert - max_insert),
    (min_insert - max_no_insert) > 0 ~ abs(min_insert - max_no_insert),
    TRUE ~ 0  # 如果没有显著差异，返回 0
  )
}
# 生成方向列
significant_direct <- diff_express_detail %>%
  mutate(
    Upstream_direct = assign_direction(
      min_No_Insertion, max_Upstream, min_Upstream, max_No_Insertion, min_diff_threshold,
      Group_Size_No_Insertion, Group_Size_Upstream, min_sample_threshold
    ),
    Inside_direct = assign_direction(
      min_No_Insertion, max_Inside, min_Inside, max_No_Insertion, min_diff_threshold,
      Group_Size_No_Insertion, Group_Size_Inside, min_sample_threshold
    ),
    Downstream_direct = assign_direction(
      min_No_Insertion, max_Downstream, min_Downstream, max_No_Insertion, min_diff_threshold,
      Group_Size_No_Insertion, Group_Size_Downstream, min_sample_threshold
    )
  ) %>%
  select(c(Gene_name, Upstream_direct, Inside_direct, Downstream_direct))

# 生成差异绝对值列
significant_diff <- diff_express_detail %>%
  mutate(
    Upstream_diff = calculate_diff(min_No_Insertion, max_Upstream, min_Upstream, max_No_Insertion, Group_Size_No_Insertion, Group_Size_Upstream, min_sample_threshold),
    Inside_diff = calculate_diff(min_No_Insertion, max_Inside, min_Inside, max_No_Insertion, Group_Size_No_Insertion, Group_Size_Inside, min_sample_threshold),
    Downstream_diff = calculate_diff(min_No_Insertion, max_Downstream, min_Downstream, max_No_Insertion, Group_Size_No_Insertion, Group_Size_Downstream, min_sample_threshold)
  )  %>%
  select(c(Upstream_diff, Inside_diff, Downstream_diff))
# 合并方向列和差异绝对值列
significant_direct_diff <- significant_direct %>%
  bind_cols(significant_diff %>% select(Upstream_diff, Inside_diff, Downstream_diff))


# 宽表转长表
long_significant_direct <- significant_direct_diff %>%
  select(c(Gene_name, Upstream_direct, Inside_direct, Downstream_direct)) %>%
  pivot_longer(
    cols = starts_with("Upstream") | starts_with("Inside") | starts_with("Downstream"),
    names_to = c("Insert_type"),
    values_to = "direct",
    names_pattern = "(.*)_direct"
  )
long_significant_diff <- significant_direct_diff %>%
  select(c(Gene_name, Upstream_diff, Inside_diff, Downstream_diff)) %>%
  pivot_longer(
    cols = starts_with("Upstream") | starts_with("Inside") | starts_with("Downstream"),
    names_to = c("Insert_type"),
    values_to = "diff",
    names_pattern = "(.*)_diff"
  )
long_significant_direct_diff <- long_significant_direct %>%
  left_join(long_significant_diff, by = c("Gene_name", "Insert_type"))

# Calculate fold changes for different positions compared to No_Insertion
fold_changes <- group_express %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = Position, values_from = mean_expression, values_fill = list(mean_expression = NA)) %>%
  mutate(
    # Calculate fold change between Upstream, Inside, Downstream and No_Insertion
    Upstream_vs_NoInsertion_logFoldChange = log2(Upstream + 1) - log2(No_Insertion + 1),
    Inside_vs_NoInsertion_logFoldChange = log2(Inside + 1) - log2(No_Insertion + 1),
    Downstream_vs_NoInsertion_logFoldChange = log2(Downstream + 1) - log2(No_Insertion + 1),
  )

long_fold_changes <- fold_changes %>%
  select(-c(Downstream, Inside, Upstream, No_Insertion)) %>%
  pivot_longer(
    cols = starts_with("Upstream_") | starts_with("Inside_") | starts_with("Downstream_"),
    names_to = c("Insert_type"),
    values_to = "fold_change",
    names_pattern = "(.*)_vs_NoInsertion_logFoldChange"
  )

valid_genes <- long_fold_changes %>%
  left_join(long_significant_direct_diff, by = c("Gene_name", "Insert_type")) %>%
  mutate(
    significant = if_else(
      (fold_change > 1 & direct == 'up') | (fold_change < -1 & direct == 'down'),
      "Significant",
      "Not Significant"
    ),
    unique_gene_name = paste(Gene_name, Insert_type, sep = "_"),
    log10_diff = log10(diff+1)
  ) %>%
  mutate(
    significant = replace_na(significant, "Not Significant") # 将 NA 替换为 "Not Significant"
  ) %>%
  filter(log10_diff != 0 & log10_diff != Inf)


# 绘制火山图
# 筛选 fold_change 最大的显著点
top_fold_change <- valid_genes %>%
  filter(significant == "Significant") %>%
  arrange(desc(abs(fold_change))) %>%
  slice_head(n = 5)  # 选择前 5 个点

# 筛选 log10_diff 最大的显著点
top_log10_diff <- valid_genes %>%
  filter(significant == "Significant") %>%
  arrange(desc(log10_diff)) %>%
  slice_head(n = 5)  # 选择前 5 个点

# 合并筛选结果
top_genes <- bind_rows(top_fold_change, top_log10_diff) %>%
  distinct()  # 去重

pdf("DE_genes_from_TEs.pdf", width = 8, height = 6)
volcano_plot <- ggplot(valid_genes, aes(x = fold_change, y = log10_diff, color = direct)) +
  geom_point(alpha = 0.6, size = 2) +  # 绘制所有点
  geom_point(
    data = subset(valid_genes, significant == "Significant"),  # 仅对显著基因
    aes(shape = significant), size = 3, stroke = 1, color = "black"  # 添加黑色边框
  ) +
  geom_text(
    data = top_genes,  # 为部分显著基因添加标签
    aes(label = Gene_name),  # 使用 Gene_name 作为标签
    vjust = -1.5, size = 3, color = "black", check_overlap = TRUE  # 调整标签位置和样式
  ) +
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "grey")) +  # 设置颜色
  scale_shape_manual(values = c("up" = 1, "down" = 1), guide = "none") +  # 设置形状并隐藏图例
  theme_minimal() +  # 使用简洁主题
  labs(
    title = "",
    x = "log2(Fold Change)",
    y = "log10(Diff+1)",
    color = "Direction",
    shape = ""
  ) +
  geom_hline(yintercept = log10(min_diff_threshold), linetype = "dotdash", color = "black") +  # 添加显著性阈值线
  geom_vline(xintercept = c(-1, 1), linetype = "dotdash", color = "black")  # 添加 Fold Change 阈值线

# 显示图形
print(volcano_plot)
dev.off()

significant_genes <- valid_genes %>%
  filter(significant == 'Significant') %>%
  arrange(desc(abs(fold_change)))

write.table(significant_genes, "DE_genes_from_TEs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(gene_position_express, "all_gene_TEs_details.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)