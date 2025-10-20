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

# Parse the command-line arguments
args <- parse_args(p)

# Print the input files (optional)
cat("Gene expression file:", args$gene_express_file, "\n")
cat("Gene-TE associations file:", args$gene_te_associations_file, "\n")
gene_express_file <- args$gene_express_file
gene_te_associations_file <- args$gene_te_associations_file

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
  mutate(
    Position = replace_na(Position, 'No_Insertion'),
    expression = as.numeric(expression)  # 转换为数值
    ) %>%
  distinct(Gene_name, Position, Species, .keep_all = TRUE)


# 定义一个函数，用于安全地进行 t 检验
safe_t_test <- function(x, y) {
  # 检查数据是否满足 t 检验的条件
  if (length(x) > 1 && length(y) > 1 && is.numeric(x) && is.numeric(y)) {  # 确保数据是数值型且每组至少有两个样本
    tryCatch({
      return(t.test(x, y)$p.value)  # 返回 p 值
    }, error = function(e) {
      return(NA)  # 如果 t 检验出错，返回 NA
    })
  } else {
    return(NA)  # 如果数据不满足条件，返回 NA
  }
}

p_value_results <- gene_position_express %>%
  group_by(Gene_name) %>%
  summarise(
    p_value_upstream = ifelse(
      length(expression[Position == "Upstream"]) > 0 && length(expression[Position == "No_Insertion"]) > 0,
      safe_t_test(expression[Position == "Upstream"], expression[Position == "No_Insertion"]),
      NA
    ),
    p_value_downstream = ifelse(
      length(expression[Position == "Downstream"]) > 0 && length(expression[Position == "No_Insertion"]) > 0,
      safe_t_test(expression[Position == "Downstream"], expression[Position == "No_Insertion"]),
      NA
    ),
    p_value_inside = ifelse(
      length(expression[Position == "Inside"]) > 0 && length(expression[Position == "No_Insertion"]) > 0,
      safe_t_test(expression[Position == "Inside"], expression[Position == "No_Insertion"]),
      NA
    )
  )

p_adjust_value_results <- p_value_results %>%
  mutate(
    fdr_Upstream = p.adjust(p_value_upstream, method = "fdr"),
    fdr_Downstream = p.adjust(p_value_downstream, method = "fdr"),
    fdr_Inside = p.adjust(p_value_inside, method = "fdr")
  )

long_p_adjust_value <- p_adjust_value_results %>%
  select(-starts_with("p_value_")) %>%  # 排除某些列
  pivot_longer(
    cols = starts_with("fdr_"),  # 选择以 "exp_" 开头的列
    names_to = "Insert_type",      # 新列名，存储原始列名
    values_to = "P_adjust_value",    # 新列名，存储原始列的值
    names_pattern = "fdr_(.*)"
  )


# Calculate fold changes for different positions compared to No_Insertion
fold_changes <- gene_position_express %>%
  group_by(Gene_name, Position) %>%
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
  inner_join(long_p_adjust_value, by = c("Gene_name", "Insert_type")) %>%
  mutate(
    significant = case_when(
      abs(fold_change) > 1 & P_adjust_value < 0.05 ~ "Significant",
      TRUE ~ "Not Significant"
    ),
    direct = case_when(
      significant == "Significant" & fold_change > 0 ~ "up",  # Significant且fold_change > 0
      significant == "Significant" & fold_change < 0 ~ "down",  # Significant且fold_change < 0
      significant == "Not Significant" ~ "ns",  # Not Significant时，direct为"ns"
      TRUE ~ "ns"  # 默认情况下，如果没有匹配条件，也设置为"ns"
    )
  ) %>%
  # 优先保留 significant 为 "Significant" 的记录
  arrange(Gene_name, desc(significant == "Significant")) %>%
  # 按照 Gene_name 去重，只保留第一条记录
  distinct(Gene_name, .keep_all = TRUE) %>%
  filter(!is.na(fold_change), !is.na(P_adjust_value), P_adjust_value > 0)


# 绘制火山图
# 筛选 fold_change 最大的显著点
top_fold_change <- valid_genes %>%
  filter(significant == "Significant") %>%
  arrange(desc(abs(fold_change))) %>%
  slice_head(n = 5)  # 选择前 5 个点

# 筛选 log10_diff 最大的显著点
top_p_adjust_value <- valid_genes %>%
  filter(significant == "Significant") %>%
  arrange(P_adjust_value) %>%
  slice_head(n = 5)  # 选择前 5 个点

# 合并筛选结果
top_genes <- bind_rows(top_fold_change, top_p_adjust_value) %>%
  distinct()  # 去重

pdf("DE_genes_from_TEs.pdf", width = 8, height = 6)
volcano_plot <- ggplot(valid_genes, aes(x = fold_change, y = -log10(P_adjust_value))) +
  geom_point(aes(color = direct), size = 2) +
  geom_text(
    data = top_genes,  # 为部分显著基因添加标签
    aes(label = Gene_name),  # 使用 Gene_name 作为标签
    vjust = -1.5, size = 5, color = "black", check_overlap = TRUE  # 调整标签位置和样式
  ) +
  geom_point(aes(color = direct), size = 2, alpha = ifelse(valid_genes$direct == "ns", 0.2, 1)) +  # "ns"点透明度降低
  scale_color_manual(values = c("up" = "red", "down" = "blue", "ns" = "gray")) +  # 上调为红色，下调为蓝色，其他为灰色
  theme_minimal() +
  labs(x = "log2(Fold Change)", y = "-log10(P-adjusted)", title = "") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotdash", color = "black") +  # 添加显著性阈值线
  geom_vline(xintercept = c(-1, 1), linetype = "dotdash", color = "black")  # 添加 Fold Change 阈值线
  theme(legend.position = "none")
print(volcano_plot)
dev.off()

significant_genes <- valid_genes %>%
  filter(significant == 'Significant') %>%
  arrange(desc(abs(fold_change)))

write.table(significant_genes, "DE_genes_from_TEs.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(gene_position_express, "all_gene_TEs_details.tsv", sep = "\t", row.names = FALSE, col.names = TRUE)
