# =============================================================================
# EMS Mutation Spectrum & Rigorous Enrichment Framework 
#
# Key Upgrades from Critical Review:
# 1. Robust Column Mapping & Frequency Parsing (Regex-based).
# 2. Pre-flight Audit: Logs NAs, consistent N-sizes, filters Line_Total < 5.
# 3. Exact Binomial Test (vs Global Mean) with FDR correction instead of GLM.
# 4. Standardized Residuals mapped to P-values & FDR corrected before highlighting.
# 5. Dynamic plot sizing and overlapping text prevention using /.pt scalar.
# 6. Session Info export for full reproducibility.
# =============================================================================

# ---------------------------
# 1. Load Packages & Init
# ---------------------------
packages <- c("dplyr", "tidyr", "ggplot2", "readr", "forcats", "ggsci", "scales", "stringr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ggplot2 font sizing standard: pt to mm conversion factor is internal `.pt`
base_font_size <- 7
geom_text_size <- base_font_size / .pt 

pub_theme <- theme_bw(base_size = base_font_size, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold")
  )

out_dir <- paste0("Mutation_Analysis_M2_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE)

# Start Audit Log
audit_log_path <- file.path(out_dir, "00_Data_Audit_Log.txt")
sink(audit_log_path)
cat("=== EMS Data Pre-flight Audit Log ===\n")
cat("Run Time:", as.character(Sys.time()), "\n\n")

# ---------------------------
# 2. Bulletproof Data Loading & Audit
# ---------------------------
file_path <- "EMS_M2_MP2603V2_Sheet3.csv"
if (!file.exists(file_path)) stop("Input file not found.")

df_raw <- read_csv(file_path, show_col_types = FALSE)

# Report Initial NAs
na_count <- sum(is.na(df_raw))
cat(sprintf("[AUDIT] Detected %d missing values in the raw dataset.\n", na_count))

# Robust Column Mapping via Regex Keywords
safe_rename <- function(df) {
  orig_cols <- colnames(df)
  lower_cols <- tolower(orig_cols)
  
  # 匹配各列
  c_cat   <- orig_cols[str_detect(lower_cols, "category|phenotype")]
  c_desc  <- orig_cols[str_detect(lower_cols, "desc")]
  c_ofreq <- orig_cols[str_detect(lower_cols, "overall.*freq")]
  c_line  <- orig_cols[str_detect(lower_cols, "plant.*line|line.*plant")]
  c_mut   <- orig_cols[str_detect(lower_cols, "mutant.*in")]
  c_tot   <- orig_cols[str_detect(lower_cols, "line.*total")]
  c_lfreq <- orig_cols[str_detect(lower_cols, "freq.*in.*line")]
  
  # 必须存在的列检查
  if (length(c_cat) == 0) stop("Cannot find Category column (looking for 'category' or 'phenotype')")
  if (length(c_line) == 0) stop("Cannot find Plant_Line column (looking for 'plant' and 'line')")
  if (length(c_mut) == 0) stop("Cannot find Mutants_in_Line column (looking for 'mutant' and 'in')")
  if (length(c_tot) == 0) stop("Cannot find Line_Total column (looking for 'line' and 'total')")
  
  # 先重命名必须列
  df <- df %>%
    rename(
      Category = !!sym(c_cat[1]),
      Plant_Line = !!sym(c_line[1]),
      Mutants_in_Line = !!sym(c_mut[1]),
      Line_Total = !!sym(c_tot[1])
    )
  
  # 可选列分别处理
  if (length(c_desc) > 0) {
    df <- df %>% rename(Original_Desc = !!sym(c_desc[1]))
  } else {
    df$Original_Desc <- NA  # 或根据需求添加空列
    warning("Description column not found, added as NA")
  }
  
  if (length(c_ofreq) > 0) {
    df <- df %>% rename(Overall_Freq_Raw = !!sym(c_ofreq[1]))
  } else {
    df$Overall_Freq_Raw <- NA
    warning("Overall_Freq_Raw column not found, added as NA")
  }
  
  if (length(c_lfreq) > 0) {
    df <- df %>% rename(Freq_in_Line_Raw = !!sym(c_lfreq[1]))
  } else {
    df$Freq_in_Line_Raw <- NA
    warning("Freq_in_Line_Raw column not found, added as NA")
  }
  
  return(df)
}

df_clean <- safe_rename(df_raw) %>% filter(!is.na(Category) & !is.na(Plant_Line))

# Robust Mixed-Format Frequency Parser
parse_mixed_freq <- function(x) {
  x_str <- as.character(x)
  has_pct <- grepl("%", x_str)
  num <- readr::parse_number(x_str)
  # If has % or is unexpectedly > 1 (e.g., '7.83' meaning 7.83%), convert to decimal
  ifelse(!is.na(num) & (has_pct | num > 1), num / 100, num)
}

df_clean <- df_clean %>%
  mutate(
    Overall_Freq = parse_mixed_freq(Overall_Freq_Raw),
    Freq_in_Line = parse_mixed_freq(Freq_in_Line_Raw),
    Plant_Line = as.character(Plant_Line)
  )

# Audit: Filter out statistically unstable low-N lines
MIN_LINE_TOTAL <- 10  # 原为 5，现提高至 10
dropped_lines <- df_clean %>% filter(Line_Total < MIN_LINE_TOTAL) %>% pull(Plant_Line) %>% unique()
if(length(dropped_lines) > 0) {
  cat(sprintf("[AUDIT] Dropped %d Plant Lines due to Line_Total < %d: %s\n", 
              length(dropped_lines), MIN_LINE_TOTAL, paste(dropped_lines, collapse = ", ")))
}
df_clean <- df_clean %>% filter(Line_Total >= MIN_LINE_TOTAL)

# Audit: Check Line_Total consistency
line_inconsistencies <- df_clean %>% group_by(Plant_Line) %>% summarise(n_vals = n_distinct(Line_Total)) %>% filter(n_vals > 1)
if(nrow(line_inconsistencies) > 0) cat("[AUDIT WARNING] Found inconsistent Line_Total within the same Plant_Line.\n")

line_n_map <- df_clean %>%
  group_by(Plant_Line) %>%
  summarise(N = max(Line_Total, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Line_Label = paste0(Plant_Line, "\n(n=", N, ")"))

df_clean <- df_clean %>% left_join(line_n_map, by = "Plant_Line")

# Audit: Calculate consensus Overall_Freq to avoid max() distortion
overall_stats <- df_clean %>%
  group_by(Category) %>%
  summarise(Overall_Freq = first(na.omit(Overall_Freq)), .groups = 'drop') %>%
  arrange(desc(Overall_Freq)) %>%
  mutate(Category = fct_reorder(Category, Overall_Freq))

cat("=== Audit Complete ===\n")
sink()

# ---------------------------
# 3. Overall Spectrum (Lollipop Plot)
# ---------------------------
p_lollipop <- ggplot(overall_stats, aes(x = Category, y = Overall_Freq)) +
  geom_segment(aes(x = Category, xend = Category, y = 0, yend = Overall_Freq), color = "gray50", linewidth = 0.8) +
  geom_point(size = 2.5, color = "#D55E00") +
  geom_text(aes(label = scales::percent(Overall_Freq, accuracy = 0.01)), 
            hjust = -0.3, size = geom_text_size, family = "sans") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.2))) +
  coord_flip() +
  labs(title = "Phenotypic Mutation Spectrum", x = "Mutation Category", y = "Global Frequency") +
  pub_theme + theme(panel.grid.major.y = element_blank())

ggsave(file.path(out_dir, "Fig1_Spectrum_Lollipop.png"), plot = p_lollipop, width = 5.0, height = 4.0, units = "in", dpi = 300)

# ---------------------------
# 4. Statistical Validation (FDR-Corrected Exact Binomial & Residuals)
# ---------------------------
# 4A. Exact Binomial Test for Line Mutation Burden
# Calculates global baseline probability safely
total_population <- sum(line_n_map$N, na.rm = TRUE)
total_mutants <- sum(df_clean$Mutants_in_Line, na.rm = TRUE)
global_p <- total_mutants / total_population

binom_results <- df_clean %>%
  group_by(Plant_Line) %>%
  summarise(
    Total_Mutants = sum(Mutants_in_Line, na.rm = TRUE),
    Line_Total = max(Line_Total, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  rowwise() %>%
  mutate(
    Binom_P = binom.test(Total_Mutants, Line_Total, p = global_p)$p.value,
    Line_Freq = Total_Mutants / Line_Total
  ) %>%
  ungroup() %>%
  mutate(
    FDR_Adj_P = p.adjust(Binom_P, method = "BH"),
    Significance = case_when(FDR_Adj_P < 0.01 ~ "**", FDR_Adj_P < 0.05 ~ "*", TRUE ~ "ns")
  )
write.csv(binom_results, file.path(out_dir, "Tab1_Line_Burden_Binomial_FDR.csv"), row.names = FALSE)

# 4B. Category-Specific Binomial Enrichment Test
# Calculate the global baseline probability for each phenotypic category 
# and examine whether the number of mutations in each line within that category is significantly higher than expected.
category_global_stats <- df_clean %>%
  group_by(Category) %>%
  summarise(
    Total_Mutants_in_Category = sum(Mutants_in_Line, na.rm = TRUE),
    Total_Individuals = sum(Line_Total, na.rm = TRUE),
    Global_Freq = Total_Mutants_in_Category / Total_Individuals,
    .groups = 'drop'
  )

category_specific_binom <- df_clean %>%
  group_by(Plant_Line, Category) %>%
  summarise(
    Observed = sum(Mutants_in_Line, na.rm = TRUE),
    Line_Total = max(Line_Total, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  left_join(category_global_stats, by = "Category") %>%
  rowwise() %>%
  mutate(
    # One-tailed test: Whether the observed value is significantly higher than the global expectation.
    Binom_P = binom.test(Observed, Line_Total, p = Global_Freq, alternative = "greater")$p.value,
    Observed_Freq = Observed / Line_Total
  ) %>%
  ungroup() %>%
  mutate(
    FDR_Adj_P = p.adjust(Binom_P, method = "BH"),
    Significance = case_when(
      FDR_Adj_P < 0.01 ~ "**",
      FDR_Adj_P < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )
write.csv(category_specific_binom, file.path(out_dir, "Tab3_Category_Specific_Binomial_FDR.csv"), row.names = FALSE)

# 4C. FDR-Corrected Standardized Residuals for Specific Hotspots
mat_df <- df_clean %>%
  group_by(Plant_Line, Category) %>%
  summarise(Count = sum(Mutants_in_Line, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Category, values_from = Count, values_fill = 0)

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$Plant_Line

# Safe matrix degradation
mat <- mat[rowSums(mat) > 0, colSums(mat) > 0, drop = FALSE]

set.seed(2026)
chisq_mc <- chisq.test(mat, simulate.p.value = TRUE, B = 10000)

# Convert StdRes to p-values and apply FDR
stdres_long <- as.data.frame(as.table(chisq_mc$stdres)) %>%
  rename(Plant_Line = Var1, Category = Var2, StdRes = Freq) %>%
  mutate(
    P_val = 2 * pnorm(abs(StdRes), lower.tail = FALSE),
    FDR_P = p.adjust(P_val, method = "BH"),
    # Only assign * if FDR-adjusted p-value is significant AND StdRes is positive (enriched)
    Significance = case_when(
      FDR_P < 0.01 & StdRes > 0 ~ "**",
      FDR_P < 0.05 & StdRes > 0 ~ "*",
      TRUE ~ ""
    ),
    Plant_Line = as.character(Plant_Line), Category = as.character(Category)
  )
write.csv(stdres_long, file.path(out_dir, "Tab2_Hotspots_StdRes_FDR.csv"), row.names = FALSE)

# ---------------------------
# 5. Enrichment Heatmap (FDR Corrected)
# ---------------------------
heatmap_data <- df_clean %>%
  left_join(stdres_long, by = c("Plant_Line", "Category")) %>%
  mutate(
    Line_Label = fct_reorder(Line_Label, Mutants_in_Line, .fun = sum, .desc = TRUE),
    Category = factor(Category, levels = levels(overall_stats$Category))
  )

# Dynamic width formulation based on N of lines
dyn_width <- max(6.0, length(unique(heatmap_data$Plant_Line)) * 0.25 + 2)

p_heatmap <- ggplot(heatmap_data, aes(x = Line_Label, y = Category)) +
  geom_tile(aes(fill = Freq_in_Line), color = "white", linewidth = 0.3) +
  geom_text(aes(label = Significance), color = "white", size = geom_text_size * 1.5, vjust = 0.7) +
  scale_fill_viridis_c(option = "magma", labels = scales::percent_format(), 
                       na.value = "gray90", direction = -1, name = "Local Freq") +
  labs(title = "Mutation Hotspots & Enrichment",
       subtitle = "Color = Local Frequency; Asterisks (*) = FDR-corrected significant enrichment (q < 0.05)",
       x = "Plant Line (Total Sample Size)", y = "Phenotype Category") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(), panel.border = element_blank())

ggsave(file.path(out_dir, "Fig2_Enrichment_Heatmap_FDR.png"), plot = p_heatmap, width = dyn_width, height = 4.5, units = "in", dpi = 300)

# ---------------------------
# 6. Proportional Stacked Bar Chart (Overlap Defense)
# ---------------------------
stacked_data <- heatmap_data %>% filter(Mutants_in_Line > 0)

num_colors <- length(unique(stacked_data$Category))
qual_palette <- colorRampPalette(ggsci::pal_npg("nrc")(min(num_colors, 10)))(num_colors)

p_stacked <- ggplot(stacked_data, aes(x = Line_Label, y = Mutants_in_Line, fill = Category)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
  # Overlap Defense: Only show text if the proportion is > 5% of the line's total mutants
  geom_text(aes(label = ifelse(Mutants_in_Line / sum(Mutants_in_Line) > 0.05, Mutants_in_Line, "")), 
            position = position_fill(vjust = 0.5), 
            color = "white", size = geom_text_size, family = "sans", fontface = "bold") +
  scale_fill_manual(values = qual_palette) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Proportional Diversity of Mutations per Line",
       subtitle = "Labels inside bars represent absolute mutant counts (hidden if proportion < 5%)",
       x = "Plant Line (Total Sample Size)", y = "Proportion of Total Mutations") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right", legend.key.size = unit(0.3, "cm"))

ggsave(file.path(out_dir, "Fig3_Diversity_StackedBar.png"), plot = p_stacked, width = dyn_width, height = 4.5, units = "in", dpi = 300)

# ---------------------------
# 7. Reproducibility Log
# ---------------------------
sink(file.path(out_dir, "00_Session_Info.txt"))
print(sessionInfo())
sink()

message("=== Pipeline Execution Completed. All strict validation logs and plots saved. ===")
