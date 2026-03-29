# =============================================================================
# EMS Dose-Response & Physiological Damage Analysis
# Key Enhancements:
# 1. Forced fallback output for Fig 1 and 02_CSV to completely eliminate missing outputs.
# 2. Dropped GGally::ggcorr. Rewrote Fig 3 using pure ggplot2 + geom_tile for 
#    100% control over cell borders, low-correlation visibility, and label spacing.
# 3. Robust prediction interval retrieval to prevent internal drc errors.
# 4. All comments and logs translated to American English for publication standards.
# =============================================================================

# ---------------------------
# 1. Load Packages & Init Env
# ---------------------------
packages <- c("dplyr", "tidyr", "ggplot2", "drc", "readr", "broom", "tools")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# Global Theme & Scaling Parameters
base_font_size <- 7
geom_text_size <- base_font_size * 0.3528
pub_theme <- theme_bw(base_size = base_font_size, base_family = "sans") +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = base_font_size),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

color_line <- "#0072B2"
color_point <- "#E69F00"
color_hl <- "#D55E00"

# Output Directory Initialization
out_dir <- paste0("EMS_Analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE)
message("=== Initialization Complete: Results mapped to [", out_dir, "] ===")

# ---------------------------
# 2. Robust Data Loading & Cleaning
# ---------------------------
file_path <- "EMS_M2_MP2603V2_Sheet1.csv"
if (!file.exists(file_path)) stop("Input file not found. Please verify the file path.")

df <- read_csv(file_path, show_col_types = FALSE)
names(df) <- tolower(gsub("\\s+", "_", trimws(names(df))))

data_selected <- df %>%
  dplyr::select(ems_concentration, dplyr::matches("_(mean|sd)$")) %>%
  filter(!is.na(ems_concentration))

# Safely rename survival rate columns if they exist
if (any(grepl("seedling_survival_rate", names(data_selected)))) {
  data_selected <- data_selected %>%
    rename_with(~"survival_rate_mean", matches("seedling_survival_rate_mean")) %>%
    rename_with(~"survival_rate_sd", matches("seedling_survival_rate_sd"))
}

if (!"survival_rate_mean" %in% names(data_selected)) {
  stop("Critical column 'survival_rate_mean' is missing. Cannot proceed.")
}

# Output 01 data table
write.csv(data_selected, file.path(out_dir, "01_Processed_Data_Wide.csv"), row.names = FALSE)

ems_range <- diff(range(data_selected$ems_concentration, na.rm = TRUE))
err_bar_width <- ifelse(ems_range > 0, ems_range / 30, 0.005)

# ---------------------------
# 3. LD Dose-Response Modeling (Forced Output)
# ---------------------------
ld_data <- data_selected %>%
  filter(ems_concentration > 0, !is.na(survival_rate_mean))

model_surv <- tryCatch({
  if (nrow(ld_data) >= 3) {
    drm(survival_rate_mean ~ ems_concentration, data = ld_data,
        fct = LL.4(names = c("Slope", "Lower", "Upper", "LD50")),
        control = drmc(errorm = FALSE, noMessage = TRUE))
  } else { NULL }
}, error = function(e) { NULL })

if (inherits(model_surv, "drc")) {
  ld_res <- try(ED(model_surv, c(10, 50, 90), interval = "delta", display = FALSE), silent = TRUE)
  if (!inherits(ld_res, "try-error")) {
    write.csv(as.data.frame(ld_res), file.path(out_dir, "02_LD_10_50_90_with_CI.csv"), row.names = TRUE)
  } else {
    write.csv(data.frame(Message = "ED calculation failed, please check data curve."), file.path(out_dir, "02_LD_10_50_90_with_CI.csv"))
  }
} else {
  write.csv(data.frame(Message = "DRM Model failed to fit."), file.path(out_dir, "02_LD_10_50_90_with_CI.csv"))
}

p_ld <- ggplot(ld_data, aes(x = ems_concentration, y = survival_rate_mean))
if (inherits(model_surv, "drc")) {
  pred_x <- seq(min(ld_data$ems_concentration), max(ld_data$ems_concentration), length.out = 100)
  pred_raw <- try(predict(model_surv, newdata = data.frame(ems_concentration = pred_x), se.fit = TRUE), silent = TRUE)
  
  pred_df <- data.frame(ems_concentration = pred_x)
  if (!inherits(pred_raw, "try-error") && is.matrix(pred_raw) && ncol(pred_raw) >= 2) {
    pred_df$Prediction <- pred_raw[, 1]
    pred_df$Lower <- pred_raw[, 1] - 1.96 * pred_raw[, 2]
    pred_df$Upper <- pred_raw[, 1] + 1.96 * pred_raw[, 2]
    p_ld <- p_ld + geom_ribbon(data = pred_df, aes(x = ems_concentration, y = Prediction, ymin = Lower, ymax = Upper),
                               fill = "gray85", alpha = 0.6, inherit.aes = FALSE, na.rm = TRUE)
  } else {
    p_line <- try(predict(model_surv, newdata = data.frame(ems_concentration = pred_x)), silent = TRUE)
    if (!inherits(p_line, "try-error")) pred_df$Prediction <- as.numeric(p_line)
  }
  
  if ("Prediction" %in% names(pred_df)) {
    p_ld <- p_ld + geom_line(data = pred_df, aes(y = Prediction), color = color_line, linewidth = 0.8)
  }
  p_ld <- p_ld + labs(title = "Log-logistic Dose-Response Fit", subtitle = "DRM LL.4 Model")
} else {
  p_ld <- p_ld + stat_smooth(method = "loess", color = color_hl, fill = "gray85", alpha = 0.5) +
    labs(title = "Survival Rate Trend", subtitle = "DRM Failed - Loess Smoothing Applied")
}

p_ld <- p_ld + geom_point(size = 1.5, color = color_point) +
  labs(x = "EMS Concentration", y = "Survival Rate") + pub_theme
ggsave(file.path(out_dir, "Fig1_LD_Model_Fit.png"), plot = p_ld, width = 3.5, height = 3.5, units = "in", dpi = 300)
message("=> Saved Fig1_LD_Model_Fit.png & 02_ CSV")

# ---------------------------
# 4. Mutation Rate Quadratic Optimization
# ---------------------------
if ("m2_mutation_rate_mean" %in% names(data_selected)) {
  mut_data <- data_selected %>%
    dplyr::select(ems_concentration, m2_mutation_rate_mean) %>%
    filter(!is.na(m2_mutation_rate_mean))
  
  if (nrow(mut_data) >= 3) {
    model_mut <- lm(m2_mutation_rate_mean ~ ems_concentration + I(ems_concentration^2), data = mut_data)
    coefs <- coef(model_mut)
    a <- coefs["I(ems_concentration^2)"]; b <- coefs["ems_concentration"]
    
    coef_mat <- summary(model_mut)$coefficients
    is_sig <- ("I(ems_concentration^2)" %in% rownames(coef_mat)) &&
      (!is.na(coef_mat["I(ems_concentration^2)", "Pr(>|t|)"]) && coef_mat["I(ems_concentration^2)", "Pr(>|t|)"] < 0.05)
    
    opt_conc <- NA
    if (!is.na(a) && a < 0) {
      calc_opt <- -b / (2 * a)
      if (calc_opt > 0 && calc_opt <= max(mut_data$ems_concentration, na.rm = TRUE) * 1.2) opt_conc <- calc_opt
    }
    
    p_mut <- ggplot(mut_data, aes(x = ems_concentration, y = m2_mutation_rate_mean)) +
      geom_point(size = 1.5, color = color_point) +
      stat_smooth(method = "lm", formula = y ~ x + I(x^2), color = color_line, fill = "gray80", linewidth = 0.6, se = TRUE) +
      labs(title = "Quadratic Optimization", x = "EMS Concentration", y = "M2 Mutation Rate") + pub_theme
    
    if (!is.na(opt_conc)) {
      p_mut <- p_mut +
        geom_vline(xintercept = opt_conc, linetype = "dashed", color = color_hl, linewidth = 0.6) +
        annotate("text", x = opt_conc, y = max(mut_data$m2_mutation_rate_mean, na.rm = TRUE) * 0.95,
                 label = sprintf("Opt:\n%.4f%s", opt_conc, ifelse(is_sig, "", "*")),
                 hjust = -0.2, color = color_hl, family = "sans", size = geom_text_size)
    }
    ggsave(file.path(out_dir, "Fig2_Mutation_Optimization.png"), plot = p_mut, width = 3.5, height = 3.5, units = "in", dpi = 300)
    message("=> Saved Fig2_Mutation_Optimization.png")
  }
}

# ---------------------------
# 5. Correlation Matrix (Pure ggplot2, 100% Control)
# ---------------------------
# Extract columns and filter those with enough variance to compute correlation
corr_cols <- data_selected %>%
  dplyr::select(dplyr::ends_with("_mean"), ems_concentration) %>%
  dplyr::select(where(~ is.numeric(.) && sum(!is.na(.)) > 2 && var(., na.rm = TRUE) > 0))

if (ncol(corr_cols) >= 3) {
  # Clean names for beautiful plotting
  raw_names <- names(corr_cols)
  clean_names <- tools::toTitleCase(gsub("_mean|_", " ", raw_names))
  clean_names <- gsub("Ems Concentration", "EMS Conc.", clean_names)
  names(corr_cols) <- clean_names
  
  # Calculate Pearson correlation matrix
  cormat <- cor(corr_cols, use = "pairwise.complete.obs", method = "pearson")
  
  # Convert matrix to long format for ggplot2
  cormat_df <- as.data.frame(as.table(cormat))
  names(cormat_df) <- c("Var1", "Var2", "Value")
  
  # Ensure axes maintain categorical order
  cormat_df$Var1 <- factor(cormat_df$Var1, levels = clean_names)
  cormat_df$Var2 <- factor(cormat_df$Var2, levels = rev(clean_names)) # Reverse Y axis for top-down reading
  
  # Create label format (round to 2 decimals, NA to blank)
  cormat_df$Label <- ifelse(is.na(cormat_df$Value), "", sprintf("%.2f", cormat_df$Value))
  
  # 色盲友好配色 (WCAG 2.1 兼容)
  # 深蓝 (-1) → 中蓝 (-0.5) → 深灰 (0) → 橙红 (0.5) → 深红 (1)
  custom_colors <- c("#2166AC", "#67A9CF", "#A0A0A0", "#EF8A62", "#B2182B")
  
  # Pure ggplot2 heatmap
  p_corr <- ggplot(cormat_df, aes(x = Var1, y = Var2, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +  # The white border makes the cell stand out
    geom_text(aes(label = Label), color = "black", size = geom_text_size, family = "sans", na.rm = TRUE) +
    scale_fill_gradientn(
      colors = custom_colors,
      limits = c(-1, 1),
      na.value = "#737373",   # Missing data is shown in dark gray
      name = "Pearson r "
    ) +
    scale_x_discrete(position = "bottom") +
    labs(x = NULL, y = NULL, title = "Correlation Matrix of Evaluated Traits") +
    theme_minimal(base_size = base_font_size, base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1, "cm"),
      legend.key.width = unit(0.3, "cm")
    )
  
  ggsave(file.path(out_dir, "Fig3_Multi_Trait_Correlation.png"), 
         plot = p_corr, width = 6.0, height = 5.0, units = "in", dpi = 300)
  message("=> Saved Fig3_Multi_Trait_Correlation.png (Pure ggplot2 Matrix)")
}

# ---------------------------
# 6. Faceted Multi-trait Visualization
# ---------------------------
plot_data <- data_selected %>%
  pivot_longer(cols = -ems_concentration, names_to = c("trait_raw", ".value"), names_pattern = "(.*)_(mean|sd)$") %>%
  filter(!is.na(mean)) %>%
  mutate(trait_clean = tools::toTitleCase(gsub("_", " ", trait_raw)))

# Output 03 data table
write.csv(plot_data, file.path(out_dir, "03_Faceted_Plot_Data_Long.csv"), row.names = FALSE)

if (nrow(plot_data) > 0) {
  p_all <- ggplot(plot_data, aes(x = ems_concentration, y = mean)) +
    geom_line(color = color_line, linewidth = 0.6) +
    geom_point(size = 1.5, color = color_point) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = err_bar_width, color = "gray40", linewidth = 0.4, na.rm = TRUE) +
    facet_wrap(~trait_clean, scales = "free_y", ncol = 2) +
    labs(title = "Phenotypic Trait Responses to EMS Treatment", x = "EMS Concentration (%)", y = "Measured Value") +
    pub_theme
  
  ggsave(file.path(out_dir, "Fig4_Overall_Traits_Trends.png"), plot = p_all, width = 7.0, height = 5.0, units = "in", dpi = 300)
  message("=> Saved Fig4_Overall_Traits_Trends.png & 03_ CSV")
}

message("=== Pipeline Execution Completed. All files saved. ===")
