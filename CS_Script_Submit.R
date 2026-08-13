# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(ggforce)      # geom_sina
library(patchwork)    # figure assembly
library(cowplot)      # plot_grid, get_legend
library(vegan)        # PCA support
library(factoextra)   # PCA visualisation
library(scales)       # rescale
library(pls)          # PLS regression
library(glmnet)       # ENET
library(stabs)        # stability selection
library(car)          # Levene's test
library(rstatix)      # KW + Dunn

# ── Output directory ──────────────────────────────────────────────────────────
#mix of output formats used for figures. set output directory here:
#out_dir <- 'x'

######################################################################################
# 0. DATA PREPARATION
######################################################################################

# ── Read data ─────────────────────────────────────────────────────────────────
#set working directory
#setwd("x")
#these are the two files that are provided in the repository
tbi.df <- read.csv("LitterDecayDat.csv")
pal.df <- read.csv("CoreDat.csv")

# ── Common variable integration ───────────────────────────────────────────────
common_vars   <- intersect(names(tbi.df), names(pal.df))
tbi.df_com    <- tbi.df %>% select(all_of(common_vars))
pal.df_com    <- pal.df %>% select(all_of(common_vars))

combined.df.com <- bind_rows(tbi.df_com, pal.df_com) %>%
  mutate(dataset = rep(c("tbi", "pal"), c(nrow(tbi.df), nrow(pal.df)))) %>%
  na.omit()

# ── Aesthetic maps (consistent across all figures) ────────────────────────────
group_colours <- c(
  "tbi_g8"      = "darkgreen",
  "tbi_g50"     = "darkgreen",
  "tbi_r8"      = "red",
  "tbi_r50"     = "red",
  "pal_0-30"    = "#9E3D22",
  "pal_30-100"  = "#9E3D22",
  "pal_100-200" = "#9E3D22",
  "pal_200-300" = "#9E3D22"
)

depth_shapes <- c(
  "tbi_g8"      = 17,
  "tbi_g50"     = 16,
  "tbi_r8"      = 17,
  "tbi_r50"     = 16,
  "pal_0-30"    = 17,
  "pal_30-100"  = 16,
  "pal_100-200" = 3,
  "pal_200-300" = 18
)

group_labels <- c(
  "tbi_g8"      = "Labile OM 8cm",
  "tbi_g50"     = "Labile OM 50cm",
  "tbi_r8"      = "Recalcitrant OM 8cm",
  "tbi_r50"     = "Recalcitrant OM 50cm",
  "pal_0-30"    = "Peat 0\u201330cm",
  "pal_30-100"  = "Peat 30\u2013100cm",
  "pal_100-200" = "Peat 100\u2013200cm",
  "pal_200-300" = "Peat 200\u2013300cm"
)

# ── Group + factor levels ─────────────────────────────────────────────────────
combined.df.com <- combined.df.com %>%
  mutate(
    group = case_when(
      grepl("tbi_g", code_all1) ~ "Labile OM",
      grepl("tbi_r", code_all1) ~ "Recalcitrant OM",
      grepl("pal",   code_all1) ~ "Peat"
    ),
    group     = factor(group, levels = c("Labile OM", "Recalcitrant OM", "Peat")),
    group3    = group,
    code_all1 = factor(code_all1, levels = names(group_colours))
  )

# ── Remove C organic outlier ──────────────────────────────────────────────────
# Visual inspection: one peat Corg value ~25% flagged as likely contamination
combined.df.com.clean <- combined.df.com %>%
  filter(!(X.Corg < 30 & dataset == "pal"))

removed <- anti_join(combined.df.com, combined.df.com.clean, by = names(combined.df.com))
cat("Outlier removed:", nrow(removed), "row(s)\n")
print(removed %>% select(code_all1, X.Corg, dataset))

# ── Variable sets ─────────────────────────────────────────────────────────────
var_labels <- c(
  X.Corg     = "A) Org. C (%)",
  X.N        = "B) N (%)",
  d15N       = "C) \u03b415N (\u2030)",
  d13C       = "D) \u03b413C (\u2030)",
  X.T        = "E) Light Trans. (%T)",
  carboArom  = "F) Carbo:Arom",
  carboAliph = "G) Carbo:Aliph",
  carboAcid  = "H) Carbo:Acid",
  carbo      = "I) Carbo (Abs)",
  aliph      = "J) Aliph (Abs)",
  acid       = "K) Acids (Abs)",
  arom       = "L) Arom (Abs)"
)
vars_to_plot <- names(var_labels)

# predictor set for ENET/PLS
pred_vars <- c("d15N", "d13C", "X.N", "X.Corg", "CN_ratio",
               "X.T", "acid", "arom", "carbo", "aliph",
               "carboArom", "carboAliph", "carboAcid")

custom_labels <- c("d15N", "d13C", "%N", "%Corg", "CN ratio", "%T",
                   "Acids", "Arom", "Carbo", "Aliph",
                   "CarboArom", "CarboAliph", "CarboAcid")

# ── Subsets ───────────────────────────────────────────────────────────────────
g8.all    <- tbi.df %>% filter(code_all1 == "tbi_g8")
g50.all   <- tbi.df %>% filter(code_all1 == "tbi_g50")
r8.all    <- tbi.df %>% filter(code_all1 == "tbi_r8")
r50.all   <- tbi.df %>% filter(code_all1 == "tbi_r50")
pal.dep30  <- pal.df %>% filter(code_all1 == "pal_0-30")
pal.dep100 <- pal.df %>% filter(code_all1 == "pal_30-100")
pal.dep200 <- pal.df %>% filter(code_all1 == "pal_100-200")
pal.dep300 <- pal.df %>% filter(code_all1 == "pal_200-300")

######################################################################################
#assumption tests: Kruskal-Wallis + Dunn post-hoc
######################################################################################
# ── 1. Simplified group variable (3 levels only) ──────────────────────────────

combined.df.com.clean <- combined.df.com.clean %>%
  mutate(
    group3 = case_when(
      grepl("tbi_g", code_all1) ~ "Labile OM",
      grepl("tbi_r", code_all1) ~ "Recalcitrant OM",
      grepl("pal",   code_all1) ~ "Peat"
    ),
    group3 = factor(group3, levels = c("Labile OM", "Recalcitrant OM", "Peat"))
  )

# ── 2. Assumption checks ──────────────────────────────────────────────────────
# Run for each variable — normality (Shapiro-Wilk) + homogeneity of variance (Levene)
# ANOVA is fairly robust to normality violations at n > 30, but Levene matters

assumption_check <- combined.df.com.clean %>%
  select(group3, all_of(vars_to_plot)) %>%
  pivot_longer(cols = all_of(vars_to_plot),
               names_to  = "variable",
               values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    # Shapiro-Wilk per group (only valid for n < 5000)
    shapiro_labile = shapiro.test(value[group3 == "Labile OM"])$p.value,
    shapiro_recalc = shapiro.test(value[group3 == "Recalcitrant OM"])$p.value,
    shapiro_peat   = shapiro.test(value[group3 == "Peat"])$p.value,
    # Levene's test across all 3 groups
    levene_p       = leveneTest(value ~ group3,
                                data = cur_data())$`Pr(>F)`[1],
    .groups = "drop"
  ) %>%
  mutate(
    normal_assumption = shapiro_labile > 0.05 &
      shapiro_recalc > 0.05 &
      shapiro_peat   > 0.05,
    equal_var         = levene_p > 0.05,
    recommended_test  = case_when(
      normal_assumption & equal_var  ~ "One-way ANOVA + Tukey HSD",
      normal_assumption & !equal_var ~ "Welch ANOVA + Games-Howell",
      TRUE                           ~ "Kruskal-Wallis + Dunn"
    )
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

print(assumption_check %>% select(variable, normal_assumption,
                                  equal_var, recommended_test), n = Inf)

# ── 3. Run all three tests — let assumption results guide which to report ──────

results_list <- vars_to_plot %>%
  set_names() %>%
  map(function(var) {
    
    dat <- combined.df.com.clean %>%
      select(group3, value = all_of(var)) %>%
      drop_na()
    
    list(
      # parametric: one-way ANOVA + Tukey
      anova     = aov(value ~ group3, data = dat) %>% summary(),
      tukey     = TukeyHSD(aov(value ~ group3, data = dat)),
      
      # parametric, unequal variance: Welch + Games-Howell
      welch     = oneway.test(value ~ group3, data = dat, var.equal = FALSE),
      gameshowell = dat %>% games_howell_test(value ~ group3),
      
      # non-parametric: Kruskal-Wallis + Dunn
      kruskal   = kruskal.test(value ~ group3, data = dat),
      dunn      = dat %>% dunn_test(value ~ group3, p.adjust.method = "bonferroni")
    )
  })

# ── 4. Tidy summary table — p-values from all three approaches ────────────────

pval_summary <- vars_to_plot %>%
  set_names() %>%
  map_dfr(function(var) {
    
    dat <- combined.df.com.clean %>%
      select(group3, value = all_of(var)) %>%
      drop_na()
    
    # overall test p-values
    anova_p   <- summary(aov(value ~ group3, data = dat))[[1]]$`Pr(>F)`[1]
    welch_p   <- oneway.test(value ~ group3, data = dat,
                             var.equal = FALSE)$p.value
    kruskal_p <- kruskal.test(value ~ group3, data = dat)$p.value
    
    # post-hoc: Dunn with Bonferroni (works regardless of which overall test used)
    dunn <- dat %>%
      dunn_test(value ~ group3, p.adjust.method = "bonferroni") %>%
      mutate(comparison = paste(group1, "vs", group2)) %>%
      select(comparison, p.adj, p.adj.signif)
    
    tibble(
      variable        = var,
      anova_p         = round(anova_p,   4),
      welch_p         = round(welch_p,   4),
      kruskal_p       = round(kruskal_p, 4),
      dunn_results    = list(dunn)
    )
  }, .id = NULL)

# print overall p-values
print(pval_summary %>% select(-dunn_results), n = Inf)

# print post-hoc comparisons for each variable
walk(vars_to_plot, function(var) {
  cat("\n──", var_labels[var], "──\n")
  pval_summary %>%
    filter(variable == var) %>%
    pull(dunn_results) %>%
    .[[1]] %>%
    print()
})

# ── 5. Export summary ─────────────────────────────────────────────────────────

# flat table of all Dunn post-hoc results
dunn_flat <- pval_summary %>%
  select(variable, dunn_results) %>%
  unnest(dunn_results) %>%
  mutate(variable = var_labels[variable])

write.csv(dunn_flat,  "stats_dunn_posthoc.csv",  row.names = FALSE)
write.csv(pval_summary %>% select(-dunn_results),
          "stats_overall_pvalues.csv", row.names = FALSE)

######################################################################################
#Figure: data visualization---------
######################################################################################
# ── 1. Aesthetic maps ─────────────────────────────────────────────────────────

group_colours <- c(
  "tbi_g8"      = "darkgreen",
  "tbi_g50"     = "darkgreen",
  "tbi_r8"      = "red",
  "tbi_r50"     = "red",
  "pal_0-30"    = "#9E3D22",
  "pal_30-100"  = "#9E3D22",
  "pal_100-200" = "#9E3D22",
  "pal_200-300" = "#9E3D22"
)

depth_shapes <- c(
  "tbi_g8"      = 17,
  "tbi_g50"     = 16,
  "tbi_r8"      = 17,
  "tbi_r50"     = 16,
  "pal_0-30"    = 17,
  "pal_30-100"  = 16,
  "pal_100-200" = 3,
  "pal_200-300" = 18
)

group_labels <- c(
  "tbi_g8"      = "Labile OM 8cm",
  "tbi_g50"     = "Labile OM 50cm",
  "tbi_r8"      = "Recalcitrant OM 8cm",
  "tbi_r50"     = "Recalcitrant OM 50cm",
  "pal_0-30"    = "Peat 0\u201330cm",
  "pal_30-100"  = "Peat 30\u2013100cm",
  "pal_100-200" = "Peat 100\u2013200cm",
  "pal_200-300" = "Peat 200\u2013300cm"
)

# ── 2. Group identifier + factor order ────────────────────────────────────────

combined.df.com <- combined.df.com %>%
  mutate(
    group = case_when(
      grepl("tbi_g", code_all1) ~ "Labile OM",
      grepl("tbi_r", code_all1) ~ "Recalcitrant OM",
      grepl("pal",   code_all1) ~ "Peat"
    ),
    group     = factor(group, levels = c("Labile OM", "Recalcitrant OM", "Peat")),
    code_all1 = factor(code_all1, levels = names(group_colours))
  )

# ── 3. Remove C organic outlier ───────────────────────────────────────────────

combined.df.com.clean <- combined.df.com %>%
  filter(!(X.Corg < 30 & dataset == "pal"))

removed <- anti_join(combined.df.com, combined.df.com.clean,
                     by = names(combined.df.com))
cat("Removed", nrow(removed), "row(s):\n")
print(removed %>% select(code_all1, X.Corg, dataset))

# ── 4. Variable labels — letter prefix kept short, units on second line ────────

var_labels <- c(
  X.Corg     = "A) Org. C (%)",
  X.N        = "B) N (%)",
  d15N       = "C) δ15N (\u2030)",
  d13C       = "D) δ13C (\u2030)",
  X.T        = "E) Light Trans. (%T)",
  carboArom  = "F) Carbo:Arom (Abs)",
  carboAliph = "G) Carbo:Aliph (Abs)",
  carboAcid  = "H) Carbo:Acid (Abs)",
  carbo      = "I) Carbo (Abs).",
  aliph      = "J) Aliph.(Abs)",
  acid       = "K) Acids (Abs)",
  arom       = "L) Arom. (Abs)"
)

# y-axis labels — full text lives here instead of strip
yaxis_labels <- c(
  X.Corg     = "Organic C (%)",
  X.N        = "N (%)",
  d13C       = "\u03b4\u00b9\u00b3C (\u2030)",
  d15N       = "\u03b4\u00b9\u2075N (\u2030)",
  X.T        = "Light transmission (%T)",
  carboArom  = "Carbohydrates:Aromatics",
  carboAliph = "Carbohydrates:Aliphatics",
  carboAcid  = "Carbohydrates:Acids",
  carbo      = "Carbohydrates (a.u.)",
  aliph      = "Aliphatics (a.u.)",
  acid       = "Acids (a.u.)",
  arom       = "Aromatics (a.u.)"
)

vars_to_plot <- names(var_labels)

# ── 5. Pivot to long ──────────────────────────────────────────────────────────

combined.long <- combined.df.com.clean %>%
  select(code_all1, group, all_of(vars_to_plot)) %>%
  pivot_longer(cols      = all_of(vars_to_plot),
               names_to  = "variable",
               values_to = "value") %>%
  mutate(
    var_strip = factor(variable, levels = vars_to_plot, labels = var_labels),
    var_ylab  = factor(variable, levels = vars_to_plot, labels = yaxis_labels),
    code_all1 = factor(code_all1, levels = names(group_colours))
  )

# ── 6. Plot ───────────────────────────────────────────────────────────────────

fig_explore <- ggplot(combined.long,
                      aes(x = group, y = value,
                          colour = code_all1, shape = code_all1)) +
  geom_sina(size = 1.8, alpha = 0.7, stroke = 0.3,
            maxwidth = 0.6, seed = 42) +
  scale_colour_manual(values = group_colours,
                      labels = group_labels,
                      name   = NULL) +
  scale_shape_manual(values  = depth_shapes,
                     labels  = group_labels,
                     name    = NULL) +
  facet_wrap(~ var_strip,
             scales = "free_y",
             ncol   = 4) +
  # y-axis label per facet via a blank labeller trick —
  # actual units shown in axis title area using tag
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1, size = 7.5),
    axis.text.y      = element_text(size = 7.5),
    axis.title.y     = element_text(size = 8, colour = "grey30"),
    strip.background = element_blank(),
    strip.text       = element_text(size = 8, face = "italic"),
    legend.position  = "bottom",
    legend.box       = "horizontal",
    legend.key.size  = unit(0.5, "cm"),
    legend.text      = element_text(size = 7.5),
    panel.spacing.x  = unit(0.8, "cm"),
    panel.spacing.y  = unit(0.7, "cm")
  ) +
  guides(colour = guide_legend(nrow = 2),
         shape  = guide_legend(nrow = 2))
fig_explore

#set output directory here
#setwd('x')
#also set file name
tiff(("x.tiff"), height = 18, width = 22, units = 'cm', compression = "lzw", res = 600)
fig_explore
dev.off()

######################################################################################
# FIGURE + TABLE 2: PCA----------
######################################################################################

# ── Prepare data ──────────────────────────────────────────────────────────────
identifier_vars  <- c("code_all", "code_all1", "dataset")
numeric_data     <- combined.df.com %>%
  select(-all_of(identifier_vars)) %>%
  select(where(is.numeric)) %>%
  na.omit()

numeric_data_norm <- as.data.frame(lapply(numeric_data, rescale))
colnames(numeric_data_norm) <- colnames(numeric_data)

pca_result <- prcomp(numeric_data_norm, center = TRUE, scale. = TRUE)

# match rows back to identifiers after na.omit
combined_matched <- combined.df.com %>%
  select(-all_of(identifier_vars)) %>%
  select(where(is.numeric)) %>%
  mutate(row_id = row_number()) %>%
  drop_na() %>%
  pull(row_id)

pca_scores           <- as.data.frame(pca_result$x)
pca_scores$code_all1 <- combined.df.com$code_all1[combined_matched]
pca_scores$dataset   <- combined.df.com$dataset[combined_matched]

# ── Symbol + colour maps ──────────────────────────────────────────────────────
symbol_map <- c("tbi_g8" = 17, "tbi_g50" = 16, "tbi_r8" = 17, "tbi_r50" = 16,
                "pal_0-30" = 17, "pal_30-100" = 16, "pal_100-200" = 3, "pal_200-300" = 18)
color_map  <- c("tbi_g8" = "darkgreen", "tbi_g50" = "darkgreen",
                "tbi_r8" = "red",       "tbi_r50" = "red",
                "pal_0-30" = "#9E3D22", "pal_30-100" = "#9E3D22",
                "pal_100-200" = "#9E3D22", "pal_200-300" = "#9E3D22")

# ── Figure 2: PCA biplot ──────────────────────────────────────────────────────
pca_var <- summary(pca_result)$importance
pc1_pct <- round(pca_var[2, 1] * 100, 0)
pc2_pct <- round(pca_var[2, 2] * 100, 0)

loadings    <- pca_result$rotation
scale_factor <- 9

#set file name
tiff(file.path(out_dir, "x.tiff"),
     width = 6, height = 6, units = "in", res = 400)
par(mfrow = c(1, 1), mar = c(4, 4, 1, 1), oma = c(1, 1, 1, 1))
plot(pca_scores$PC1, pca_scores$PC2,
     col = color_map[as.character(pca_scores$code_all1)],
     pch = symbol_map[as.character(pca_scores$code_all1)],
     xlab = paste0("PC1: ", pc1_pct, "%"),
     ylab = paste0("PC2: ", pc2_pct, "%"),
     main = " ",
     xlim = c(min(pca_scores$PC1) - 1.0, max(pca_scores$PC1) + 2.0),
     ylim = c(min(pca_scores$PC2) - 1.0, max(pca_scores$PC2) + 1.0))
arrows(0, 0, loadings[, 1] * scale_factor, loadings[, 2] * scale_factor,
       col = "black", length = 0.1)
text(loadings[, 1] * scale_factor, loadings[, 2] * scale_factor,
     labels = rownames(loadings), pos = 4, cex = 0.8)
legend("topright", legend = names(symbol_map),
       pch = symbol_map, col = unname(color_map),
       title = "Group", cex = 0.8)
dev.off()

# ── Table 2: PCA loadings ─────────────────────────────────────────────────────
loadings_df <- as.data.frame(round(loadings[, 1:4], 3))
loadings_df$variable <- rownames(loadings_df)
loadings_df <- loadings_df %>% select(variable, everything())

write.csv(loadings_df, file.path(out_dir, "Table2_PCA_loadings.csv"), row.names = FALSE)
cat("\nTable 2 (PCA loadings, first 4 PCs):\n")
print(loadings_df)

######################################################################################
# TABLES: ENET — performance + stability--------
######################################################################################

# ── Helper: run ENET + stabsel for one dataset ────────────────────────────────
run_enet <- function(df, response_var, label) {
  df_clean <- df %>% select(all_of(c(pred_vars, response_var))) %>% drop_na()
  x_mat    <- as.matrix(df_clean[, pred_vars])
  y_vec    <- df_clean[[response_var]]
  x_sc     <- scale(x_mat)
  y_sc     <- scale(y_vec)
  
  set.seed(123)
  cv_mod    <- cv.glmnet(x_sc, y_sc, alpha = 0.5, standardize = FALSE)
  best_lam  <- cv_mod$lambda.min
  best_mod  <- glmnet(x_sc, y_sc, alpha = 0.5, lambda = best_lam)
  y_pred    <- predict(best_mod, s = best_lam, newx = x_sc)
  
  sst  <- sum((y_sc - mean(y_sc))^2)
  sse  <- sum((y_pred - y_sc)^2)
  rsq  <- 1 - sse / sst
  mincv <- min(cv_mod$cvm)
  
  stab <- stabsel(x = x_sc, y = y_sc, fitfun = glmnet.lasso, cutoff = 0.6, PFER = 1)
  stab_probs <- stab$max
  
  list(
    label     = label,
    rsq       = round(rsq, 3),
    min_cv    = round(mincv, 3),
    coefs     = coef(best_mod),
    stab_probs = stab_probs
  )
}

# ── Run all subsets ───────────────────────────────────────────────────────────
enet_results <- list(
  g8   = run_enet(g8.all,    "ml_mean", "Labile OM 8cm"),
  g50  = run_enet(g50.all,   "ml_mean", "Labile OM 50cm"),
  r8   = run_enet(r8.all,    "ml_mean", "Recalcitrant OM 8cm"),
  r50  = run_enet(r50.all,   "ml_mean", "Recalcitrant OM 50cm"),
  dep30  = run_enet(pal.dep30,  "BD.", "Peat 0-30cm"),
  dep100 = run_enet(pal.dep100 %>% drop_na(), "BD.", "Peat 30-100cm"),
  dep200 = run_enet(pal.dep200 %>% drop_na(), "BD.", "Peat 100-200cm"),
  dep300 = run_enet(pal.dep300, "BD.", "Peat 200-300cm")
)

# ── Table 3: ENET performance ─────────────────────────────────────────────────
enet_perf <- map_dfr(enet_results, function(r) {
  tibble(group = r$label, R2 = r$rsq, min_CV_error = r$min_cv)
})
write.csv(enet_perf, file.path(out_dir, "Table3_ENET_performance.csv"), row.names = FALSE)
cat("\nTable 3 (ENET performance):\n"); print(enet_perf)

# ── Table 4: Stability selection probabilities ────────────────────────────────
stab_df <- map_dfr(enet_results, function(r) {
  probs <- r$stab_probs
  as_tibble(t(probs)) %>% mutate(group = r$label)
}) %>%
  select(group, everything()) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# add mean ± sd across groups
stab_df <- stab_df %>%
  bind_rows(
    stab_df %>%
      summarise(across(where(is.numeric), mean), group = "Mean") %>%
      bind_rows(
        stab_df %>%
          summarise(across(where(is.numeric), sd), group = "SD")
      )
  )

write.csv(stab_df, file.path(out_dir, "Table4_ENET_stability.csv"), row.names = FALSE)
cat("\nTable 4 (Stability selection probabilities):\n"); print(stab_df)

######################################################################################
# FIGURE + TABLES: PLS — %VAR, VIP, RMSEP-----
######################################################################################

# ── Helper: run PLS for one dataset ──────────────────────────────────────────
run_pls <- function(df, response_var, label) {
  df_clean <- df %>% select(all_of(c(pred_vars, response_var))) %>% drop_na()
  formula  <- as.formula(paste("cbind(",
                               paste(pred_vars, collapse = ", "),
                               ") ~", response_var))
  mod <- plsr(formula, data = df_clean, validation = "CV", segments = 10, scale = TRUE)
  vip <- varImp(mod, scale = FALSE)$Overall
  names(vip) <- pred_vars
  rmsep_val  <- RMSEP(mod)$val[1, 1, 2]   # 1 component CV RMSEP
  
  list(
    label    = label,
    mod      = mod,
    vip      = round(vip, 3),
    rmsep    = round(rmsep_val, 4)
  )
}

pls_results <- list(
  tbi    = run_pls(tbi.df,    "ml_mean", "TBI all"),
  g8     = run_pls(g8.all,    "ml_mean", "Labile OM 8cm"),
  g50    = run_pls(g50.all,   "ml_mean", "Labile OM 50cm"),
  r8     = run_pls(r8.all,    "ml_mean", "Recalcitrant OM 8cm"),
  r50    = run_pls(r50.all,   "ml_mean", "Recalcitrant OM 50cm"),
  pal    = run_pls(pal.df,    "BD.",      "Peat all"),
  dep30  = run_pls(pal.dep30, "BD.",      "Peat 0-30cm"),
  dep100 = run_pls(pal.dep100 %>% drop_na(), "BD.", "Peat 30-100cm"),
  dep200 = run_pls(pal.dep200 %>% drop_na(), "BD.", "Peat 100-200cm"),
  dep300 = run_pls(pal.dep300, "BD.",     "Peat 200-300cm")
)

# ── Table 5: VIP scores ───────────────────────────────────────────────────────
vip_table <- map_dfr(pls_results, function(r) {
  as_tibble(t(r$vip)) %>% mutate(group = r$label)
}) %>% select(group, all_of(pred_vars))

write.csv(vip_table, file.path(out_dir, "Table5_PLS_VIP.csv"), row.names = FALSE)
cat("\nTable 5 (PLS VIP scores):\n"); print(vip_table)

# ── Table 6: RMSEP ────────────────────────────────────────────────────────────
rmsep_table <- map_dfr(pls_results, function(r) {
  tibble(group = r$label, RMSEP_1comp = r$rmsep)
})
write.csv(rmsep_table, file.path(out_dir, "Table6_PLS_RMSEP.csv"), row.names = FALSE)
cat("\nTable 6 (RMSEP):\n"); print(rmsep_table)

# ── Figure 3: ENET diagnostic + PLS summary (4-panel) ────────────────────────
enet_diag  <- read.csv("enet_diag.csv")
df_stabprob <- read.csv("enet_stab.csv")
df_varvip2  <- read.csv("varvip.csv")

library(forcats)

enet_diagLong <- enet_diag %>%
  pivot_longer(-var, names_to = "model", values_to = "selected") %>%
  mutate(model = factor(model, levels = c("g8","g50","r8","r50",
                                          "pal_0.30","pal_30.100",
                                          "pal_100.200","pal_200.300")))

p3a <- ggplot(enet_diagLong, aes(x = model, y = fct_rev(var))) +
  geom_point(aes(shape = factor(selected), fill = factor(selected)), size = 4) +
  scale_shape_manual(values = c(21, 16)) +
  scale_fill_manual(values  = c("white", "black")) +
  theme_minimal(base_size = 12) +
  labs(title = "(A) ENET variable selection", x = "Model", y = "Variable",
       shape = "Selected", fill = "Selected") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

p4a <- ggplot(df_stabprob, aes(x = factor(var, levels = df_stabprob$var), y = avg)) +
  geom_col(fill = "black") +
  geom_errorbar(aes(ymin = pmax(0, avg - stdev), ymax = avg + stdev), width = 0.2) +
  #geom_hline(yintercept = 0.6, linetype = "dashed", colour = "red", linewidth = 0.5) +
  #annotate("text", x = 1, y = 0.62, label = "threshold = 0.6",
  #         hjust = 0, size = 3, colour = "red") +
  labs(title = "(B) ENET stability", x = "Variable", y = "Mean selection probability") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p5a <- ggplot(df_varvip2, aes(x = factor(X, levels = df_varvip2$X), y = VAR)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = 0, ymax = VAR + VAR_stdev), width = 0.2) +
  labs(title = "(C) PLS % variance explained", x = "Variable", y = "Variance explained (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p6a <- ggplot(df_varvip2, aes(x = factor(X, levels = df_varvip2$X), y = VIP)) +
  geom_col(fill = "firebrick") +
  geom_errorbar(aes(ymin = pmax(0, VIP - VIP_stdev), ymax = VIP + VIP_stdev), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black", linewidth = 0.5) +
  annotate("text", x = 1, y = 1.05, label = "VIP = 1", hjust = 0, size = 3) +
  labs(title = "(D) PLS VIP scores", x = "Variable", y = "VIP score") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#set file name
tiff(file.path(out_dir, "x.tiff"),
     width = 8, height = 8, units = "in", res = 400)
plot_grid(p3a, p4a, p5a, p6a, ncol = 2)
dev.off()

######################################################################################
# FIGURE: %T and d15N regression relationships (2x2 combined)-----
######################################################################################
#set output directory here
#setwd('x')
#also set file name
tiff(file = "x",
     width = 6, height = 8, units = "in", res = 400)
# 1. Use layout() to compress the 3rd row vertically
layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE), heights = c(1, 1, 0.45))

# Set tighter margins for main plot panels A-D
par(mar = c(3.5, 4, 1.8, 1), oma = c(1, 1, 1, 1), mgp = c(2.2, 0.7, 0))

# --- PANEL (A) ---
plot(g8.all$X.T, g8.all$ml_mean, col = "darkgreen", pch = 16, cex = 0.7,
     ylim = c(30,100), xlim = c(0,50), xlab = "% Light Transmission", ylab = "Mass Loss %")
abline(lm(g8.all$ml_mean ~ g8.all$X.T), col = "darkgreen", lwd = 2)

par(new = TRUE)
plot(g50.all$X.T, g50.all$ml_mean, col = "darkgreen", pch = 17, cex = 0.7,
     ylim = c(30,100), xlim = c(0,50), xlab = "", ylab = "", axes = FALSE)
abline(lm(g50.all$ml_mean ~ g50.all$X.T), col = "darkgreen", lwd = 2, lty = 2)

par(new = TRUE)
plot(r8.all$X.T, r8.all$ml_mean, col = "red", pch = 16, cex = 0.7,
     ylim = c(30,100), xlim = c(0,50), xlab = "", ylab = "", axes = FALSE)
abline(lm(r8.all$ml_mean ~ r8.all$X.T), col = "red", lwd = 2)

par(new = TRUE)
plot(r50.all$X.T, r50.all$ml_mean, col = "red", pch = 17, cex = 0.7,
     ylim = c(30,100), xlim = c(0,50), xlab = "", ylab = "", axes = FALSE)
abline(lm(r50.all$ml_mean ~ r50.all$X.T), col = "red", lwd = 2, lty = 2)
title("(A)", line = 0.5)

# --- PANEL (B) ---
plot(pal.dep30$X.T, pal.dep30$BD., col = "#9E3D22", pch = 16, cex = 0.7,
     ylim = c(0, 0.20), xlim = c(0, 5), xlab = "% Light Transmission", 
     ylab = expression(paste("Bulk Density (g/cm"^3, ")")))
abline(lm(pal.dep30$BD. ~ pal.dep30$X.T), col = "#9E3D22", lwd = 2)

par(new = TRUE)
plot(pal.dep100$X.T, pal.dep100$BD., col = "#9E3D22", pch = 17, cex = 0.7,
     ylim = c(0, 0.20), xlim = c(0, 5), xlab = "", ylab = "", axes = FALSE)
abline(lm(pal.dep100$BD. ~ pal.dep100$X.T), col = "#9E3D22", lwd = 2, lty = 2)

par(new = TRUE)
plot(pal.dep200$X.T, pal.dep200$BD., col = "#9E3D22", pch = 18, cex = 0.7,
     ylim = c(0, 0.20), xlim = c(0, 5), xlab = "", ylab = "", axes = FALSE)
abline(lm(pal.dep200$BD. ~ pal.dep200$X.T), col = "#9E3D22", lwd = 2, lty = 3)

par(new = TRUE)
plot(pal.dep300$X.T, pal.dep300$BD., col = "#9E3D22", pch = 15, cex = 0.7,
     ylim = c(0, 0.20), xlim = c(0, 5), xlab = "", ylab = "", axes = FALSE)
abline(lm(pal.dep300$BD. ~ pal.dep300$X.T), col = "#9E3D22", lwd = 2, lty = 4)
title("(B)", line = 0.5)

# --- PANEL (C) ---
plot(g8.all$d15N, g8.all$ml_mean, col = "darkgreen", pch = 16, cex = 0.7,
     ylim = c(0, 100), xlim = c(-1.5, 5.5), 
     xlab = expression(δ^15 * "N (‰)"), ylab = "Mass Loss %")
abline(lm(g8.all$ml_mean ~ g8.all$d15N), col = "darkgreen", lwd = 2)

par(new = TRUE)
plot(g50.all$d15N, g50.all$ml_mean, col = "darkgreen", pch = 17, cex = 0.7,
     ylim = c(0, 100), xlim = c(-1.5, 5.5), xlab = "", ylab = "", axes = FALSE)
abline(lm(g50.all$ml_mean ~ g50.all$d15N), col = "darkgreen", lwd = 2, lty = 2)

par(new = TRUE)
plot(r8.all$d15N, r8.all$ml_mean, col = "red", pch = 16, cex = 0.7,
     ylim = c(0, 100), xlim = c(-1.5, 5.5), xlab = "", ylab = "", axes = FALSE)
abline(lm(r8.all$ml_mean ~ r8.all$d15N), col = "red", lwd = 2)

par(new = TRUE)
plot(r50.all$d15N, r50.all$ml_mean, col = "red", pch = 17, cex = 0.7,
     ylim = c(0, 100), xlim = c(-1.5, 5.5), xlab = "", ylab = "", axes = FALSE)
abline(lm(r50.all$ml_mean ~ r50.all$d15N), col = "red", lwd = 2, lty = 2)
title("(C)", line = 0.5)

# --- PANEL (D) ---
plot(pal.dep30$d15N, pal.dep30$BD., col = "#9E3D22", pch = 16, cex = 0.7,
     ylim = c(0, 0.3), xlim = c(-1.5, 4), 
     xlab = expression(δ^15 * "N (‰)"), 
     ylab = expression(paste("Bulk Density (g/cm"^3, ")")))
abline(lm(pal.dep30$BD. ~ pal.dep30$d15N), col = "#9E3D22", lwd = 2)

par(new = TRUE)
plot(pal.dep100$d15N, pal.dep100$BD., col = "#9E3D22", pch = 17, cex = 0.7,
     ylim = c(0, 0.3), xlim = c(-1.5, 4), xlab = "", ylab = "", axes = FALSE)
abline(lm(pal.dep100$BD. ~ pal.dep100$d15N), col = "#9E3D22", lwd = 2, lty = 2)

par(new = TRUE)
plot(pal.dep200$d15N, pal.dep200$BD., col = "#9E3D22", pch = 18, cex = 0.7,
     ylim = c(0, 0.3), xlim = c(-1.5, 4), xlab = "", ylab = "", axes = FALSE)
abline(lm(pal.dep200$BD. ~ pal.dep200$d15N), col = "#9E3D22", lwd = 2, lty = 3)

par(new = TRUE)
plot(pal.dep300$d15N, pal.dep300$BD., col = "#9E3D22", pch = 15, cex = 0.7,
     ylim = c(0, 0.3), xlim = c(-1.5, 4), xlab = "", ylab = "", axes = FALSE)
abline(lm(pal.dep300$BD. ~ pal.dep300$d15N), col = "#9E3D22", lwd = 2, lty = 4)
title("(D)", line = 0.5)

# --- LEGENDS (ROW 3) ---
# Reset margins to 0 for legend panels so they take full space
par(mar = c(0, 0, 0, 0))

# Tea Bag Legend
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
legend("center", 
       title = "Tea Bags", 
       legend = c("Green tea bags buried at 8cm (tbi_g8)", 
                  "Green tea bags buried at 50cm (tbi_g50)", 
                  "Rooibos tea bags buried at 8cm (tbi_r8)", 
                  "Rooibos tea bags buried at 50cm (tbi_r50)"), 
       col = c("darkgreen", "darkgreen", "red", "red"), 
       lty = c(1, 2, 1, 2), 
       lwd = 3, 
       cex = 1.05,       # Larger text size
       title.cex = 1.1,
       seg.len = 2,      # Shorten sample lines to save space
       x.intersp = 0.8,  # Tighten text distance to line
       bty = "n")

# Peat Core Legend
plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
legend("center", 
       title = "Peat Cores", 
       legend = c("Peat core 0-30cm (pal_0-30)", 
                  "Peat core 30-100cm (pal_30-100)", 
                  "Peat core 100-200cm (pal_100-200)", 
                  "Peat core 200-300cm (pal_200-300)"), 
       col = "#9E3D22", 
       lty = 1:4, 
       lwd = 3, 
       cex = 1.05,       # Larger text size
       title.cex = 1.1,
       seg.len = 2, 
       x.intersp = 0.8,
       bty = "n")

# Reset layout
layout(1)

dev.off()


######################################################################################
# FIGURE: d15N vs CN ratio (combined dataset)-----
######################################################################################

log_model <- lm(d15N ~ log(X.N), data = combined.df.com)
x_seq <- seq(min(combined.df.com$X.N, na.rm = TRUE),
             max(combined.df.com$X.N, na.rm = TRUE), length.out = 200)
y_seq <- coef(log_model)[1] + coef(log_model)[2] * log(x_seq)
rsq_log <- round(summary(log_model)$r.squared, 3)

#set file name
tiff(file.path(out_dir, "x.tiff"),
     width = 7, height = 5, units = "in", res = 400)
par(mar = c(4, 4, 1, 1))
plot(combined.df.com$X.N, combined.df.com$d15N,
     col  = color_map[as.character(combined.df.com$code_all1)],
     pch  = symbol_map[as.character(combined.df.com$code_all1)],
     cex  = 1.2,
     xlab = "N (%)", ylab = "\u03b415N (\u2030)")
lines(x_seq, y_seq, col = "grey30", lwd = 2, lty = 2)
text(x = max(x_seq) * 0.6, y = min(y_seq) + 0.5,
     labels = paste0("R\u00b2 = ", rsq_log, "; p < 0.001"),
     cex = 0.9)
legend("bottomright", legend = names(symbol_map),
       pch = symbol_map, col = unname(color_map), title = "Group", cex = 0.75)
dev.off()

######################################################################################
# TABLE: Regression diagnostics — %T and d15N vs. mass loss / bulk density----
######################################################################################

run_lm_table <- function(x, y, group_label, xname, yname) {
  mod <- lm(y ~ x)
  s   <- summary(mod)
  tibble(
    group      = group_label,
    x_var      = xname,
    y_var      = yname,
    n          = length(x),
    intercept  = round(coef(mod)[1], 3),
    slope      = round(coef(mod)[2], 3),
    R2         = round(s$r.squared, 3),
    p_value    = round(s$coefficients[2, 4], 4)
  )
}

reg_table7 <- bind_rows(
  # %T vs mass loss
  run_lm_table(g8.all$X.T,  g8.all$ml_mean,  "Labile OM 8cm",         "%T", "Mass loss"),
  run_lm_table(g50.all$X.T, g50.all$ml_mean, "Labile OM 50cm",        "%T", "Mass loss"),
  run_lm_table(r8.all$X.T,  r8.all$ml_mean,  "Recalcitrant OM 8cm",   "%T", "Mass loss"),
  run_lm_table(r50.all$X.T, r50.all$ml_mean, "Recalcitrant OM 50cm",  "%T", "Mass loss"),
  # %T vs bulk density
  run_lm_table(pal.dep30$X.T,  pal.dep30$BD.,  "Peat 0-30cm",   "%T", "Bulk density"),
  run_lm_table(pal.dep100$X.T, pal.dep100$BD., "Peat 30-100cm", "%T", "Bulk density"),
  run_lm_table(pal.dep200$X.T, pal.dep200$BD., "Peat 100-200cm","%T", "Bulk density"),
  run_lm_table(pal.dep300$X.T, pal.dep300$BD., "Peat 200-300cm","%T", "Bulk density"),
  # d15N vs mass loss
  run_lm_table(g8.all$d15N,  g8.all$ml_mean,  "Labile OM 8cm",        "d15N", "Mass loss"),
  run_lm_table(g50.all$d15N, g50.all$ml_mean, "Labile OM 50cm",       "d15N", "Mass loss"),
  run_lm_table(r8.all$d15N,  r8.all$ml_mean,  "Recalcitrant OM 8cm",  "d15N", "Mass loss"),
  run_lm_table(r50.all$d15N, r50.all$ml_mean, "Recalcitrant OM 50cm", "d15N", "Mass loss"),
  # d15N vs bulk density
  run_lm_table(pal.dep30$d15N,  pal.dep30$BD.,  "Peat 0-30cm",   "d15N", "Bulk density"),
  run_lm_table(pal.dep100$d15N, pal.dep100$BD., "Peat 30-100cm", "d15N", "Bulk density"),
  run_lm_table(pal.dep200$d15N, pal.dep200$BD., "Peat 100-200cm","d15N", "Bulk density"),
  run_lm_table(pal.dep300$d15N, pal.dep300$BD., "Peat 200-300cm","d15N", "Bulk density")
)

write.csv(reg_table7, file.path(out_dir, "Table7_regression_pT_d15N.csv"), row.names = FALSE)
cat("\nTable 7 (Regression diagnostics, %T and d15N):\n"); print(reg_table7)

######################################################################################
# TABLE: Regression diagnostics — d15N, %N, %T inter-relationships-----
######################################################################################

datasets_all <- list(
  "Labile OM 8cm"        = g8.all,
  "Labile OM 50cm"       = g50.all,
  "Recalcitrant OM 8cm"  = r8.all,
  "Recalcitrant OM 50cm" = r50.all,
  "Peat 0-30cm"          = pal.dep30,
  "Peat 30-100cm"        = pal.dep100,
  "Peat 100-200cm"       = pal.dep200,
  "Peat 200-300cm"       = pal.dep300
)

reg_table8 <- map_dfr(names(datasets_all), function(nm) {
  d <- datasets_all[[nm]]
  bind_rows(
    run_lm_table(d$X.N,  d$d15N, nm, "%N",  "d15N"),
    run_lm_table(d$X.T,  d$d15N, nm, "%T",  "d15N"),
    run_lm_table(d$X.T,  d$X.N,  nm, "%T",  "%N")
  )
})

write.csv(reg_table8, file.path(out_dir, "Table8_regression_N_T_d15N.csv"), row.names = FALSE)
cat("\nTable 8 (N-humification relationships):\n"); print(reg_table8)

cat("\n\n=== All outputs saved to:", out_dir, "===\n")

######################################################################################
# FIGURE: Bulk density profile in SM-----------
######################################################################################

ggplot(pal.df, aes(x = BD., y = depth.x)) +
  geom_line()+
  geom_point() +
  scale_y_reverse()

ggplot(df, aes(x = BD., y = depth.x)) +
  geom_line(orientation = "y", linewidth = 0.8, color = "grey40") +
  geom_point(size = 2.5, color = "firebrick") +
  scale_y_reverse() +  # Puts depth 0 at the top and increases downward
  labs(
    x = "Bulk Density (g/cm³)",
    y = "Depth (cm)"
  ) +
  theme_bw()

######################################################################################
# END
######################################################################################
