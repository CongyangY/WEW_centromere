suppressPackageStartupMessages({
  library(stats)
})

# =========================
# 路径设置
# =========================
input_dir <- "../../domesication/Langdon/cluster_results/tables"
output_file <- file.path(input_dir, "all_inversion_fisher_tests.tsv")

# 只读 sample assignment 文件
files <- list.files(
  input_dir,
  pattern = "_sample_assignment\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("没有找到 *_sample_assignment.tsv 文件")
}

# =========================
# 辅助函数
# =========================

# 安全版 Fisher 检验
safe_fisher <- function(mat) {
  out <- list(
    p.value = NA,
    odds.ratio = NA
  )
  
  # 必须是 2x2 才能提 OR
  if (nrow(mat) < 2 || ncol(mat) < 2) return(out)
  
  ft <- tryCatch(
    fisher.test(mat),
    error = function(e) NULL
  )
  
  if (is.null(ft)) return(out)
  
  out$p.value <- ft$p.value
  
  # 只有 2x2 时有 odds ratio
  if (!is.null(ft$estimate)) {
    out$odds.ratio <- unname(ft$estimate)
  }
  
  return(out)
}

# 定义 cultivated-like cluster：
# 取 (DEM + DUR) 在该 cluster 中比例最高的那个
get_cultivated_like_cluster <- function(df) {
  tab <- table(df$cluster, df$group)
  
  # 保证列存在
  for (g in c("DEM", "DUR")) {
    if (!g %in% colnames(tab)) {
      tab <- cbind(tab, setNames(rep(0, nrow(tab)), g))
    }
  }
  
  tab <- tab[, sort(colnames(tab)), drop = FALSE]
  row_sum <- rowSums(tab)
  cultivated_score <- (tab[, "DEM"] + tab[, "DUR"]) / ifelse(row_sum == 0, NA, row_sum)
  
  cultivated_like_cluster <- names(which.max(cultivated_score))
  return(cultivated_like_cluster)
}

# 2x2 检验：某个 groupA vs groupB，在 cultivated-like cluster 上是否有富集差异
run_binary_fisher <- function(df, cultivated_like_cluster, groupA, groupB) {
  sub <- df[df$group %in% c(groupA, groupB), , drop = FALSE]
  
  if (nrow(sub) == 0) {
    return(list(p.value = NA, odds.ratio = NA, nA = 0, nB = 0))
  }
  
  sub$group2 <- factor(sub$group, levels = c(groupA, groupB))
  sub$is_cultivated_like <- factor(
    ifelse(sub$cluster == cultivated_like_cluster, "yes", "no"),
    levels = c("yes", "no")
  )
  
  mat <- table(sub$group2, sub$is_cultivated_like)
  
  # 如果两组里某一组完全没有样本，没法做有效检验
  if (sum(mat[groupA, ]) == 0 || sum(mat[groupB, ]) == 0) {
    return(list(
      p.value = NA,
      odds.ratio = NA,
      nA = sum(sub$group == groupA),
      nB = sum(sub$group == groupB)
    ))
  }
  
  ft <- tryCatch(
    fisher.test(mat),
    error = function(e) NULL
  )
  
  if (is.null(ft)) {
    return(list(
      p.value = NA,
      odds.ratio = NA,
      nA = sum(sub$group == groupA),
      nB = sum(sub$group == groupB)
    ))
  }
  
  return(list(
    p.value = ft$p.value,
    odds.ratio = if (!is.null(ft$estimate)) unname(ft$estimate) else NA,
    nA = sum(sub$group == groupA),
    nB = sum(sub$group == groupB)
  ))
}

# =========================
# 主循环
# =========================
res_list <- list()

for (f in files) {
  dat <- read.table(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # 列名标准化
  colnames(dat)[1:5] <- c("sample", "group", "cluster", "PC1", "PC2")
  
  # 去掉 OTHER 之外还是保留 global；targeted 里单独处理
  dat$group <- as.character(dat$group)
  dat$cluster <- as.character(dat$cluster)
  
  event_id <- sub("_sample_assignment\\.tsv$", "", basename(f))
  
  # 1. global Fisher: group × cluster
  tab_global <- table(dat$group, dat$cluster)
  global_p <- tryCatch(
    fisher.test(tab_global)$p.value,
    error = function(e) NA
  )
  
  # 2. 定义 cultivated-like cluster
  cultivated_like_cluster <- get_cultivated_like_cluster(dat)
  
  # 3. 二分类 Fisher
  res_wew_dem <- run_binary_fisher(dat, cultivated_like_cluster, "WEW", "DEM")
  res_wew_dur <- run_binary_fisher(dat, cultivated_like_cluster, "WEW", "DUR")
  
  # 合并 DEM + DUR
  dat2 <- dat
  dat2$group2 <- ifelse(dat2$group %in% c("DEM", "DUR"), "CULT", dat2$group)
  dat2 <- dat2[dat2$group2 %in% c("WEW", "CULT"), , drop = FALSE]
  dat2$group <- dat2$group2
  
  res_wew_cult <- run_binary_fisher(dat2, cultivated_like_cluster, "WEW", "CULT")
  
  # 4. 频率
  freq_tab <- table(dat$group, dat$cluster)
  freq_prop <- prop.table(freq_tab, 1)
  
  get_freq <- function(g, c) {
    if (g %in% rownames(freq_prop) && c %in% colnames(freq_prop)) {
      return(as.numeric(freq_prop[g, c]))
    } else {
      return(NA)
    }
  }
  
  res_list[[event_id]] <- data.frame(
    event_id = event_id,
    n_samples = nrow(dat),
    n_clusters = length(unique(dat$cluster)),
    cultivated_like_cluster = cultivated_like_cluster,
    global_fisher_p = global_p,
    
    freq_WEW = get_freq("WEW", cultivated_like_cluster),
    freq_DEM = get_freq("DEM", cultivated_like_cluster),
    freq_DUR = get_freq("DUR", cultivated_like_cluster),
    freq_OTHER = get_freq("OTHER", cultivated_like_cluster),
    
    p_WEW_vs_DEM = res_wew_dem$p.value,
    OR_WEW_vs_DEM = res_wew_dem$odds.ratio,
    
    p_WEW_vs_DUR = res_wew_dur$p.value,
    OR_WEW_vs_DUR = res_wew_dur$odds.ratio,
    
    p_WEW_vs_CULT = res_wew_cult$p.value,
    OR_WEW_vs_CULT = res_wew_cult$odds.ratio,
    
    stringsAsFactors = FALSE
  )
}

res_df <- do.call(rbind, res_list)

# =========================
# FDR 校正
# =========================
res_df$FDR_global <- p.adjust(res_df$global_fisher_p, method = "BH")
res_df$FDR_WEW_vs_DEM <- p.adjust(res_df$p_WEW_vs_DEM, method = "BH")
res_df$FDR_WEW_vs_DUR <- p.adjust(res_df$p_WEW_vs_DUR, method = "BH")
res_df$FDR_WEW_vs_CULT <- p.adjust(res_df$p_WEW_vs_CULT, method = "BH")

# 排序：先看 WEW vs CULT
res_df <- res_df[order(res_df$FDR_WEW_vs_CULT, res_df$p_WEW_vs_CULT), ]

write.table(
  res_df,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done. Output written to:\n", output_file, "\n")