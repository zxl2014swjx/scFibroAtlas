################################## START ##################################
#' analysis_subtypes.R
#' 使用多种解卷积方法计算成纤维细胞亚型比例，并提供完整的可视化
#' @param preprocessed 输入的表达矩阵文件路径或矩阵对象
#' @param method 解卷积方法: "cibersort", "GSVA", "ssGSEA", "xCell", "MCPcounter", "ESTIMATE"
#' @param reference 参考数据文件路径或矩阵对象
#' @param core 使用的核心数
#' @param output_dir 输出目录
#' @param create_plots 是否创建可视化图形
#' @param plot_format 图形格式: "pdf", "png", "both"
#' @param verbose 是否显示详细信息
#' @return 包含亚型分数、模型信息和可视化结果的对象
analysis_subtypes <- function(preprocessed, method = "cibersort",reference = NULL,core = 5,output_dir = "1.Subtypes_deconvolution",create_plots = TRUE,plot_format = "both",verbose = TRUE) {
  start_time <- Sys.time()
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  if (is.null(core)) {core <- parallel::detectCores() - 1}
  if (verbose) message(paste("Using", core, "cores for deconvolution"))
  if (verbose) message("Step 1: Loading and preprocessing input data...")
  expr_data <- load_expression_data(preprocessed)
  if (verbose) message("Step 2: Loading fibroblast reference data...")
  fibro_ref <- load_reference_data(reference, method = method)
  if (verbose) message(paste("Step 3: Running deconvolution with", method, "..."))
  deconv_results <- run_deconvolution(expr_data = expr_data,reference = fibro_ref,method = method,cores = core)
  if (verbose) message("Step 4: Post-processing results...")
  processed_results <- process_deconvolution_results(deconv_results = deconv_results,method = method)
  if (verbose) message("Step 5: Performing quality control...")
  qc_results <- perform_quality_control(expr_data = expr_data,deconv_results = processed_results,method = method)
  if (create_plots) {
    if (verbose) message("Step 6: Generating visualizations...")
    plots <- generate_deconvolution_visualizations(deconv_results = processed_results,expr_data = expr_data,reference = fibro_ref,method = method,output_dir = output_dir,plot_format = plot_format)
  } else {
    plots <- NULL
  }
  if (verbose) message("Step 7: Saving results...")
  save_results(deconv_results = processed_results,qc_results = qc_results,plots = plots,output_dir = output_dir)
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  if (verbose) message(paste("Deconvolution analysis completed in", round(runtime, 2), "minutes!"))
  result <- list(deconv_scores = processed_results$scores,deconv_stats = processed_results$stats,quality_control = qc_results,plots = plots,parameters = list(method = method,n_samples = ncol(expr_data),n_celltypes = ncol(processed_results$scores),runtime = runtime,timestamp = Sys.time()),
    input_data_summary = list(n_genes = nrow(expr_data),n_samples = ncol(expr_data),data_range = range(expr_data, na.rm = TRUE),data_mean = mean(expr_data, na.rm = TRUE)),
    output_dir = output_dir)
  class(result) <- "FibroblastDeconvolution"
  return(result)
}
load_expression_data <- function(data_input) {
  if (is.character(data_input)) {
    if (!file.exists(data_input)) {stop(paste("Input file not found:", data_input))}
    file_ext <- tools::file_ext(data_input)
    if (file_ext %in% c("txt", "tsv", "tab")) {
      data <- read.table(data_input, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    } else if (file_ext == "csv") {
      data <- read.csv(data_input, header = TRUE, row.names = 1, check.names = FALSE)
    } else if (file_ext == "rds") {
      data <- readRDS(data_input)
    } else {
      stop(paste("Unsupported file format:", file_ext))
    }
  } else if (is.matrix(data_input) || inherits(data_input, "Matrix")) {
    data <- data_input
  } else {
    stop("data_input must be a file path or a matrix")
  }
  if (!is.matrix(data)) {data <- as.matrix(data)}
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop("Input data is empty")
  }
  zero_rows <- rowSums(data, na.rm = TRUE) == 0
  if (sum(zero_rows) > 0) {
    warning(paste("Removing", sum(zero_rows), "genes with all zero expression"))
    data <- data[!zero_rows, ]
  }
  if (max(data, na.rm = TRUE) > 100) {
    message("Data appears to be raw counts, applying log2(CPM+1) transformation...")
    lib_sizes <- colSums(data, na.rm = TRUE)
    cpm <- sweep(data, 2, lib_sizes / 1e6, "/")
    data <- log2(cpm + 1)
  }
  return(data)
}
load_reference_data <- function(reference_input, method = "cibersort") {
  if (is.null(reference_input)) {
    message("Using built-in fibroblast reference data...")
    if (method == "cibersort") {
      data(inbuilt_fibroblast_signature)
      return(inbuilt_fibroblast_signature)
    } else if (method == "GSVA" || method == "ssGSEA") {
      data(inbuilt_fibroblast_genesets)
      return(inbuilt_fibroblast_genesets)
    } else {
      stop(paste("No built-in reference for method:", method))
    }
  }
  if (is.character(reference_input)) {
    if (!file.exists(reference_input)) {
      stop(paste("Reference file not found:", reference_input))
    }
    file_ext <- tools::file_ext(reference_input)
    if (file_ext %in% c("txt", "tsv", "tab")) {
      data <- read.table(reference_input, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
    } else if (file_ext == "csv") {
      data <- read.csv(reference_input, header = TRUE, row.names = 1, check.names = FALSE)
    } else if (file_ext == "rds") {
      data <- readRDS(reference_input)
    } else if (file_ext == "gmt") {
      data <- read_gmt_file(reference_input)
    } else {
      stop(paste("Unsupported reference file format:", file_ext))
    }
  } else if (is.list(reference_input) || is.matrix(reference_input)) {
    data <- reference_input
  } else {
    stop("reference must be a file path, matrix, or list")
  }
  return(data)
}
#' 读取GMT文件
read_gmt_file <- function(gmt_file) {
  con <- file(gmt_file, "r")
  genesets <- list()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) >= 3) {
      geneset_name <- fields[1]
      geneset_desc <- fields[2]
      genes <- fields[3:length(fields)]
      genesets[[geneset_name]] <- genes
    }
  }
  close(con)
  return(genesets)
}
#' 运行解卷积算法
run_deconvolution <- function(expr_data, reference, method = "cibersort", cores = 1) {
  method_functions <- list(cibersort = run_cibersort,GSVA = run_gsva,ssGSEA = run_ssgsea,xCell = run_xcell,MCPcounter = run_mcpcounter,ESTIMATE = run_estimate)
  if (!method %in% names(method_functions)) {
    stop(paste("Unsupported method:", method, "\nSupported methods:", paste(names(method_functions), collapse = ", ")))
  }
  result <- method_functions[[method]](expr_data, reference, cores)
  return(result)
}
#' 运行CIBERSORT
run_cibersort <- function(expr_data, reference, cores = 1) {
  if (!require("CIBERSORT", quietly = TRUE)) {
    message("Installing CIBERSORT from GitHub...")
    devtools::install_github("Moonerss/CIBERSORT", quiet = TRUE)}
  library(CIBERSORT)
  if (!is.matrix(reference)) {
    reference <- as.matrix(reference)
  }
  common_genes <- intersect(rownames(expr_data), rownames(reference))
  if (length(common_genes) < 100) {
    warning(paste("Only", length(common_genes), "common genes between expression data and signature"))
  }
  expr_subset <- expr_data[common_genes, , drop = FALSE]
  ref_subset <- reference[common_genes, , drop = FALSE]
  message("Running CIBERSORT (this may take several minutes)...")
  results <- CIBERSORT::cibersort(sig_matrix = ref_subset,mixture_file = expr_subset,perm = 100,QN = FALSE)
  return(results)
}
#' 运行GSVA
run_gsva <- function(expr_data, genesets, cores = 1) {
  if (!require("GSVA", quietly = TRUE)) {BiocManager::install("GSVA", ask = FALSE)}
  library(GSVA)
  if (cores > 1) {
    if (!require("BiocParallel", quietly = TRUE)) {BiocManager::install("BiocParallel", ask = FALSE)}
    library(BiocParallel)
    param <- MulticoreParam(workers = cores)
    bpparam <- param
  } else {
    bpparam <- NULL
  }
  results <- gsva(expr = as.matrix(expr_data),gset.idx.list = genesets,method = "gsva",kcdf = "Gaussian",verbose = FALSE,BPPARAM = bpparam)
  return(results)
}
#' 运行ssGSEA
run_ssgsea <- function(expr_data, genesets, cores = 1) {
  if (!require("GSVA", quietly = TRUE)) {BiocManager::install("GSVA", ask = FALSE)}
  library(GSVA)
  if (cores > 1) {
    if (!require("BiocParallel", quietly = TRUE)) {BiocManager::install("BiocParallel", ask = FALSE)}
    library(BiocParallel)
    param <- MulticoreParam(workers = cores)
    bpparam <- param
  } else {
    bpparam <- NULL
  }
  # 运行ssGSEA
  results <- gsva(expr = as.matrix(expr_data),gset.idx.list = genesets,method = "ssgsea",verbose = FALSE,BPPARAM = bpparam)
  return(results)
}
#' 运行xCell
run_xcell <- function(expr_data, reference, cores = 1) {
  if (!require("xCell", quietly = TRUE)) {
    message("Installing xCell...")
    devtools::install_github("dviraran/xCell", quiet = TRUE)}
  library(xCell)
  results <- xCellAnalysis(expr_data)
  return(results)
}
#' 运行MCPcounter
run_mcpcounter <- function(expr_data, reference, cores = 1) {
  if (!require("MCPcounter", quietly = TRUE)) {
    install_github("ebecht/MCPcounter",ref="master", subdir="Source")}
  library(MCPcounter)
  results <- MCPcounter.estimate(expr_data, featuresType = "HUGO_symbols")
  return(results)
}
#' 运行ESTIMATE
run_estimate <- function(expr_data, reference, cores = 1) {
  if (!require("estimate", quietly = TRUE)) {
    install.packages("estimate", repos="http://r-forge.r-project.org", dependencies=TRUE)}
  library(estimate)
  temp_input <- tempfile(fileext = ".txt")
  write.table(expr_data, temp_input, sep = "\t", quote = FALSE, col.names = NA)
  filterCommonGenes(temp_input, temp_input, id = "GeneSymbol")
  estimateScore(temp_input, temp_input, platform = "affymetrix")
  results <- read.table(paste0(temp_input, "_scores.gct"), header = TRUE, sep = "\t", skip = 2, row.names = 1)
  unlink(temp_input)
  unlink(paste0(temp_input, "_scores.gct"))
  return(t(results))
}
#' 处理解卷积结果
process_deconvolution_results <- function(deconv_results, method = "cibersort") {
  scores <- NULL
  stats <- NULL
  if (method == "cibersort") {
    if (is.data.frame(deconv_results)) {
      stat_cols <- c("P-value", "Correlation", "RMSE")
      score_cols <- setdiff(colnames(deconv_results), stat_cols)
      scores <- deconv_results[, score_cols, drop = FALSE]
      if (all(stat_cols %in% colnames(deconv_results))) {
        stats <- deconv_results[, stat_cols, drop = FALSE]
      }
    } else {
      scores <- deconv_results
    }    
  } else if (method %in% c("GSVA", "ssGSEA")) {
    if (is.matrix(deconv_results)) {
      scores <- t(deconv_results)
    } else {
      scores <- deconv_results
    }
  } else if (method == "xCell") {
    if (is.matrix(deconv_results)) {
      scores <- t(deconv_results)
    } else {
      scores <- deconv_results
    }
  } else {
    scores <- deconv_results
  }
  if (!is.data.frame(scores) && !is.null(scores)) {
    scores <- as.data.frame(scores)
  }
  return(list(scores = scores, stats = stats))
}
#' 执行质量控制
perform_quality_control <- function(expr_data, deconv_results, method = "cibersort") {
  scores <- deconv_results$scores
  qc_metrics <- list()
  qc_metrics$score_range <- list(min = min(scores, na.rm = TRUE),max = max(scores, na.rm = TRUE),mean = mean(as.matrix(scores), na.rm = TRUE),sd = sd(as.matrix(scores), na.rm = TRUE))
  qc_metrics$missing_values <- list(total_na = sum(is.na(scores)),percent_na = mean(is.na(scores)) * 100)
  if (method == "cibersort") {
    negative_scores <- scores < 0
    qc_metrics$negative_scores <- list(
      total_negative = sum(negative_scores, na.rm = TRUE),
      percent_negative = mean(negative_scores, na.rm = TRUE) * 100
    )
  }
  if (method == "cibersort") {
    row_sums <- rowSums(scores, na.rm = TRUE)
    qc_metrics$row_sums <- list(
      min = min(row_sums),
      max = max(row_sums),
      mean = mean(row_sums),
      sd = sd(row_sums)
    )
  }
  if (ncol(scores) > 1) {
    cor_matrix <- cor(scores, use = "complete.obs")
    qc_metrics$correlation <- cor_matrix
  }
  return(qc_metrics)
}
#' 生成解卷积可视化
generate_deconvolution_visualizations <- function(deconv_results, expr_data, reference, method = "cibersort",output_dir = "deconvolution_results",plot_format = "both") {
  scores <- deconv_results$scores
  plots <- list()
  plots$heatmap <- plot_deconvolution_heatmap(scores = scores,output_dir = output_dir,plot_format = plot_format)
  plots$barplot <- plot_deconvolution_barplot(scores = scores,output_dir = output_dir,plot_format = plot_format)
  if (ncol(scores) > 2) {
    plots$correlation_heatmap <- plot_correlation_heatmap(
      scores = scores,
      output_dir = output_dir,
      plot_format = plot_format
    )
  }
  plots$pca <- plot_deconvolution_pca(scores = scores,output_dir = output_dir,plot_format = plot_format)
  plots$boxplot <- plot_deconvolution_boxplot(scores = scores,output_dir = output_dir,plot_format = plot_format)
  plots$violin <- plot_deconvolution_violin(scores = scores,output_dir = output_dir,plot_format = plot_format)
  plots$stacked_area <- plot_stacked_area(scores = scores,output_dir = output_dir,plot_format = plot_format)
  if (nrow(scores) <= 20) {  # 雷达图不适合太多样本
    plots$radar <- plot_radar_chart(
      scores = scores,
      output_dir = output_dir,
      plot_format = plot_format
    )
  }
  if (!is.null(deconv_results$stats)) {
    plots$qc_plots <- plot_quality_control(
      stats = deconv_results$stats,
      output_dir = output_dir,
      plot_format = plot_format
    )
  }
  return(plots)
}
plot_deconvolution_heatmap <- function(scores, output_dir, plot_format) {
  if (!require("pheatmap", quietly = TRUE) || !require("RColorBrewer", quietly = TRUE)) {
    install.packages(c("pheatmap", "RColorBrewer"))}
  library(pheatmap)
  library(RColorBrewer)
  plot_data <- as.matrix(scores)
  color_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
  p <- pheatmap(t(plot_data),main = "Fibroblast Subtype Scores Heatmap",color = color_palette,scale = "row",show_rownames = TRUE,show_colnames = ncol(plot_data) <= 50,cluster_rows = TRUE,cluster_cols = TRUE,fontsize_row = 8,fontsize_col = 8,angle_col = 45)
  save_plot(p, filename = file.path(output_dir, "deconvolution_heatmap"),format = plot_format,width = 12, height = 8)
  return(p)
}
#' 绘制解卷积条形图
plot_deconvolution_barplot <- function(scores, output_dir, plot_format) {
  if (!require("ggplot2", quietly = TRUE) || !require("reshape2", quietly = TRUE)) {
    install.packages(c("ggplot2", "reshape2"))}
  library(ggplot2)
  library(reshape2)
  plot_data <- scores
  plot_data$Sample <- rownames(plot_data)
  plot_data_melted <- melt(plot_data, id.vars = "Sample", variable.name = "Subtype", value.name = "Score")
  p <- ggplot(plot_data_melted, aes(x = Sample, y = Score, fill = Subtype)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),legend.position = "right",plot.title = element_text(hjust = 0.5)) +
  labs(title = "Fibroblast Subtype Composition",x = "Sample", y = "Proportion") +
  scale_fill_brewer(palette = "Set3")
  save_plot(p, filename = file.path(output_dir, "deconvolution_barplot"),format = plot_format, width = 14, height = 6)
  return(p)
}
#' 绘制相关性热图
plot_correlation_heatmap <- function(scores, output_dir, plot_format) {
  if (!require("corrplot", quietly = TRUE)) {
    install.packages("corrplot")
  }
  library(corrplot)
  cor_matrix <- cor(scores, use = "complete.obs")
  if (plot_format %in% c("pdf", "both")) {
    pdf(file.path(output_dir, "correlation_heatmap.pdf"), width = 10, height = 8)
    corrplot(cor_matrix, method = "color",type = "upper",order = "hclust",tl.col = "black",tl.srt = 45,addCoef.col = "black",number.cex = 0.7,main = "Subtype Correlation Matrix",mar = c(0, 0, 2, 0))
    dev.off()
  }
  if (plot_format %in% c("png", "both")) {
    png(file.path(output_dir, "correlation_heatmap.png"), width = 1000, height = 800, res = 150)
    corrplot(cor_matrix, method = "color",type = "upper",order = "hclust",tl.col = "black",tl.srt = 45,addCoef.col = "black",number.cex = 0.7,main = "Subtype Correlation Matrix",mar = c(0, 0, 2, 0))
    dev.off()
  }
  return(cor_matrix)
}
#' 绘制解卷积PCA图
plot_deconvolution_pca <- function(scores, output_dir, plot_format) {
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)
  pca_result <- prcomp(scores, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$Sample <- rownames(scores)
  variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(hjust = 0, vjust = 0, size = 3, check_overlap = TRUE) +
    theme_bw() +
    labs(title = "PCA of Fibroblast Subtype Scores",x = paste0("PC1 (", variance_explained[1], "%)"),y = paste0("PC2 (", variance_explained[2], "%)"))
  save_plot(p, filename = file.path(output_dir, "deconvolution_pca"),format = plot_format,width = 10, height = 8)
  return(list(plot = p, pca_result = pca_result))
}
#' 绘制解卷积箱线图
plot_deconvolution_boxplot <- function(scores, output_dir, plot_format) {
  if (!require("ggplot2", quietly = TRUE) || !require("reshape2", quietly = TRUE)) {
    install.packages(c("ggplot2", "reshape2"))}
  library(ggplot2)
  library(reshape2)
  plot_data <- melt(as.matrix(scores))
  colnames(plot_data) <- c("Sample", "Subtype", "Score")
  p <- ggplot(plot_data, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_boxplot(outlier.shape = 19, outlier.size = 1) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",plot.title = element_text(hjust = 0.5)) +
    labs(title = "Distribution of Subtype Scores",x = "Fibroblast Subtype",y = "Score") +
    scale_fill_brewer(palette = "Set3")
    save_plot(p, filename = file.path(output_dir, "deconvolution_boxplot"),format = plot_format,width = 10, height = 6)
  return(p)
}
#' 绘制解卷积小提琴图
plot_deconvolution_violin <- function(scores, output_dir, plot_format) {
  if (!require("ggplot2", quietly = TRUE) || !require("reshape2", quietly = TRUE)) {
    install.packages(c("ggplot2", "reshape2"))}
  library(ggplot2)
  library(reshape2)
  plot_data <- melt(as.matrix(scores))
  colnames(plot_data) <- c("Sample", "Subtype", "Score")
  p <- ggplot(plot_data, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none",plot.title = element_text(hjust = 0.5)) +
    labs(title = "Violin Plot of Subtype Scores",x = "Fibroblast Subtype",y = "Score") +
    scale_fill_brewer(palette = "Set3")
  save_plot(p, filename = file.path(output_dir, "deconvolution_violin"),format = plot_format,width = 10, height = 6) 
  return(p)
}
#' 绘制堆积面积图
plot_stacked_area <- function(scores, output_dir, plot_format) {
  if (!require("ggplot2", quietly = TRUE) || !require("reshape2", quietly = TRUE)) {
    install.packages(c("ggplot2", "reshape2"))}
  library(ggplot2)
  library(reshape2)
  plot_data <- scores
  plot_data$Sample <- factor(rownames(plot_data), levels = rownames(plot_data))
  plot_data_melted <- melt(plot_data, id.vars = "Sample", variable.name = "Subtype", value.name = "Score")
  plot_data_melted <- plot_data_melted[order(plot_data_melted$Sample), ]

  p <- ggplot(plot_data_melted, aes(x = Sample, y = Score, fill = Subtype)) +
    geom_area(position = "stack", alpha = 0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),legend.position = "right",plot.title = element_text(hjust = 0.5)) +
    labs(title = "Stacked Area Plot of Subtype Composition",x = "Sample", y = "Proportion") +
    scale_fill_brewer(palette = "Set3")
  save_plot(p, filename = file.path(output_dir, "stacked_area_plot"),format = plot_format,width = 14, height = 6)
  return(p)
}
#' 绘制雷达图
plot_radar_chart <- function(scores, output_dir, plot_format) {
  if (!require("fmsb", quietly = TRUE) || !require("scales", quietly = TRUE)) {
    install.packages(c("fmsb", "scales"))}
  library(fmsb)
  library(scales)
  plot_data <- as.data.frame(scores)
  if (nrow(plot_data) > 10) {
    message("Too many samples for radar chart, plotting first 10 only")
    plot_data <- plot_data[1:10, ]
  }
  plot_data_scaled <- as.data.frame(lapply(plot_data, rescale))
  max_min <- data.frame(max = rep(1, ncol(plot_data_scaled)),min = rep(0, ncol(plot_data_scaled)))
  colnames(max_min) <- colnames(plot_data_scaled)
  rownames(max_min) <- c("max", "min")
  plot_data_radar <- rbind(max_min, plot_data_scaled)
  colors_border <- scales::hue_pal()(nrow(plot_data_scaled))
  colors_in <- scales::alpha(colors_border, 0.3)
  if (plot_format %in% c("pdf", "both")) {
    pdf(file.path(output_dir, "radar_chart.pdf"), width = 10, height = 8)
    op <- par(mar = c(1, 2, 2, 1))
    radarchart(plot_data_radar, axistype = 1,pcol = colors_border,pfcol = colors_in,plwd = 2,plty = 1,cglcol = "grey",cglty = 1,axislabcol = "grey",caxislabels = seq(0, 1, 0.25),cglwd = 0.8,vlcex = 0.8,title = "Radar Chart of Subtype Scores")
    legend(x = 1.3, y = 1.3, legend = rownames(plot_data_scaled),bty = "n", pch = 20, col = colors_border,text.col = "black", cex = 0.8, pt.cex = 1.5)
    par(op)
    dev.off()
  }
  if (plot_format %in% c("png", "both")) {
    png(file.path(output_dir, "radar_chart.png"), 
        width = 1000, height = 800, res = 150)
    op <- par(mar = c(1, 2, 2, 1))
    radarchart(plot_data_radar, axistype = 1,pcol = colors_border,pfcol = colors_in,plwd = 2,plty = 1,cglcol = "grey",cglty = 1,axislabcol = "grey",caxislabels = seq(0, 1, 0.25),cglwd = 0.8,vlcex = 0.8,title = "Radar Chart of Subtype Scores")
    legend(x = 1.3, y = 1.3, legend = rownames(plot_data_scaled),bty = "n", pch = 20, col = colors_border,text.col = "black", cex = 0.8, pt.cex = 1.5)
    par(op)
    dev.off()
  }
  return(plot_data_radar)
}
#' 绘制质量控制图
plot_quality_control <- function(stats, output_dir, plot_format) {
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")}
  library(ggplot2)
  plots <- list()
  if ("P-value" %in% colnames(stats)) {
    p1 <- ggplot(stats, aes(x = `P-value`)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "white") +
      geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
      theme_bw() +
      labs(title = "Distribution of P-values",x = "P-value", y = "Count")
    plots$pvalue_distribution <- p1
    save_plot(p1, filename = file.path(output_dir, "pvalue_distribution"),format = plot_format,width = 8, height = 6)
  }
  # 2. 相关性分布
  if ("Correlation" %in% colnames(stats)) {
    p2 <- ggplot(stats, aes(x = Correlation)) +
      geom_histogram(bins = 30, fill = "darkgreen", color = "white") +
      theme_bw() +
      labs(title = "Distribution of Correlations",x = "Correlation", y = "Count")
    plots$correlation_distribution <- p2
    save_plot(p2, filename = file.path(output_dir, "correlation_distribution"),format = plot_format,width = 8, height = 6)
  }
  # 3. RMSE分布
  if ("RMSE" %in% colnames(stats)) {
    p3 <- ggplot(stats, aes(x = RMSE)) +
      geom_histogram(bins = 30, fill = "darkorange", color = "white") +
      theme_bw() +
      labs(title = "Distribution of RMSE",x = "RMSE", y = "Count")
    plots$rmse_distribution <- p3
    save_plot(p3, filename = file.path(output_dir, "rmse_distribution"),format = plot_format,width = 8, height = 6)
  }
  return(plots)
}
save_plot <- function(plot_obj, filename, format = "both", width = 10, height = 8) {
  if (format %in% c("pdf", "both")) {
    pdf(paste0(filename, ".pdf"), width = width, height = height)
    print(plot_obj)
    dev.off()
  }
  if (format %in% c("png", "both")) {
    png(paste0(filename, ".png"), width = width * 100, height = height * 100, res = 150)
    print(plot_obj)
    dev.off()
  }
}
save_results <- function(deconv_results, qc_results, plots, output_dir) {
  write.table(deconv_results$scores,file.path(output_dir, "deconvolution_scores.txt"),sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  if (!is.null(deconv_results$stats)) {
    write.table(deconv_results$stats,file.path(output_dir, "deconvolution_stats.txt"),sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  }
  saveRDS(qc_results, file.path(output_dir, "quality_control.rds"))
  summary_df <- data.frame(
    Metric = c("Number of samples", "Number of subtypes", "Mean score", "SD of scores"),
    Value = c(nrow(deconv_results$scores), ncol(deconv_results$scores),mean(as.matrix(deconv_results$scores), na.rm = TRUE),sd(as.matrix(deconv_results$scores), na.rm = TRUE))
  )
  write.csv(summary_df, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)
  message(paste("Results saved to:", output_dir))
}
################################## END ##################################
