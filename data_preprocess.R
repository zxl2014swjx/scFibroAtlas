################################## START ##################################
#' data_preprocess.R
#' @param expr_matrix 表达矩阵列表
#' @param method 标准化方法: "cpm", "logcpm", "tpm", "deseq2"
#' @param scale scale类型: "none", "gene", "sample"
#' @param log_transform 是否使用log转换
#' @param pseudo_count 标准化时当原始表达值为0时，使用log(raw +1)
#' @param gene_lengths 当输入数据是gene count时，标准化数据需要输入基因长度 gene length
#' @param cores 使用的总核心数
#' @param backend 并行后端: "foreach"
data_preprocess <- function(expr_matrix , method = "cpm",scale = "gene",log_transform = TRUE,pseudo_count = 1,gene_lengths = NULL,cores = NULL,backend = "foreach",outdir="0.data_preprocess") {
dir.create("0.Data_preprocess", showWarnings = FALSE)
expr_matrix <- read.table(expr_matrix,header=T,sep="\t",row.names=1,check.names=F)
  # 设置核心数
  if (is.null(cores)) {cores <- 1} 
  # 1. 标准化
  if (method == "cpm") {
    # CPM标准化
    lib_sizes <- colSums(expr_matrix)
    norm_matrix <- sweep(expr_matrix, 2, lib_sizes, FUN = "/") * 1e6  
  } else if (method == "logcpm") {
    # log2(CPM+1)
    lib_sizes <- colSums(expr_matrix)
    cpm_matrix <- sweep(expr_matrix, 2, lib_sizes, FUN = "/") * 1e6
    norm_matrix <- log2(cpm_matrix + pseudo_count) 
  } else if (method == "tpm") {
    # TPM标准化
    if (is.null(gene_lengths)) {
      stop("TPM normalization requires gene_lengths vector")
    }
    rpkm <- sweep(expr_matrix, 1, gene_lengths/1000, FUN = "/")
    rpkm <- sweep(rpkm, 2, colSums(expr_matrix)/1e6, FUN = "/")
    norm_matrix <- sweep(rpkm, 2, colSums(rpkm)/1e6, FUN = "/")  
  } else if (method == "deseq2") {
    # DESeq2标准化
    if (!require("DESeq2", quietly = TRUE)) {BiocManager::install("DESeq2")}
    library(DESeq2)
    col_data <- data.frame(sample = colnames(expr_matrix))
    rownames(col_data) <- col_data$sample
    dds <- DESeqDataSetFromMatrix(
      countData = expr_matrix,
      colData = col_data,
      design = ~ 1)
    dds <- estimateSizeFactors(dds)
    norm_matrix <- counts(dds, normalized = TRUE) 
  } else {
    stop(paste("Unsupported method:", method))
  }
  # 2. log转换（如果还没有）
  if (log_transform && !method %in% c("logcpm")) {
    norm_matrix <- log2(norm_matrix + pseudo_count)}
  # 3. scale
  if (scale != "none") {
    norm_matrix <- scale_matrix(
      norm_matrix, 
      by = scale, 
      cores = cores, 
      backend = backend)}
  # 保持行名和列名
  rownames(norm_matrix) <- rownames(expr_matrix)
  colnames(norm_matrix) <- colnames(expr_matrix)
  return(norm_matrix)
}

#' 通用的scale函数
scale_matrix <- function(expr_matrix, by = "gene", cores = NULL, backend = "none") {
  if (is.null(cores)) {cores <- 1}
  if (by == "gene") {
    # 按基因scale
    if (cores <= 1 || backend == "none" || nrow(expr_matrix) < 1000) {
      return(t(scale(t(expr_matrix))))
    }
    if (backend == "foreach") {
      if (!require("foreach", quietly = TRUE) || !require("doParallel", quietly = TRUE)) {
        install.packages(c("foreach", "doParallel"))
      }
      library(foreach)
      library(doParallel)
      cl <- makeCluster(cores)
      registerDoParallel(cl)
      # 分割基因
      n_genes <- nrow(expr_matrix)
      block_size <- ceiling(n_genes / cores)
      gene_indices <- split(1:n_genes, ceiling(seq_along(1:n_genes) / block_size))
      scaled_blocks <- foreach(i = 1:length(gene_indices), .combine = rbind) %dopar% {
        idx <- gene_indices[[i]]
        chunk <- expr_matrix[idx, , drop = FALSE]
        t(scale(t(chunk)))
      }
      stopCluster(cl)
      rownames(scaled_blocks) <- rownames(expr_matrix)
      colnames(scaled_blocks) <- colnames(expr_matrix)
      return(scaled_blocks)
    } else {
      return(t(scale(t(expr_matrix))))
    }
  } else if (by == "sample") {
    # 按样本scale
    return(scale(expr_matrix))
  } else {
    stop(paste("Unsupported scale type:", by))
  }
}
################################## END ##################################
