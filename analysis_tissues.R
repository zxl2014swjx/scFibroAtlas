################################## START ##################################
# analysis_tissues.R
# 输入: data_preprocess.csv, 用户自定义sample_metadata.csv
# sample_metadata的数据格式为
# 创建样本元数据
# sample_metadata <- data.frame(
#  tissue = sample(paste0("Tissue_", 1:n_tissues), n_samples, replace = TRUE),
#  age = sample(20:80, n_samples, replace = TRUE),
#  sex = sample(c("M", "F"), n_samples, replace = TRUE))
# rownames(sample_metadata) <- colnames(score_matrix)
# 输出: 2.Tissue_enrichment文件夹及其所有分析结果
# 设置工作目录和加载包
cat("设置工作环境...\n")
setwd(".")  # 设置工作目录
# 检查并安装必要的包
required_packages <- c("reshape2", "pheatmap", "ggplot2","RColorBrewer","patchwork","corrplot")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("安装缺失的包:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(reshape2)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
  library(patchwork)
  library(corrplot)
})
cat("\n创建输出目录...\n")
dir.create("2.Tissue_enrichment", showWarnings = FALSE)
cat("\n1. 正在加载数据...\n")
score_matrix<-read.table("0.Data_preprocess/data_preprocessed.txt",header=T,sep="\t",row.names=1,check.names=F)

tissue_enrichment_analysis_fixed <- function(score_matrix,sample_metadata,n_components = 10,clustering_method = "kmeans",n_clusters = NULL) {
  validate_inputs_fixed(score_matrix, sample_metadata)
  result <- list()
  result$parameters <- list(n_components = n_components,clustering_method = clustering_method,n_clusters = n_clusters)
  message("Step 1: Data preprocessing...")
  processed_data <- preprocess_scores_fixed(score_matrix = score_matrix,sample_metadata = sample_metadata)
  result$processed_data <- processed_data
  message("Step 2: Tissue-specific enrichment analysis...")
  tissue_specificity <- analyze_tissue_specificity_fixed(score_matrix = processed_data$scores,metadata = processed_data$metadata)
  result$tissue_specificity <- tissue_specificity
  message("Step 3: Dimensionality reduction...")
  dim_reduction_results <- perform_dimensionality_reduction_fixed(score_matrix = processed_data$scores,n_components = n_components)
  result$dim_reduction <- dim_reduction_results
  if (is.null(dim_reduction_results$embedding) || nrow(dim_reduction_results$embedding) == 0) {
    warning("Dimensionality reduction failed, using original data for clustering")
    dim_reduction_results$embedding <- t(processed_data$scores)[, 1:min(10, ncol(processed_data$scores))]}
  message("Step 4: Clustering analysis...")
  clustering_results <- perform_clustering_fixed(reduced_data = dim_reduction_results$embedding, method = clustering_method,n_clusters = n_clusters)
  result$clustering <- clustering_results
  message("Step 5: Enrichment analysis...")
  enrichment_results <- perform_enrichment_analysis_fixed(score_matrix = processed_data$scores,clusters = clustering_results$clusters,metadata = processed_data$metadata)
  result$enrichment <- enrichment_results
  message("Step 6: Generating visualizations...")
  plots <- generate_visualizations_fixed(result)
  result$plots <- plots
  message("Analysis completed successfully!")
  return(result)
}
#' 组织特异性分析
analyze_tissue_specificity_fixed <- function(score_matrix, metadata) {
  tissues <- unique(metadata$tissue)
  celltypes <- rownames(score_matrix)
  tissue_mean_scores <- matrix(0, nrow = length(celltypes), ncol = length(tissues))
  rownames(tissue_mean_scores) <- celltypes
  colnames(tissue_mean_scores) <- tissues
  for (tissue in tissues) {
    tissue_samples <- rownames(metadata)[metadata$tissue == tissue]
    if (length(tissue_samples) > 0) {
      tissue_scores <- score_matrix[, tissue_samples, drop = FALSE]
      tissue_mean_scores[, tissue] <- rowMeans(tissue_scores, na.rm = TRUE)
    }
  }
  # 2. 计算组织特异性指数 (Tau)
  tau_scores <- calculate_tau_index(tissue_mean_scores)
  # 3. 识别组织特异性细胞类型
  tissue_specific_celltypes <- identify_tissue_specific_celltypes(
    tissue_mean_scores = tissue_mean_scores,
    tau_scores = tau_scores
  )
  # 4. 组织间相关性
  tissue_correlation <- cor(tissue_mean_scores, use = "complete.obs")
  return(list(
    tissue_mean_scores = tissue_mean_scores,
    tau_scores = tau_scores,
    tissue_specific_celltypes = tissue_specific_celltypes,
    tissue_correlation = tissue_correlation
  ))
}
#' 计算Tau组织特异性指数
calculate_tau_index <- function(tissue_mean_scores) {
  # Tau = (1 - x_i/max(x)) / (n - 1) 的和
  n_tissues <- ncol(tissue_mean_scores)
  tau <- numeric(nrow(tissue_mean_scores))
  names(tau) <- rownames(tissue_mean_scores)
  for (i in 1:nrow(tissue_mean_scores)) {
    x <- tissue_mean_scores[i, ]
    max_val <- max(x, na.rm = TRUE)
    if (max_val > 0) {
      tau[i] <- sum(1 - x/max_val, na.rm = TRUE) / (n_tissues - 1)
    } else {
      tau[i] <- 0
    }
  }
  return(tau)
}
#' 识别组织特异性细胞类型
identify_tissue_specific_celltypes <- function(tissue_mean_scores, tau_scores, threshold = 0.8) {
  celltypes <- rownames(tissue_mean_scores)
  tissues <- colnames(tissue_mean_scores)
  # 高Tau值表示高组织特异性
  specific_celltypes <- celltypes[tau_scores > threshold]
  # 对于每个高Tau的细胞类型，确定其主要组织
  tissue_assignments <- list()
  for (celltype in specific_celltypes) {
    scores <- tissue_mean_scores[celltype, ]
    max_tissue <- names(which.max(scores))
    tissue_assignments[[celltype]] <- max_tissue
  }
  return(list(
    specific_celltypes = specific_celltypes,
    tissue_assignments = tissue_assignments,
    tau_threshold = threshold
  ))
}
#' 富集分析
perform_enrichment_analysis_fixed <- function(score_matrix, clusters, metadata) {
  tissues <- unique(metadata$tissue)
  cluster_ids <- unique(clusters)
  # 1. 聚类与组织的关系
  cluster_tissue_enrichment <- matrix(0, nrow = length(cluster_ids), ncol = length(tissues))
  rownames(cluster_tissue_enrichment) <- paste0("Cluster", cluster_ids)
  colnames(cluster_tissue_enrichment) <- tissues
  for (cluster_id in cluster_ids) {
    cluster_samples <- names(clusters)[clusters == cluster_id]
    for (tissue in tissues) {
      tissue_samples <- rownames(metadata)[metadata$tissue == tissue]
      overlap <- intersect(cluster_samples, tissue_samples)
      # 超几何检验
      p_value <- phyper(
        q = length(overlap) - 1,
        m = length(tissue_samples),
        n = nrow(metadata) - length(tissue_samples),
        k = length(cluster_samples),
        lower.tail = FALSE
      )
      cluster_tissue_enrichment[paste0("Cluster", cluster_id), tissue] <- -log10(p_value)
    }
  }
  # 2. 细胞类型在聚类中的富集
  celltypes <- rownames(score_matrix)
  cluster_celltype_enrichment <- matrix(0, nrow = length(celltypes), ncol = length(cluster_ids))
  rownames(cluster_celltype_enrichment) <- celltypes
  colnames(cluster_celltype_enrichment) <- paste0("Cluster", cluster_ids)
  for (cluster_id in cluster_ids) {
    cluster_samples <- names(clusters)[clusters == cluster_id]
    for (celltype in celltypes) {
      # 比较该细胞类型在聚类内外的得分
      in_cluster <- score_matrix[celltype, cluster_samples]
      out_cluster <- score_matrix[celltype, setdiff(colnames(score_matrix), cluster_samples)]
      
      if (length(in_cluster) > 1 && length(out_cluster) > 1) {
        test_result <- wilcox.test(in_cluster, out_cluster)
        cluster_celltype_enrichment[celltype, paste0("Cluster", cluster_id)] <- 
          -log10(test_result$p.value) * sign(mean(in_cluster) - mean(out_cluster))
      }
    }
  }
  return(list(
    cluster_tissue_enrichment = cluster_tissue_enrichment,
    cluster_celltype_enrichment = cluster_celltype_enrichment
  ))
}
#' 修复的验证输入函数
validate_inputs_fixed <- function(score_matrix, sample_metadata) {
  if (!is.matrix(score_matrix) && !inherits(score_matrix, "Matrix")) {
    score_matrix <- as.matrix(score_matrix)}
  if (is.null(rownames(score_matrix))) {
    stop("score_matrix must have rownames (cell types)")}
  if (is.null(colnames(score_matrix))) {
    stop("score_matrix must have colnames (sample names)")}
  if (!"tissue" %in% colnames(sample_metadata)) {
    stop("sample_metadata must contain 'tissue' column")}
  if (!is.data.frame(sample_metadata)) {
    sample_metadata <- as.data.frame(sample_metadata)}
  if (is.null(rownames(sample_metadata))) {
    rownames(sample_metadata) <- paste0("Sample", 1:nrow(sample_metadata))}
  common_samples <- intersect(colnames(score_matrix), rownames(sample_metadata))
  if (length(common_samples) == 0) {
    stop("No common samples between score_matrix and sample_metadata")}
  message(paste("Analysis will be performed on", length(common_samples), "samples"))
  return(TRUE)
}
#' 修复的数据预处理函数
preprocess_scores_fixed <- function(score_matrix, sample_metadata) {
  common_samples <- intersect(colnames(score_matrix), rownames(sample_metadata))
  score_matrix <- score_matrix[, common_samples, drop = FALSE]
  sample_metadata <- sample_metadata[common_samples, , drop = FALSE]
  if (!is.data.frame(sample_metadata)) {
    sample_metadata <- as.data.frame(sample_metadata)
  }
  sample_quality <- colSums(score_matrix, na.rm = TRUE)
  if (all(sample_quality == 0)) {
    warning("All samples have zero scores")
    keep_samples <- rep(TRUE, ncol(score_matrix))
  } else {
    keep_samples <- sample_quality > quantile(sample_quality, 0.05, na.rm = TRUE)
  }
  score_matrix <- score_matrix[, keep_samples, drop = FALSE]
  sample_metadata <- sample_metadata[keep_samples, , drop = FALSE]
  if (nrow(score_matrix) > 1) {
    score_matrix_norm <- t(scale(t(score_matrix)))
  } else {
    score_matrix_norm <- score_matrix
  }
  score_matrix_norm[is.na(score_matrix_norm)] <- 0
  score_matrix_norm[is.infinite(score_matrix_norm)] <- 0
  return(list(
    scores = score_matrix_norm,
    metadata = sample_metadata,
    original_scores = score_matrix
  ))
}
#' 修复的降维函数
perform_dimensionality_reduction_fixed <- function(score_matrix, n_components = 10) {
  results <- list()
  # 1. PCA (主成分分析) - 始终可用
  message("  Performing PCA...")
  tryCatch({
    # 转置：样本为行，基因/细胞类型为列
    data_for_pca <- t(score_matrix)
    # 检查数据是否适合PCA
    if (nrow(data_for_pca) > 1 && ncol(data_for_pca) > 1) {
      pca_result <- prcomp(data_for_pca, center = TRUE, scale. = TRUE)
      n_components_actual <- min(n_components, ncol(pca_result$x))
      results$pca <- list(
        embedding = pca_result$x[, 1:n_components_actual, drop = FALSE],
        variance = pca_result$sdev^2,
        loadings = pca_result$rotation[, 1:n_components_actual, drop = FALSE],
        sdev = pca_result$sdev)
      results$embedding <- results$pca$embedding
    } else {
      warning("Data not suitable for PCA, using original data")
      results$embedding <- data_for_pca
    }
  }, error = function(e) {
    warning(paste("PCA failed:", e$message))
    results$embedding <- t(score_matrix)
  })
  # 2. t-SNE (可选)
  if (require("Rtsne", quietly = TRUE)) {
    message("  Performing t-SNE...")
    tryCatch({
      tsne_result <- Rtsne::Rtsne(t(score_matrix), dims = 2, perplexity = min(30, ncol(score_matrix)/3),verbose = FALSE,max_iter = 1000)
      results$tsne <- tsne_result$Y
      rownames(results$tsne) <- colnames(score_matrix)
    }, error = function(e) {
      warning(paste("t-SNE failed:", e$message))
    })
  } else {
    message("  Rtsne package not available, skipping t-SNE")
  }
  # 3. UMAP (可选)
  if (require("umap", quietly = TRUE)) {
    message("  Performing UMAP...")
    tryCatch({
      umap_config <- umap::umap.defaults
      umap_config$n_components <- 2
      umap_config$n_neighbors <- min(15, ncol(score_matrix) - 1)
      umap_result <- umap::umap(t(score_matrix), config = umap_config)
      results$umap <- umap_result$layout
      rownames(results$umap) <- colnames(score_matrix)
    }, error = function(e) {
      warning(paste("UMAP failed:", e$message))
    })
  } else {
    message("  umap package not available, skipping UMAP")
  }
  # 4. MDS (多维尺度分析)
  message("  Performing MDS...")
  tryCatch({
    dist_matrix <- dist(t(score_matrix))
    mds_result <- cmdscale(dist_matrix, k = 2)
    results$mds <- mds_result
    rownames(results$mds) <- colnames(score_matrix)
  }, error = function(e) {
    warning(paste("MDS failed:", e$message))
  })
  
  return(results)
}
#' 修复的聚类函数
perform_clustering_fixed <- function(reduced_data, method = "kmeans", n_clusters = NULL) {
  if (is.null(reduced_data) || nrow(reduced_data) == 0) {
    stop("reduced_data is empty or NULL")}
  if (is.list(reduced_data) && "embedding" %in% names(reduced_data)) {
    reduced_data <- reduced_data$embedding}
  if (!is.matrix(reduced_data)) {
    reduced_data <- as.matrix(reduced_data)}
  if (is.null(n_clusters)) {
    n_clusters <- suggest_optimal_clusters_fixed(reduced_data)
    message(paste("  Auto-detected optimal clusters:", n_clusters))}
  if (n_clusters == 1) {
    clusters <- rep(1, nrow(reduced_data))
    names(clusters) <- rownames(reduced_data)
    return(list(
      clusters = clusters,
      n_clusters = 1,
      method = "none",
      quality_metrics = list(silhouette = 0, davies_bouldin = 0)
    ))
  }
  clusters <- switch(method,
    "kmeans" = {
      tryCatch({
        kmeans_result <- kmeans(reduced_data, centers = n_clusters, nstart = 25)
        clusters <- kmeans_result$cluster
        names(clusters) <- rownames(reduced_data)
        clusters
      }, error = function(e) {
        warning(paste("K-means failed:", e$message))
        # 返回单聚类
        rep(1, nrow(reduced_data))
      })
    },
    "hierarchical" = {
      tryCatch({
        dist_matrix <- dist(reduced_data)
        hclust_result <- hclust(dist_matrix)
        clusters <- cutree(hclust_result, k = n_clusters)
        names(clusters) <- rownames(reduced_data)
        clusters
      }, error = function(e) {
        warning(paste("Hierarchical clustering failed:", e$message))
        rep(1, nrow(reduced_data))
      })
    },
    "dbscan" = {
      if (require("dbscan", quietly = TRUE)) {
        tryCatch({
          # 自动确定eps
          k_dist <- sort(kNNdist(reduced_data, k = 5))
          eps <- k_dist[which.max(diff(k_dist))] + 0.1
          dbscan_result <- dbscan::dbscan(reduced_data, eps = eps, minPts = 5)
          clusters <- dbscan_result$cluster
          names(clusters) <- rownames(reduced_data)
          clusters
        }, error = function(e) {
          warning(paste("DBSCAN failed:", e$message))
          rep(1, nrow(reduced_data))
        })
      } else {
        warning("dbscan package not installed. Using kmeans.")
        kmeans_result <- kmeans(reduced_data, centers = n_clusters, nstart = 25)
        kmeans_result$cluster
      }
    },
    "louvain" = {
      if (require("igraph", quietly = TRUE)) {
        tryCatch({
          # 构建KNN图
          knn_graph <- build_knn_graph_fixed(reduced_data, k = min(20, nrow(reduced_data) - 1))
          louvain_result <- igraph::cluster_louvain(knn_graph)
          clusters <- louvain_result$membership
          names(clusters) <- rownames(reduced_data)
          clusters
        }, error = function(e) {
          warning(paste("Louvain clustering failed:", e$message))
          rep(1, nrow(reduced_data))
        })
      } else {
        warning("igraph package not installed. Using kmeans.")
        kmeans_result <- kmeans(reduced_data, centers = n_clusters, nstart = 25)
        kmeans_result$cluster
      }
    },
    "leiden" = {
      if (require("leiden", quietly = TRUE) && require("igraph", quietly = TRUE)) {
        tryCatch({
          knn_graph <- build_knn_graph_fixed(reduced_data, k = min(20, nrow(reduced_data) - 1))
          leiden_result <- leiden::leiden(knn_graph)
          clusters <- leiden_result
          names(clusters) <- rownames(reduced_data)
          clusters
        }, error = function(e) {
          warning(paste("Leiden clustering failed:", e$message))
          rep(1, nrow(reduced_data))
        })
      } else {
        warning("leiden or igraph package not installed. Using kmeans.")
        kmeans_result <- kmeans(reduced_data, centers = n_clusters, nstart = 25)
        kmeans_result$cluster
      }
    },
    {
      warning(paste("Method", method, "not recognized. Using kmeans."))
      kmeans_result <- kmeans(reduced_data, centers = n_clusters, nstart = 25)
      kmeans_result$cluster
    }
  )
  if (is.null(names(clusters))) {
    names(clusters) <- rownames(reduced_data)
  }
  quality_metrics <- calculate_clustering_quality_fixed(reduced_data, clusters)
  return(list(
    clusters = clusters,
    n_clusters = length(unique(clusters)),
    method = method,
    quality_metrics = quality_metrics
  ))
}
#' 修复的构建KNN图函数
build_knn_graph_fixed <- function(data, k = 20) {
  if (require("FNN", quietly = TRUE)) {
    k <- min(k, nrow(data) - 1)
    if (k < 1) {
      # 如果样本太少，创建完全图
      edges <- t(combn(1:nrow(data), 2))
      graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
      return(graph)
    }
    knn_result <- FNN::get.knn(data, k = k)
    edges <- matrix(0, nrow = 0, ncol = 2)
    for (i in 1:nrow(data)) {
      neighbors <- knn_result$nn.index[i, ]
      for (j in neighbors) {
        edges <- rbind(edges, c(i, j))
      }
    }
    graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
    return(graph)
  } else {
    dist_matrix <- as.matrix(dist(data))
    threshold <- quantile(dist_matrix, 0.1)  # 连接最近的10%
    adj_matrix <- dist_matrix <= threshold
    diag(adj_matrix) <- FALSE
    graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    return(graph)
  }
}
#' 修复的建议最优聚类数
suggest_optimal_clusters_fixed <- function(data, max_clusters = 10) {
  n_samples <- nrow(data)
  if (n_samples <= 3) {return(1)}
  max_clusters <- min(max_clusters, n_samples - 1)
  if (max_clusters < 2) {return(1)}
  if (require("cluster", quietly = TRUE)) {
    tryCatch({
      wss <- numeric(max_clusters)
      for (i in 1:max_clusters) {
        kmeans_result <- kmeans(data, centers = i, nstart = 10)
        wss[i] <- kmeans_result$tot.withinss
      }
      diff_wss <- diff(wss)
      relative_reduction <- diff_wss / wss[-length(wss)]
      optimal_k <- which.min(relative_reduction) + 1
      return(min(max(optimal_k, 2), max_clusters))
    }, error = function(e) {
      return(min(3, max_clusters))
    })
  } else {
    return(min(floor(sqrt(n_samples)), max_clusters))
  }
}
#' 修复的聚类质量计算
calculate_clustering_quality_fixed <- function(data, clusters) {
  metrics <- list()
  n_clusters <- length(unique(clusters))
  if (n_clusters == 1) {
    metrics$silhouette <- 0
    metrics$davies_bouldin <- 0
    return(metrics)
  }
  if (require("cluster", quietly = TRUE)) {
    tryCatch({
      # Silhouette score
      sil_score <- cluster::silhouette(clusters, dist(data))
      metrics$silhouette <- mean(sil_score[, 3], na.rm = TRUE)
      # Davies-Bouldin index
      metrics$davies_bouldin <- cluster::index.DB(data, clusters)$DB
    }, error = function(e) {
      metrics$silhouette <- NA
      metrics$davies_bouldin <- NA
    })
  }
  return(metrics)
}
#' 生成可视化
generate_visualizations_fixed <- function(result) {
  plots <- list()
  # 1. 组织特异性热图
  plots$tissue_specificity_heatmap <- plot_tissue_specificity_heatmap(tissue_mean_scores = result$tissue_specificity$tissue_mean_scores,tau_scores = result$tissue_specificity$tau_scores)
  # 2. 降维图
  plots$dim_reduction <- plot_dimension_reduction(dim_data = result$dim_reduction,metadata = result$processed_data$metadata)
  # 3. 聚类图
  plots$clustering <- plot_clustering_results(dim_data = result$dim_reduction,clusters = result$clustering$clusters,metadata = result$processed_data$metadata)
  # 4. 富集热图
  plots$enrichment_heatmaps <- plot_enrichment_heatmaps(tissue_enrichment = result$enrichment$cluster_tissue_enrichment,celltype_enrichment = result$enrichment$cluster_celltype_enrichment)
  # 5. 组织间相关性
  plots$tissue_correlation <- plot_tissue_correlation(correlation_matrix = result$tissue_specificity$tissue_correlation)
  return(plots)
}
#' 绘制组织特异性热图
plot_tissue_specificity_heatmap <- function(tissue_mean_scores, tau_scores) {
  if (!require("pheatmap", quietly = TRUE) || !require("RColorBrewer", quietly = TRUE)) {
    install.packages(c("pheatmap", "RColorBrewer"))}
  library(pheatmap)
  library(RColorBrewer)
  # 添加Tau分数作为行注释
  row_annotation <- data.frame(Tau = tau_scores)
  rownames(row_annotation) <- names(tau_scores)
  p <- pheatmap(tissue_mean_scores,main = "Tissue-specific Cell Type Scores",color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),annotation_row = row_annotation,show_rownames = TRUE,show_colnames = TRUE,fontsize_row = 8,fontsize_col = 10,cluster_rows = TRUE,cluster_cols = TRUE)
  return(p)
}
#' 绘制降维图
plot_dimension_reduction <- function(dim_data, metadata) {
  if (!require("ggplot2", quietly = TRUE) || !require("patchwork", quietly = TRUE)) {
    install.packages(c("ggplot2", "patchwork"))}
  library(ggplot2)
  library(patchwork)
  plots <- list()
  # PCA
  if (!is.null(dim_data$pca$embedding)) {
    pca_df <- as.data.frame(dim_data$pca$embedding[, 1:2])
    colnames(pca_df) <- c("PC1", "PC2")
    pca_df$Tissue <- metadata$tissue[match(rownames(pca_df), rownames(metadata))]
    plots$pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
      geom_point(size = 2, alpha = 0.7) +
      theme_bw() +
      ggtitle("PCA Plot by Tissue")
  }
  # t-SNE
  if (!is.null(dim_data$tsne)) {
    tsne_df <- as.data.frame(dim_data$tsne)
    colnames(tsne_df) <- c("tSNE1", "tSNE2")
    tsne_df$Tissue <- metadata$tissue[match(rownames(tsne_df), rownames(metadata))]
    plots$tsne <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Tissue)) +
      geom_point(size = 2, alpha = 0.7) +
      theme_bw() +
      ggtitle("t-SNE Plot by Tissue")
  }
  # UMAP
  if (!is.null(dim_data$umap)) {
    umap_df <- as.data.frame(dim_data$umap)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$Tissue <- metadata$tissue[match(rownames(umap_df), rownames(metadata))]
    plots$umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Tissue)) +
      geom_point(size = 2, alpha = 0.7) +
      theme_bw() +
      ggtitle("UMAP Plot by Tissue")
  }
  return(plots)
}
#' 绘制聚类结果
plot_clustering_results <- function(dim_data, clusters, metadata) {
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
  }
  library(ggplot2)
  plots <- list()
  # 使用UMAP坐标
  if (!is.null(dim_data$umap)) {
    plot_df <- as.data.frame(dim_data$umap)
    colnames(plot_df) <- c("UMAP1", "UMAP2")
    plot_df$Cluster <- as.factor(clusters[match(rownames(plot_df), names(clusters))])
    plot_df$Tissue <- metadata$tissue[match(rownames(plot_df), rownames(metadata))]
    # 按聚类着色
    p1 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
      geom_point(size = 2, alpha = 0.7) +
      theme_bw() +
      ggtitle("Clustering Results")
    # 按组织着色
    p2 <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = Tissue)) +
      geom_point(size = 2, alpha = 0.7) +
      theme_bw() +
      ggtitle("Tissue Distribution")
    plots$clustering_umap <- list(by_cluster = p1, by_tissue = p2)
  }
  return(plots)
}
#' 绘制富集热图
plot_enrichment_heatmaps <- function(tissue_enrichment, celltype_enrichment) {
  if (!require("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")}
  library(pheatmap)
  plots <- list()
  # 组织富集热图
  if (!is.null(tissue_enrichment) && nrow(tissue_enrichment) > 1) {
    plots$tissue_enrichment <- pheatmap(tissue_enrichment,main = "Cluster-Tissue Enrichment (-log10 p-value)",color = colorRampPalette(c("white", "red"))(100),show_rownames = TRUE,show_colnames = TRUE,cluster_rows = TRUE,cluster_cols = TRUE)
  }
  # 细胞类型富集热图
  if (!is.null(celltype_enrichment) && nrow(celltype_enrichment) > 1) {
    plots$celltype_enrichment <- pheatmap(celltype_enrichment,main = "Cell Type Enrichment in Clusters",color = colorRampPalette(c("blue", "white", "red"))(100),show_rownames = TRUE,show_colnames = TRUE,cluster_rows = TRUE,cluster_cols = TRUE)
  }
  return(plots)
}
#' 绘制组织间相关性
plot_tissue_correlation <- function(correlation_matrix) {
  if (!require("corrplot", quietly = TRUE)) {
    install.packages("corrplot")}
  library(corrplot)
  p <- corrplot(correlation_matrix,method = "color",type = "upper",order = "hclust",tl.col = "black",tl.srt = 45,main = "Tissue Correlation Matrix",mar = c(0, 0, 2, 0))
  return(p)
}
#' 保存分析结果
save_analysis_results <- function(result, output_dir = "tissue_enrichment_results") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)}
  saveRDS(result, file.path(output_dir, "tissue_enrichment_analysis.rds"))
  write.csv(result$tissue_specificity$tissue_mean_scores,file.path(output_dir, "tissue_mean_scores.csv"))
  write.csv(result$tissue_specificity$tissue_correlation,file.path(output_dir, "tissue_correlation.csv"))
  write.csv(data.frame(CellType = names(result$tissue_specificity$tau_scores), Tau = result$tissue_specificity$tau_scores),file.path(output_dir, "tau_scores.csv"))
  pdf(file.path(output_dir, "visualizations.pdf"), width = 12, height = 10)
  for (plot_type in names(result$plots)) {
    if (inherits(result$plots[[plot_type]], "list")) {
      for (subplot in names(result$plots[[plot_type]])) {
        print(result$plots[[plot_type]][[subplot]])
      }
    } else {
      print(result$plots[[plot_type]])
    }
  }
  dev.off()
  message(paste("Results saved to:", output_dir))
}
result_fixed <- tissue_enrichment_analysis_fixed(
  score_matrix = score_matrix,
  sample_metadata = sample_metadata,
  n_components = 5,
  clustering_method = "kmeans",
  n_clusters = 5)
save_analysis_results(result_fixed, output_dir = "2.Tissue_enrichment")
################################## END ##################################
