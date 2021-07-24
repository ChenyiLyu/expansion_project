capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

# visualization settings

FeaturePlotScoring<- function(obj, feature, metadata_column, ...){
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal<- min(obj@meta.data[, feature])
  maximal<- max(obj@meta.data[, feature])
  p<- FeaturePlot(obj, features = feature,  coord.fixed = 1, min.cutoff = 0, 
                  pt.size = 1, repel = TRUE,
                  split.by =metadata_column, order = TRUE) & scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"),
                                                                                    limits=c(0, maximal)) & theme_bw(base_line_size = 0)
  return(p)
}

theme_hm_DE <-theme(axis.text.y = element_text(size=5),
                    plot.title = element_text(color="black", size=14, face="bold"),
                    legend.text = element_text(size = 14),
                    legend.title = element_text(size=14))

theme_fp_col <- scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlGn")))
theme_fp_score_col <- scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Blues"))

theme_umap_col <- brewer.pal(n=12,(palette = 'Paired'))
theme_bp <-theme(axis.text = element_text(size=14),
                 legend.text = element_text(size = 14),
                 legend.title = element_text(size=14))

theme_fp <- theme_bw() & NoLegend() & NoAxes()
# function
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}
get_conserved <- function(cluster){
  # Create function to get conserved markers for any given cluster
  FindConservedMarkers(seurat_integrated,
                       ident.1 = cluster,
                       grouping.var = "sample",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

venn_plot <- function(set1, set2, set3, reg, file){
  if (reg=='up') {
    col = alpha("#440154ff",0.3)
  }else{
    col = alpha('#fde725ff',0.3)
  }
  # Chart
  venn <- venn.diagram(
    x = list(set1, set2, set3),
    na = 'remove',
    category.names = c(deparse(substitute(set1)), deparse(substitute(set2)), deparse(substitute(set3))),
    filename = NULL,
    fill= col,
    
    units = "px",
    height = 100,
    width = 100,
    resolution = 100,
    output=TRUE,
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
  pdf(file)
  grid.draw(venn)
  dev.off()
  overlap <- intersect(set1,intersect(set2,set3))
  overlap <- overlap[!is.na(overlap)]
  return(overlap)
}

# parameter
position <- c('Spinous I','Spinous II','Spinous III','Prolif_basal','Basal','SG','Infundibulum','HF III', 'HF IV', 'HF V', 'HFSC','Immune_cell')
# @source https://github.com/statOmics/tradeSeq/blob/master/R/evaluateK.R
#' Evaluate an appropriate number of knots.
#'
#' @param aicMat The output from \code{\link[tradeSeq]{evaluateK}}
#' @param k The range of knots to evaluate. `3:10` by default. Extracted from 
#' the column names by default
#' @param aicDiff Used for selecting genes with significantly varying AIC values
#' over the range of evaluated knots to make the barplot output. Default is set
#' to 2, meaning that only genes whose AIC range is larger than 2 will be used
#' to check for the optimal number of knots through the barplot visualization
#' that is part of the output of this function.
#' @examples
#' ## This is an artificial example, please check the vignette for a realistic one.
#' set.seed(8)
#' data(sds, package="tradeSeq")
#' loadings <- matrix(runif(2000*2, -2, 2), nrow = 2, ncol = 2000)
#' counts <- round(abs(t(slingshot::reducedDim(sds) %*% loadings))) + 100
#' aicK <- evaluateK(counts = counts, sds = sds, nGenes = 100,
#'                   k = 3:5, verbose = FALSE, plot = FALSE)
#' plot_evalutateK_results(aicK, k = 3:5)
#' @export
plot_evalutateK_results <- function(aicMat, k = NULL, aicDiff = 2) {
  if (is.null(k)) {
    k <- as.numeric(sub("k: ", "", colnames(aicMat)))
  }
  
  init_shape <- graphics::par()$mfrow
  graphics::par(mfrow = c(1, 4))
  # boxplots of AIC
  devs <- matrix(NA, nrow = nrow(aicMat), ncol = length(k))
  for (ii in seq_len(nrow(aicMat))) {
    devs[ii, ] <- aicMat[ii, ] - mean(aicMat[ii, ])
  }
  graphics::boxplot(devs, ylab = "Deviation from genewise average AIC",
                    xlab = "Number of knots", xaxt = "n")
  graphics::axis(1, at = seq_len(length(k)), labels = k)
  # scatterplot of average AIC
  plot(x = k, y = colMeans(aicMat, na.rm = TRUE), type = "b",
       ylab = "Average AIC", xlab = "Number of knots")
  # scatterplot of relative AIC
  plot(x = k, y = colMeans(aicMat / aicMat[, 1], na.rm = TRUE), type = "b",
       ylab = "Relative AIC", xlab = "Number of knots")
  # barplot of optimal AIC for genes with at least a difference of 2.
  aicRange <- apply(apply(aicMat, 1, range), 2, diff)
  varID <- which(aicRange > aicDiff)
  if (length(varID) > 0) {
    aicMatSub <- aicMat[varID, ]
    tab <- table(k[apply(aicMatSub, 1, which.min)])
    graphics::barplot(tab, xlab = "Number of knots",
                      ylab = "# Genes with optimal k")
  }
  graphics::par(mfrow = init_shape)
  return()
}













