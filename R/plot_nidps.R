#' NIDP summary plot of NeuroimaGene object
#'
#' Generate overview plot of the neuroimagene object according to nidps
#' @param ng_obj NeuroimaGene Object
#' @param maxNidps maximum number of NIDPs to visualize. default=30
#' @param title optional title tag for the plot
#' @param verbose print runtime messages to R console. Default to FALSE
#' @keywords neuroimaging
#' @export
#' @import data.table ggplot2
#' @returns a ggplot class object detailing mean effect size magnitude per NIDP, colored by brain region
#' @examples
#' gene_list <- c('TRIM35', 'PROSER3', 'EXOSC6', 'PICK1', 'UPK1A', 'ESPNL', 'ZIC4')
#' ng <- neuroimaGene(gene_list, atlas = NA, mtc = 'BH', vignette = TRUE)
#' plot_nidps(ng)
#'
plot_nidps <- function(ng_obj, maxNidps = 30, title = NA, verbose = FALSE) {
  # initialize column names as null variables
  zscore <- meanZ <- gwas_phenotype <- secondary <- NULL

  if(is.na(title)){
    tag <- ''
  } else {
    tag <- paste(' ',as.character(title))
  }
  if(!(is.integer(maxNidps) || is.double(maxNidps))) {
    stop('maxNidps must be of data.type: integer or double', call. = F)
  }
  ng_summ <- ng_obj[, list(meanZ = mean(zscore), sign = sign(mean(zscore))),
                by = c('gwas_phenotype')]

  ng <- data.table::setDT(merge(ng_summ, anno, by = 'gwas_phenotype'))

  if (length(unique(ng$gwas_phenotype)) > maxNidps ) {
    if(verbose){message(paste('WARNING: Greater than', maxNidps, 'NIDPs detected in input data. Plot will only show the top', maxNidps, 'NIDPs ranked by effect size magnitude'))}
    nidps <- ng[order(-abs(meanZ))][1:maxNidps,]$gwas_phenotype
    ng <- ng[gwas_phenotype %in% nidps,]
  }

  gn_plot <- ggplot2::ggplot(ng, aes(x = gwas_phenotype, y = abs(meanZ), color= secondary, group = as.character(sign))) +
    geom_point(aes(shape=as.character(sign)), size = 4) + #size=as.character(gn_ct))) +
    theme_light()+
    ggtitle(paste0("Mean Effect Size per NIDP across\nall genes", tag)) +
    xlab('NIDPs') +
    ylab('Normalized effect size magnitude') +
    theme_light() +
    theme(axis.text.x = element_text(angle = 0, size = 11, hjust = 0.5, vjust =0.5),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          plot.title.position = "plot",
    ) +
    scale_shape_discrete(name  ="sign",
                         breaks=c("1", "-1"),
                         labels=c("pos", "neg")) +
    coord_flip()

  return(gn_plot)
}
