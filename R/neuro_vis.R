#' 2D visualization plot of a neuroimaGene object
#'
#' Generates a 2D visualization plot of the neuroimaGene object. Neuroimaging
#' regions are defined by the atlas parameter and colored according to the
#' magnitude and direction of the aggregate effect from each gene in the
#' NeuroimaGene object. Colors can be defined by the user.
#'
#' @param ng_obj NeuroimaGene object produced by neuroimaGene() function
#' @param atlas desired atlas for visualization. Desikan (default), Subcortex, DKT, Destrieux.
#' @param lowcol color for low end of Zscore spectrum. Default is red
#' @param midcol color for middle of Zscore spectrum. Default is white
#' @param highcol color for top end of Zscore spectrum. Default is blue
#' @param title optional title tag for the plot
#' @keywords neuroimaging
#' @export
#' @import data.table ggplot2 stringr ggseg sf
#' @importFrom utils write.table
#' @importFrom stats na.omit
#' @returns class: ggplot object depicting 2D visualization of the NIDPs from the neuroimaGene object portrayed on the brain and shaded by mean effect size.
#' @examples
#' gene_list <- c('TRIM35', 'PROSER3', 'EXOSC6', 'PICK1', 'UPK1A', 'ESPNL', 'ZIC4')
#' ng <- neuroimaGene(gene_list, atlas = NA, mtc = 'BH', vignette = TRUE)
#' neuro_vis(ng, atlas = 'DKT')
#'
#'
neuro_vis <- function(ng_obj, atlas = 'Desikan', lowcol = 'red2', midcol = 'white', highcol = 'royalblue2', title = NA) {
  # initialize column names as null variables
  zscore <- atlasnm <- atl <- gwas_phenotype <- meanZ <- measurement <- x <- y <- label <- NULL
  
  # load required local data from package
  dkt_atl <- readRDS(system.file("extdata","dkt_atlas_new.rda", package = "neuroimaGene"))
  dest_atl <- readRDS(system.file("extdata","dest_atlas_new.rda", package = "neuroimaGene"))
  
  if(is.na(title)){
    tag <- ''
  } else {
    tag <- paste(' ',as.character(title))
  }
  ng_summ <- ng_obj[, list(meanZ = mean(zscore), sign = sign(mean(zscore))),
                    by = c('gwas_phenotype')]
  ng <- data.table::setDT(merge(ng_summ, anno[, c('gwas_phenotype', 'measurement')], by = 'gwas_phenotype'))
  ng <- data.table::setDT(merge(ng, fs_anno, by = 'gwas_phenotype'))
  
  
  atldir <- data.table::data.table(atlasnm = c('Desikan', 'DKT', 'Destrieux', 'Subcortex',
                                               'desikan', 'dkt', 'destrieux', 'subcortex'),
                                   realname =c('Desikan', 'DKT', 'Destrieux', 'Subcortex',
                                               'Desikan', 'DKT', 'Destrieux', 'Subcortex'),
                                   fsatl = c('ggseg::dk()', 'dkt_atl',
                                             'dest_atl','ggseg::aseg()',
                                             'ggseg::dk()', 'dkt_atl',
                                             'dest_atl','ggseg::aseg()'),
                                   fsnm = c('dk', 'dkt', 'destrieux', 'aseg',
                                            'dk', 'dkt', 'destrieux', 'aseg'))
  atlname = atldir[atlasnm == atlas,]$realname
  fs = atldir[atlasnm == atlas,]$fsatl
  fs2 = atldir[atlasnm == atlas,]$fsnm
  stat_vis <- stats::na.omit(ng[atl == atlname,])
  
  if (dim(stat_vis)[1] == 0) {
    stop(paste0('No nidps from the',atlas,'atlas detected.'))
  }
  
  if (atlas == 'Subcortex' ) {
    # plot aseg volumes
    aseg_vol <- stat_vis[gwas_phenotype %like% 'volume' & atl == 'Subcortex',]
    brain_plot2 <- ggplot(aseg_vol) +
      geom_brain(atlas = aseg(), aes(fill = meanZ), position = position_brain(nrow = 2)) +
      ggtitle(paste0('Subcortical NIDPs (aseg atlas)', tag))+
      scale_fill_gradient2(low = lowcol, mid = midcol, high = highcol, na.value = "lightgrey", name = "meanZ") +
      annotate_brain(atlas = aseg(), position = position_brain(nrow = 2), colour = "grey30",
                     family = "sans") +
      theme_void()
    
  } else if (atlas == 'Desikan'){
    
    measures = unique(stat_vis$measurement[!is.na(stat_vis$measurement)])
    stat_vis_final <- data.table::data.table()
    for(msr in measures) {
      temp1 <- data.table::as.data.table(merge(stat_vis[measurement == msr,],
                                               data.table::as.data.table(dk()$data$vertices)[,c('label', 'vertices')],
                                               by = c('label'),
                                               all.y = TRUE))
      
      temp1$measurement <- msr
      stat_vis_final = rbind(stat_vis_final, temp1)
    }
    stat_vis_final$atlas <- fs2
    
    pos <- position_brain(hemi ~ view)
    brain_plot <- ggplot(stat_vis_final) +
      geom_brain(atlas = dk(), aes(fill = meanZ), position = pos) +
      scale_fill_gradient2(low = lowcol, mid = midcol, high = highcol, na.value = "lightgrey", name = "meanZ") +
      facet_wrap(~measurement, nrow = length(unique(stat_vis_final$measurement))) +
      theme_void()
    
    bp <- ggplot_build(brain_plot)
    top_x <- max(bp$layout$panel_params[[1]]$x$breaks)
    top_y <- max(bp$layout$panel_params[[1]]$y$breaks)
    
    col_labels <- data.frame(
      x = c((top_x/8), (3*top_x)/8, (5*top_x)/8, (7*top_x)/8),  
      y = -top_y/20,                  
      label = c("inferior", "lateral", "medial", "superior")
    )
    row_labels <- data.frame(
      x = 0 ,#-top_x/20,                 
      y = c(top_y/4, (3*top_y)/4),        
      label = c("right", "left")
    )
    
    brain_plot2 <- ggplot(stat_vis_final) +
      geom_brain(atlas = dk(), aes(fill = meanZ), position = pos) +
      geom_text(data = col_labels, aes(x = x, y = y, label = label),
                size = 3, colour = "grey30",
                family = "sans") +
      geom_text(data = row_labels, aes(x = x, y = y, label = label),
                size = 3,  colour = "grey30",
                family = "sans", angle = 90, vjust =-1) +
      ggtitle(paste0(atlas, ' atlas NIDPs', tag))+
      scale_fill_gradient2(low = lowcol, mid = midcol, high = highcol, na.value = "lightgrey", name = "meanZ") +
      facet_wrap(~measurement, nrow = length(unique(stat_vis_final$measurement)), strip.position="left") +
      theme_void() +
      theme(strip.text = element_text(face = "bold", family = 'sans', size = 10, angle = 90),
            strip.background = element_rect(fill="grey90"))
    
  } else {
    
    measures = unique(stat_vis$measurement[!is.na(stat_vis$measurement)])
    atlas_joined <- sf::st_sf(geometry = sf::st_sfc())
    for(msr in measures) {
      temp <- data.table::as.data.table(merge(stat_vis[measurement == msr,],
                                              data.table::as.data.table(eval(parse(text = fs))$data$sf),
                                              by = c('label'),
                                              all.y = TRUE))
      temp$measurement <- msr
      atlas_joined = rbind(atlas_joined, sf::st_as_sf(temp))
    }
    atlas_joined$atlas <- fs2
    
    brain_plot <- ggplot(atlas_joined) +
      geom_sf(aes(fill = meanZ)) +
      scale_fill_gradient2(low = lowcol, mid = midcol, high = highcol,
                           na.value = "lightgrey", name = "meanZ") +
      facet_wrap(~measurement, nrow = length(unique(measures)), strip.position="left") +
      theme_void()
    
    bp <- ggplot_build(brain_plot)
    top_x <- max(bp$layout$panel_params[[1]]$x$breaks)
    top_y <- max(bp$layout$panel_params[[1]]$y$breaks)
    
    col_labels <- data.frame(
      x = c((top_x/8), (3*top_x)/8, (5*top_x)/8, (7*top_x)/8),  # x center of each column
      y = -top_y/10,                  # just above the top row
      label = c("L lateral", "L medial", "R medial", "R Lateral")
    )
    
    brain_plot2 <- ggplot(atlas_joined) +
      geom_sf(aes(fill = meanZ)) +
      geom_text(data = col_labels, aes(x = x, y = y, label = label),
                size = 3, colour = "grey30",
                family = "sans") +
      ggtitle(paste0(atlas, ' atlas NIDPs', tag))+
      scale_fill_gradient2(low = lowcol, mid = midcol, high = highcol,
                           na.value = "lightgrey", name = "meanZ") +
      facet_wrap(~measurement, nrow = length(unique(measures)), strip.position="left") +
      theme_void() +
      theme(strip.text = element_text(face = "bold", family = 'sans', size = 10, angle = 90),
            strip.background = element_rect(fill="grey90"))
    
    
  }
  return(brain_plot2)
}
