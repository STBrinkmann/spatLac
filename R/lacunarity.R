#' Lacunarity
#' @description Computes the Lacunarity of a Spatial Raster object or a list of Spatial Raster objects.
#'
#' @param r SpatRaster or list of SpatRaster objects; Raster image can be loaded with \code{\link[terra]{rast}}
#' @param box character; Either SQUARE for a square neighborhood window, or CIRCLE for a round neighborhood window
#' @param plot logical; Should the summary Lacunarity figure, showing Lacunarity of all SpatRasters be printed?
#' @param save_plot FALSE or folder path; If not FALSE, a folder path to save Lacunarity plots (see details)
#' @param progress logical; Show progress bar?
#' @param ncores numeric; The number of cores to use
#'
#' @details
#' \code{lacunarity} is based on the algorithm for binary images provided by \href{https://doi.org/10.1007/BF00125351}{Plotnick et al. 1993}.
#' \href{https://doi.org/10.1016/j.ecocom.2011.01.001}{Hoechstetter et al. 2011} further applyed this algorithm on continuous raster images.
#' If \code{r} is a binary raster (e.g. Land Use) Lacunarity is beening calculated using the Plotnicks algorithm.
#'
#' @return \code{\link[dplyr]{tibble}} containing all Lacunarity values
#' @export
#'
#' @examples
#' # Create a SpatRast as input
#' mat_sample <- matrix(data = c(
#'    1,1,0,1,1,1,0,1,0,1,1,0,
#'    0,0,0,0,0,1,0,0,0,1,1,1,
#'    0,1,0,1,1,1,1,1,0,1,1,0,
#'    1,0,1,1,1,0,0,0,0,0,0,0,
#'    1,1,0,1,0,1,0,0,1,1,0,0,
#'    0,1,0,1,1,0,0,1,0,0,1,0,
#'    0,0,0,0,0,1,1,1,1,1,1,1,
#'    0,1,1,0,0,0,1,1,1,1,0,0,
#'    0,1,1,1,0,1,1,0,1,0,0,1,
#'    0,1,0,0,0,0,0,0,0,1,1,1,
#'    0,1,0,1,1,1,0,1,1,0,1,0,
#'    0,1,0,0,0,1,0,1,1,1,0,1
#'    ), nrow = 12, ncol = 12, byrow = TRUE
#' )
#' r <- terra::rast(mat_sample)
#'
#' lacunarity(r)
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom dplyr tibble
#' @importFrom dplyr bind_rows
#' @importFrom dplyr arrange
#' @importFrom terra unique
#' @importFrom terra as.matrix
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggsave
#' @importFrom magrittr %>%
#' @useDynLib spatLac, .registration = TRUE
lacunarity <- function(r, box = "SQUARE", plot = FALSE, save_plot = FALSE, progress = FALSE, ncores = 1L) {
  
  # 1. Check input ----------------------------------------------------------
  # r
  if (class(r) == "list") {
    if (!all((lapply(r, class) == "SpatRaster"))) {
      stop("all elements of r must  be SpatRaster objects. Make sure to use the terra package for loading raster images.")
    }
  } else {
    if (class(r) != "SpatRaster") {
      stop("r must  be a SpatRaster object. Make sure to use the terra package for loading raster images.")
    }
  }
  
  #box
  if (!is.character(box)) {
    stop("box must be character")
  } else if (toupper(box) == "SQUARE") {
    box <- 1
  } else if (toupper(box) == "CIRCLE") {
    box <- 2
  } else {
    stop("box must be either SQUARE or CIRCLE")
  }
  
  # plot
  if (!is.logical(plot)) {
    stop("plot must be logical (TRUE/FALSE)")
  }
  
  # save_plot
  if(is.logical(save_plot) && save_plot) {
    stop("if you wish to save the output plots, provide a folder path")
  } else if(is.character(save_plot)) {
    # Tempdir
    temp_path <- tempfile(paste0("Lacunarity_", gsub("-", "_", Sys.Date()), "__"), save_plot)
    dir.create(temp_path, recursive = TRUE)
  }
  
  # progress
  if (!is.logical(progress)) {
    stop("progress must be logical (TRUE/FALSE)")
  }
  
  # ncores
  if (!is.numeric(ncores)) {
    stop("ncores must be numeric")
  } else {
    ncores <- ifelse(ncores<1, 1L, floor(ncores))
  }
  
  
  # 2. Main loop ------------------------------------------------------------
  
  # Output tibble
  out <- dplyr::tibble(
    name = as.character(),
    i = as.integer(),
    r = as.integer(),
    "ln(r)" = as.integer(),
    Lac = as.numeric(),
    "ln(Lac)" = as.numeric()
  )
  
  # Convert r to list
  if(!is.list(r)) r <- list(r)
  
  for (i in seq_along(r)) {
    # Get current r
    this_r <- r[[i]]
    if (progress) {
      cat(paste("Computing Lacunarity for", names(this_r)[1]))
      cat("\n")
    }
    
    # Calculate r vector from raster dimension
    max_r <- 1
    while (2^(max_r+1)+1 < min(dim(this_r)[1:2])/2) {
      max_r <- max_r + 1
    }
    
    r_vec <- c(2^(1:max_r)+1, floor(min(dim(this_r)[1:2])/2))
    
    # Is r binary? (e.g. Greenspace raster)
    lac_fun <- as.integer(nrow(terra::unique(this_r)) <= 2)
    
    # Convert raster to matrix
    rast_mat <- terra::as.matrix(this_r, wide = TRUE)
    
    # Calculate Lacunarity for all w
    this_lac <- rcpp_lacunarity(mat = rast_mat,
                                w_vec = r_vec,
                                fun = lac_fun,
                                mode = box,
                                ncores = ncores,
                                display_progress = progress)
    
    this_out <- dplyr::tibble(
      name = rep(names(this_r)[1], length(this_lac)),
      i = i,
      r = r_vec,
      "ln(r)" = log(r_vec),
      Lac = this_lac,
      "ln(Lac)" = log(this_lac)
    )
    
    out <- dplyr::bind_rows(out, this_out)
    if(progress) cat("\n")
  }
  
  
  # Plot output of all r
  if (plot) {
    if (length(unique(out$name)) != length(unique(out$i))) {
      r_names <- as.character(out$i)
    } else {
      r_names <- as.character(out$name)
    }
    
    max_x <- ceiling(max(out$`ln(r)`, na.rm = TRUE))
    max_y <- ceiling(max(out$`ln(Lac)`, na.rm = TRUE))
    
    p <- ggplot2::ggplot(data = out,
                         mapping = ggplot2::aes(x = `ln(r)`,
                                                y = `ln(Lac)`,
                                                colour = r_names),
                         environment = environment()) +
      ggplot2::geom_line(lwd = 0.8) +
      ggplot2::labs(x = "ln(r)",
                    y = "ln(lacunarity)",
                    colour = "Legend") +
      ggplot2::scale_x_continuous(breaks = seq(0, max_x, 1),
                                  limits = c(1, max_x)) +
      ggplot2::scale_y_continuous(breaks = seq(0, max_y, 0.1),
                                  limits = c(0, max_y)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                     axis.text = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 12),
                     legend.title = ggplot2::element_text(size = 14))
    
    print(p)
  }
  
  
  # Save plots and table of out
  if (is.character(save_plot)) {
    # Save out as CSV
    write.csv(out, file = file.path(temp_path, "Lacunarity.csv"), row.names = FALSE)
    
    # Save summary ggplot
    if (length(unique(out$name)) != length(unique(out$i))) {
      r_names <- as.character(out$i)
    } else {
      r_names <- as.character(out$name)
    }
    
    max_x <- ceiling(max(out$`ln(r)`, na.rm = TRUE))
    max_y <- ceiling(max(out$`ln(Lac)`, na.rm = TRUE))
    
    p <- ggplot2::ggplot(data = out,
                         mapping = ggplot2::aes(x = `ln(r)`,
                                                y = `ln(Lac)`,
                                                colour = r_names),
                         environment = environment()) +
      ggplot2::geom_line(lwd = 0.8) +
      ggplot2::labs(x = "ln(r)", y = "ln(lacunarity)",
                    colour = "Legend") +
      ggplot2::scale_x_continuous(breaks = seq(0, max_x, 1), limits = c(1, max_x)) +
      ggplot2::scale_y_continuous(breaks = seq(0, max_y, 0.1), limits = c(0, max_y)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                     axis.text = ggplot2::element_text(size = 12),
                     legend.text = ggplot2::element_text(size = 12),
                     legend.title = ggplot2::element_text(size = 14))
    
    ggplot2::ggsave(filename = file.path(temp_path, "Lacunarity_summary.svg"), 
                    plot = p, device = "svg", 
                    width = 7, height = 6.5, units = "in")
    ggplot2::ggsave(filename = file.path(temp_path, "Lacunarity_summary.png"), 
                    plot = p, device = "png", 
                    width = 7, height = 6.5, units = "in", dpi = 500)
    
    # Save unique plots
    if (length(r) > 1) {
      for (i in seq_along(r)) {
        this_out <- out[out$i == i, ]
        
        if (length(unique(this_out$name)) != length(unique(this_out$i))) {
          r_names <- as.character(this_out$i)
        } else {
          r_names <- as.character(this_out$name)
        }
        
        max_x <- ceiling(max(this_out$`ln(r)`, na.rm = TRUE))
        max_y <- ceiling(max(this_out$`ln(Lac)`, na.rm = TRUE))
        
        p <- ggplot2::ggplot(data = this_out,
                             mapping = ggplot2::aes(x = `ln(r)`,
                                                    y = `ln(Lac)`,
                                                    colour = r_names),
                             environment = environment()) +
          ggplot2::geom_line(lwd = 0.8) +
          ggplot2::labs(x = "ln(r)", y = "ln(lacunarity)",
                        colour = "Legend") +
          ggplot2::scale_x_continuous(breaks = seq(0, max_x, 1), limits = c(1, max_x)) +
          ggplot2::scale_y_continuous(breaks = seq(0, max_y, 0.1), limits = c(0, max_y)) +
          ggplot2::theme_minimal() +
          ggplot2::theme(axis.title = ggplot2::element_text(size = 14),
                         axis.text = ggplot2::element_text(size = 12),
                         legend.text = ggplot2::element_text(size = 12),
                         legend.title = ggplot2::element_text(size = 14))
        
        ggplot2::ggsave(filename = file.path(temp_path, paste0("Lacunarity_", i, ".svg")), 
                        plot = p, device = "svg", 
                        width = 7, height = 6.5, units = "in")
        ggplot2::ggsave(filename = file.path(temp_path, paste0("Lacunarity_", i, ".png")), 
                        plot = p, device = "png", 
                        width = 7, height = 6.5, units = "in", dpi = 500)
      } 
    }
  }
  return(dplyr::arrange(out, i))
}