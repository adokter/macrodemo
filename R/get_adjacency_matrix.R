#' Get adjacency matrix
#' @param data a vector of cell indices
#' @param grid a ddgridR grid
#' @param n_small_min minimum number of small cells required to calculate an abundance ratio for a large cell
#' @param quiet if TRUE, suppress informational messages (warnings/errors still returned)
#' @export
#' @return a square matrix indicating which cells are direct neighbors (value 1) and which cells are not (value 0)
get_adjacency_matrix <- function(data, grid_large){
  # tolerance for when cells are neighbors
  spacing_tolerance=1.2
  # get the grid spacing in m
  grid_spacing <- 1000*(dggetres(grid_large) %>% filter(res==grid_large$res) %>% pull(spacing_km))
  # initialize adjacency matrix
  adjacency_mat <- matrix(data = 0L, nrow = length(data), ncol = length(data))
  row.names(adjacency_mat) <- data
  for (i in 1:(length(data) - 1)) {
    coord_i <- dggridR::dgSEQNUM_to_GEO(grid_large, data[i])
    for (j in (i+1):length(data)) {
      coord_j <- dggridR::dgSEQNUM_to_GEO(grid_large, data[j])
      cell_dist <- geosphere::distm(c(coord_i$lon_deg, coord_i$lat_deg),
                                    c(coord_j$lon_deg, coord_j$lat_deg),
                                    fun = geosphere::distHaversine)
      if(cell_dist < spacing_tolerance*grid_spacing) {
        adjacency_mat[i,j] <- adjacency_mat[j,i] <- 1L
      }
    }
  }
  adjacency_mat
}
