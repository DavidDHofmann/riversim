################################################################################
#### Load Dependencies
################################################################################
#' @import terra igraph
#' @importFrom deldir deldir
NULL

################################################################################
#### Level 1 Functions
################################################################################
#' Simulate a River Network
#'
#' This function is used to simulate a river network.
#' @export
#' @param n_points integer indicating the number of randomly generated
#' points/nodes. The larger the number of points/nodes, the more complex the
#' resulting rivernetwork. Set to 1000 by default.
#' @param xmn numeric minimum x-coordiante (left border)
#' @param ymn numeric minimum y-coordiante (bottom border)
#' @param scl numeric scaling factor by which the final rivernetwork should be
#' scaled. By default, the network will have a square extent of size 1 x 1.
#' @param progressbar logical indicating if a progressbar should be printed. Set
#' to TRUE by default.
#' @param raster logical indicating if a raster should be returned (i.e. if the
#' river network should be rasterized)
#' @param nrow integer number of rows in case a raster should be returned. 250
#' by default.
#' @param ncol integer number of columns in case a raster should be returned.
#' 250 by default.
#' @return \code{SpatialPolygonsDataFrame} of the river network
#' @examples
#' # Generate rivernetwork
#' set.seed(123)
#' riv <- rivernetwork(n_points = 500)
#'
#' # Visualize it
#' plot(riv, col = "cornflowerblue", border = NA)
#' axis(1)
#' axis(2)
#'
#' # Repeat, but returning a raster instead of SPDF
#' riv <- rivernetwork(n_points = 500, raster = T)
#' plot(riv, col = c("white", "cornflowerblue"))
rivernetwork <- function(n_points = 1000, xmn = 0, ymn = 0, scl = 1
  , progressbar = T, raster = F, nrow = 250, ncol = 250) {

  # Testing
  # library(igraph)
  # library(deldir)
  # library(raster)
  # library(rgeos)
  # n_points <- 1000
  # xmn <- 0
  # ymn <- 0
  # scl <- 1
  # progressbar <- T
  # raster <- T
  # nrow <- 250
  # ncol <- 250
  # scl <- 3

  # Sample points
  pts <- .samplePoints(n = n_points, deg = 4)

  # Initiate a couple of vectors
  flowq    <- order(pts$z, decreasing = T)
  qtop     <- flowq[1]
  qend     <- flowq[length(flowq)]
  rechecks <- c()
  visited  <- c() # Can be deleted here?
  vetos    <- c() # Can be deleted here?

  # Compute the delaunay triangulation / voronoi polygons
  deltri <- deldir(pts)
  hull   <- chull(pts)

  # Generate empty network
  net <- make_empty_graph(directed = T)
  E(net)$outflow_x <- NA
  E(net)$outflow_y <- NA
  E(net)$outflow_z <- NA
  E(net)$weight    <- NA

  # Prepare progress bar
  if (progressbar) {
    n_max <- length(flowq)
    pb <- txtProgressBar(max = n_max, style = 3)
    progress <- function(n){setTxtProgressBar(pb, n)}
  }

  # Loop that runs until all nodes are done
  while (length(flowq) > 0) {
    if (progressbar) {
      progress(n_max - length(flowq))
    }

    # In case there are rechecks to do...
    if (length(rechecks) > 0) {

      # Sort them by their height
      rechecks <- rechecks[order(pts$z[rechecks])]
      rechecks <- unique(rechecks)
      thisnode <- rechecks[1]
      rechecks <- rechecks[-1]

    # In case there are no rechecks to do...
    } else {
      thisnode <- flowq[1]
      flowq <- flowq[-1]
      if (length(flowq) > 0) {
        qtop <- flowq[1]
      }
      vetos <- c()
      visited <- c()
      if (!thisnode %in% names(V(net))) {
        net <- net + vertices(thisnode)
      }
    }
    if (thisnode %in% visited) {
      for (pred in .predecessor(net, thisnode)) {
        if (pred %in% visited & pts$z[as.numeric(pred)] > pts$z[as.numeric(thisnode)]) {
          vetos <- c(vetos, pred)
        }
      }
    }
    oldouts  <- .successor(net, thisnode)
    outf     <- .outflow(thisnode, net, deltri, veto = vetos, pts = pts)
    nextnode <- outf$nbrmax
    p_out    <- outf$outflow
    m        <- outf$m
    for (oldout in oldouts) {
      net <- net - E(net, as.character(c(thisnode, oldout)))
      if (pts$z[as.numeric(oldout)] > pts$z[qtop]) {
        rechecks <- c(rechecks, oldout)
      }
    }
    visited <- c(visited, thisnode)

    # If there is no next node
    if (is.na(nextnode)) {
      if (thisnode %in% hull & length(.predecessor(net, thisnode)) > 0) {
        next
      } else {
        dtmax <- -Inf
        vf <- c(0, 0, 0)
        nmax <- NA
        for (n in .neighbor(deltri, thisnode)) {
          if (pts$z[n] < pts$z[as.numeric(thisnode)]) {
            flow <- .flowrate(pts, thisnode, n)
            dt   <- flow$dt
            vout <- flow$vf
            if (dt > dtmax) {
              dtmax <- dt
              vf <- vout
              nmax <- n
            }
          }
        }
        if (!is.na(nmax)) {
          outfl <- vf * m
          net <- net + edge(as.character(thisnode), as.character(nmax)
            , outflow_x = outfl[1]
            , outflow_y = outfl[2]
            , outflow_z = outfl[3]
            , weight    = m
          )
        }
      }

    # If there is a next node
    } else {
      if (!nextnode %in% names(V(net))) {
        net <- net + vertices(as.character(nextnode))
      }
      net <- net + edge(as.character(thisnode), as.character(nextnode)
        , outflow_x = p_out[1]
        , outflow_y = p_out[2]
        , outflow_z = p_out[3]
        , weight    = m
      )
      if (pts$z[nextnode] > pts$z[qtop]) {
        rechecks <- c(rechecks, nextnode)
      }
    }
    # print(counter)
  }

  # Convert the network to a spatial polygons dataframe
  net <- .net2spdf(net, lay = pts, xmn = xmn, ymn = ymn, scl = scl)

  # If a raster is deisred, rasterize the river network
  if (raster) {
    r <- rast(net, nrow = nrow, ncol = ncol)
    net <- rasterize(net, r, field = 1, background = 0)
  }

  # Return the network
  return(net)

}

################################################################################
#### Level 2 Functions
################################################################################
# Function to generate the height of a node
.height <- function(x, y, deg) {
  f <- 0
  for (i in 1:deg) {
    for (j in 1:deg) {
      coef <- runif(n = 1, min = -1, max = +1)
      f <- f + coef * x ** i * y ** j / (i + j + 1)
    }
  }
  return(f)
}

# Function to sample random points in space with a random height
.samplePoints <- function(n, deg = 1) {
  x <- runif(n, min = 0, max = 1)
  y <- runif(n, min = 0, max = 1)
  z <- .height(x, y, deg)
  pts <- data.frame(x, y, z)
  pts <- pts[order(pts$z, decreasing = T), ]
  pts$ID <- 1:nrow(pts)
  rownames(pts) <- NULL
  return(pts)
}

# Function that returns the neighbors the i-th node based on a delaunay
# triangulation
.neighbor <- function(deltri, i = NULL) {
  neigh <- subset(deltri$delsgs, ind1 == i | ind2 == i)
  neigh <- c(neigh$ind1, neigh$ind2)
  neigh <- unique(neigh[neigh != i])
  return(neigh)
}

# Functions to find the predecessors and successors of a node
.successor <- function(net, id) {
  as_ids(adjacent_vertices(net, v = as.character(id), mode = c("out"))[[1]])
}
.predecessor <- function(net, id) {
  as_ids(adjacent_vertices(net, v = as.character(id), mode = c("in"))[[1]])
}


# Function to find the flowrate between two nodes
.flowrate <- function(pts, p1, p2, vel = c(0, 0, 0)) {
  dx <- pts[p2, ]$x - pts[p1, ]$x
  dy <- pts[p2, ]$y - pts[p1, ]$y
  dz <- pts[p2, ]$z - pts[p1, ]$z
  dist <- sqrt(dx ** 2 + dy ** 2 + dz ** 2)
  vi <- (dx * vel[1] + dy * vel[2] + dz * vel[3]) / dist
  vfsq <- vi ** 2 - 2 * dz
  vf <- c(0, 0, 0)
  dt <- -Inf
  if (vfsq >= 0) {
    dt <- (vi + sqrt(vfsq)) / (2 * dist)
    vhat <- sqrt(vfsq) / dist
    vf <- c(vhat * dx, vhat * dy, vhat * dz)
  }
  return(list(dt = dt, vf = vf))
}

# Function that finds the outflow of a node
.outflow <- function(pt, net, deltri, veto, vetouphill = F, pts = pts) {
  px <- py <- pz <- 0
  m <- 1
  for (pred in .predecessor(net, pt)) {
    id <- get.edge.ids(net, as.character(c(pred, pt)))
    px <- px + E(net)$outflow_x[id]
    py <- py + E(net)$outflow_y[id]
    pz <- pz + E(net)$outflow_z[id]
    m <- m + E(net)$weight[id]
    id2 <- get.edge.ids(net, as.character(c(pt, pred)))
    if (pred %in% .successor(net, pt)) {
      if (E(net)$weight[id] > E(net)$weight[id2]) {
        px <- px - E(net)$outflow_x[id2]
        py <- py - E(net)$outflow_y[id2]
        pz <- pz - E(net)$outflow_z[id2]
        m <- m - E(net)$weight[id2]
      }
    }
  }
  vi     <- c(px / m, py / m, pz / m)
  dtmax  <- 0
  vf     <- c(0, 0, 0)
  nbrmax <- NA
  for (nbr in .neighbor(deltri, pt)) {
    if (length(veto) > 0 & nbr %in% veto) {
      next
    }
    if (vetouphill & pts$z[nbr] > pts$z[as.numeric(pt)]) {
      next
    }
    flow <- .flowrate(pts, as.numeric(pt), as.numeric(nbr), vi)
    dt <- flow$dt
    vn <- flow$vf
    if (dt > dtmax) {
      dtmax  <- dt
      vf     <- vn
      nbrmax <- nbr
    }
  }
  outflow <- vf * m
  return(list(nbrmax = nbrmax, outflow = outflow, m = m))
}

# Create Spatial Lines DataFrame
.net2spdf <- function(net, lay, xmn = 0, ymn = 0, scl = 1) {

  # Generate edgelist from network
  edges <- as_edgelist(net)
  edges <- as.data.frame(edges, stringsAsFactors = F)
  names(edges) <- c("From", "To")
  edges$Weight <- E(net)$weight

  # Replace node ids with actual node coordinates
  edges$x1 <- lay$x[match(edges$From, lay$ID)]
  edges$y1 <- lay$y[match(edges$From, lay$ID)]
  edges$x2 <- lay$x[match(edges$To, lay$ID)]
  edges$y2 <- lay$y[match(edges$To, lay$ID)]

  # Convert to spatial lines
  lines <- lapply(1:nrow(edges), function(x) {
    l <- as.lines(vect(rbind(
        c(edges$x1[x], edges$y1[x])
      , c(edges$x2[x], edges$y2[x])
    )))
    return(l)
  })
  lines <- do.call(rbind, lines)

  # Buffer lines according to their computed width
  width <- sqrt(edges$Weight) / sqrt(nrow(edges)) * 8 / 500
  lines_b <- buffer(lines, width = width)
  lines_b <- aggregate(lines_b, dissolve = T)

  # If required, stretch
  if (xmn != 0 | ymn != 0 | scl != 1) {
    lines_b <- rescale(lines_b, fx = scl)
    lines_b <- shift(lines_b, dx = xmn - xmin(lines_b), dy = ymn - ymin(lines_b))
  }

  # Return the final lines
  return(lines_b)
}
