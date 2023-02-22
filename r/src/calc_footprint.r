#' calc_footprint generates upstream influence footprint
#' @author Ben Fasoli
#'
#' Aggregates the upstream particle trajectories into a time integrated
#' footprint, expanding particle influence using variable 2d gaussian kernels
#' with bandwidths proportional to the mean pairwise distance between all
#' particles at each time step. Requires compiled permute.so to build the
#' gaussian kernels with fortran.
#'
#' @param p data frame containing particle trajectories, typically obtained
#'   from PARTICLE.rds but can be derived from PARTICLE.dat with column names
#'   equivalent to \code{varsiwant}. Must contain colums specifying long, lati,
#'   indx, foot, and time.
#' @param output filename argument. Must end with .nc (for ncdf output,
#'   preferred), .rds (for serialized R data output), .csv (for comma separated
#'   value table), or NULL to return the footprint object for programatic use.
#'   .nc files are saved in the CF-1.4 (Climate and Forcast Metadata) convention
#'   for native use with raster::brick() and raster::raster(). rds files do not
#'   require any additional libraries and have better compression
#' @param r_run_time receptor run time as a POSIXct object. Can be NULL
#'   resulting in NULL timestamp outputs
#' @param projection proj4 string defining the map projection of the footprint
#'   netCDF output
#' @param time_integrate logical indicating whether to integrate footprint over
#'   time or retain discrete time steps
#' @param smooth_factor factor by which to linearly scale footprint smoothing;
#'   0 to disable, defaults to 1
#' @param xmn sets grid start longitude
#' @param xmx sets grid end longitude
#' @param xres resolution for longitude grid
#' @param ymn sets grid start latitude
#' @param ymx sets grid end latitude
#' @param yres resolution for latitude grid
#'
#' @import dplyr, proj4, raster
#' @export

# Gaussian kernel weighting calculation
make_gauss_kernel <- function (rs, sigma, projection) {
  # Modified from raster:::.Gauss.weight()
  if (sigma == 0) {
    return(array(1, c(1, 1)))
  }
  require(raster)
  d <- 3 * sigma
  nx <- 1 + 2 * floor(d/rs[1])
  ny <- 1 + 2 * floor(d/rs[2])
  m <- matrix(ncol = nx, nrow = ny)
  xr <- (nx * rs[1])/2
  yr <- (ny * rs[2])/2
  r <- raster(m, xmn = -xr[1], xmx = xr[1], ymn = -yr[1], ymx = yr[1],
              crs = projection)
  p <- xyFromCell(r, 1:ncell(r))^2
  m <- 1/(2 * pi * sigma^2) * exp(-(p[, 1] + p[, 2])/(2 * sigma^2))
  m <- matrix(m, ncol = nx, nrow = ny, byrow = TRUE)
  w <- m/sum(m)
  w[is.na(w)] <- 1
  w
}

# Wrap longitudes into -180:180 range
wrap_longitude_meridian <- function(x) {
  (x %% 360 + 540) %% 360 - 180
}

# Wrap longitudes into 0:360 range
wrap_longitude_antimeridian <- function(x) {
  x <- wrap_longitude_meridian(x)
  ifelse(x < 0, x + 360, x)
}

calc_footprint <- function(p, output = NULL, r_run_time,
                           projection = '+proj=longlat',
                           smooth_factor = 1, time_integrate = F,
                           xmn, xmx, xres, ymn, ymx, yres = xres) {

  message("<<<  Inside calc_footprint.r ")
  begin_time = Sys.time();

  require(dplyr)
  require(raster)
  
  np <- length(unique(p$indx))
  time_sign <- sign(median(p$time))
  is_longlat <- grepl('+proj=longlat', projection, fixed = T)

  message("<<<  Loaded R libraries at time: ", Sys.time() - begin_time)
  
  # Determine longitude wrapping behavior for grid extents containing anti
  # meridian, including partial wraps (e.g. 20deg from 170:-170) and global
  # coverage (e.g. 360deg from -180:180)
  if (is_longlat) {
    xdist <- ((180 - xmn) - (-180 - xmx)) %% 360
    if (xdist == 0) {
      xdist <- 360
      xmn <- -180
      xmx <- 180
    } else if ((xmx < xmn) || (xmx > 180)){
      p$long <- wrap_longitude_antimeridian(p$long)
      xmn <- wrap_longitude_antimeridian(xmn)
      xmx <- wrap_longitude_antimeridian(xmx)
    }
  }

  message("<<<  Determine longitude wrapping behavior at time ", Sys.time() - begin_time)
  
  # Interpolate particle locations during first 100 minutes of simulation if
  # median distance traveled per time step is larger than grid resolution
  distances <- p %>%
    dplyr::filter(abs(time) < 100) %>%
    group_by(indx) %>%
    summarize(dx = median(abs(diff(long))),
              dy = median(abs(diff(lati)))) %>%
    ungroup()
  
  should_interpolate <- (median(distances$dx, na.rm = T) > xres) || 
    (median(distances$dy, na.rm = T) > yres)
  if (should_interpolate) {
    times <- c(seq(0, 10, by = 0.1),
               seq(10.2, 20, by = 0.2),
               seq(20.5, 100, by = 0.5)) * time_sign
    
    # Preserve total field prior to split-interpolating particle positions
    aptime <- abs(p$time)
    foot_0_10_sum <- sum(p$foot[aptime <= 10], na.rm = T)
    foot_10_20_sum <- sum(p$foot[aptime > 10 & aptime <= 20], na.rm = T)
    foot_20_100_sum <- sum(p$foot[aptime > 20 & aptime <= 100], na.rm = T)
    
    # Split particle influence along linear trajectory to sub-minute timescales
    p <- p %>%
      full_join(expand.grid(time = times,
                            indx = unique(p$indx)),
                by = c('indx', 'time')) %>%
      arrange(indx, -time) %>%
      group_by(indx) %>%
      mutate(long = na_interp(long, x = time),
             lati = na_interp(lati, x = time),
             foot = na_interp(foot, x = time)) %>%
      ungroup() %>%
      na.omit() %>%
      mutate(time = round(time, 1))
    
    # Scale interpolated values to retain total field
    aptime <- abs(p$time)
    mi <- aptime <= 10
    p$foot[mi] <- p$foot[mi] / (sum(p$foot[mi], na.rm = T) / foot_0_10_sum)
    mi <- aptime > 10 & aptime <= 20
    p$foot[mi] <- p$foot[mi] / (sum(p$foot[mi], na.rm = T) / foot_10_20_sum)
    mi <- aptime > 20 & aptime <= 100
    p$foot[mi] <- p$foot[mi] / (sum(p$foot[mi], na.rm = T) / foot_20_100_sum)
  }

  message("<<<  Interpolate particle location at time: ", Sys.time() - begin_time)
  
  # Preserve time relative to individual particle release as rtime
  p <- p %>%
    group_by(indx) %>%
    mutate(rtime = time - (time_sign) * min(abs(time))) %>%
    ungroup()

  message("<<<  Preserve time relative to individual particle release as rtime at time: ", Sys.time() - begin_time)
  
  # Translate x, y coordinates into desired map projection
  if (!is_longlat) {
    require(proj4)
    p[, c('long', 'lati')] <- project(p[, c('long', 'lati')], projection)
    grid_lims <- project(list(c(xmn, xmx), c(ymn, ymx)), projection)
    xmn <- min(grid_lims$x)
    xmx <- max(grid_lims$x)
    ymn <- min(grid_lims$y)
    ymx <- max(grid_lims$y)
  }

  message("<<<  Translate x, y coordinates into desired map projection at time: ", Sys.time() - begin_time)
  
  # Set footprint grid breaks using lower left corner of each cell
  glong <- head(seq(xmn, xmx, by = xres), -1)
  glati <- head(seq(ymn, ymx, by = yres), -1)

  message("<<<  Set footprint grid breaks using lower left corner of each cell at time: ", Sys.time() - begin_time)
  
  # Gaussian kernel bandwidth scaling by summed variances, elapsed time, and
  # average latitude of the ensemble
  kernel <- p %>%
    group_by(rtime) %>%
    dplyr::summarize(varsum = var(long) + var(lati),
                     lati = mean(lati)) %>%
    ungroup() %>%
    na.omit()
  
  di <- kernel$varsum^(1/4)
  ti <- abs(kernel$rtime/1440)^(1/2)
  grid_conv <- if (is_longlat) cos(kernel$lati * pi/180) else 1
  w <- smooth_factor * 0.06 * di * ti / grid_conv
  
  message("<<<  Gaussian kernel bandwidth scaling at time: ", Sys.time() - begin_time)

  # Determine maximum kernel size -> 0.4 sec
  xyres <- c(xres, yres)
  max_k <- make_gauss_kernel(xyres, max(w), projection)

  message("<<<  Determined maximum kernel size at time: ", Sys.time() - begin_time)
  
  # Expand grid extent using maximum kernel size
  xbuf <- ncol(max_k)
  xbufh <- (xbuf - 1) / 2
  ybuf <- nrow(max_k)
  ybufh <- (ybuf - 1) / 2
  
  glong_buf <- seq(xmn - (xbuf*xres), xmx + ((xbuf - 1)*xres), by = xres)
  glati_buf <- seq(ymn - (ybuf*yres), ymx + ((ybuf - 1)*yres), by = yres)

  message("<<<  Expand grid extent using maximum kernel size at time: ", Sys.time() - begin_time)
  
  # Remove zero influence particles and positions outside of domain
  p <- p %>%
    dplyr::filter(foot > 0,
                  long >= (xmn - xbufh*xres), long < (xmx + xbufh*xres),
                  lati >= (ymn - ybufh*yres), lati < (ymx + ybufh*yres))
  if (nrow(p) == 0) return(NULL)

  message("<<<  Remove zero influence particles and positions outside of domain at time: ", Sys.time() - begin_time)
  
  # Pre grid particle locations
  p <- p %>%
    transmute(loi = as.integer(findInterval(long, glong_buf)),
              lai = as.integer(findInterval(lati, glati_buf)),
              foot = foot,
              time = time,
              rtime) %>%
    group_by(loi, lai, time, rtime) %>%
    dplyr::summarize(foot = sum(foot, na.rm = T)) %>%
    ungroup()
  
  message("<<<  Pre grid particle locations at time: ", Sys.time() - begin_time)

  # Dimensions in accordance with CF convention (x, y, t)
  nx <- length(glati_buf)
  ny <- length(glong_buf)
  grd <- matrix(0, nrow = ny, ncol = nx)

  message("<<<  Dimensions in accordance with CF convention (x, y, t) at time: ", Sys.time() - begin_time)

  
  # Split particle data by footprint layer
  interval <- 3600
  interval_mins <- interval / 60
  p$layer <- if (time_integrate) 0 else floor(p$time / interval_mins)
  
  layers <- sort(unique(p$layer))
  nlayers <- length(layers)

  message("<<<  Split particle data by footprint layer at time: ", Sys.time() - begin_time)
  
  perf0 <- 0
  perf1 <- 0
  perf2 <- 0
  perf3 <- 0
  perf4 <- 0
  perf5 <- 0
  perf6 <- 0
  perf7 <- 0
  
  # Allocate and fill footprint output array
  foot <- array(grd, dim = c(dim(grd), nlayers))
  message("nlayers is ", nlayers, ", grid is ", nx, " x ", ny)
  for (i in 1:nlayers) {
    .C("create_footprint", nrow=as.integer(ny), ncol=as.integer(nx))
    perf0 <- perf0 - as.numeric(Sys.time())
    layer_subset <- dplyr::filter(p, layer == layers[i])
    perf0 <- perf0 + as.numeric(Sys.time())
    
      # TODO: could we maybe loop once and slot all particles into one list per timestep?
      # fast way to do this:  https://stackoverflow.com/a/29870770/3304288
      #list_ = {
      #      a <- list(0)
      #      for(i in 1:n) {a <- list(a, list(i))}
      #  },
    perf1 <- perf1 - as.numeric(Sys.time())
    rtimes <- unique(layer_subset$rtime)
    perf1 <- perf1 + as.numeric(Sys.time())
      # message("length(rtimes) is ", length(rtimes))
      # message("")
      # message("length(layer_subset$rtime) is ", length(layer_subset$rtime))
      # message("length(layer_subset) is ", length(layer_subset))
      # message("typeof(layer_subset) is ", typeof(layer_subset))
      # message("layer_subset[1] is ", layer_subset[1][1:10])
      # message("layer_subset[2] is ", layer_subset[2][1:10])
      
      # # q: how do you slice a vector
     
    by_time <- split(layer_subset, f = layer_subset$rtime)

      # message("by_time ", by_time)
      # message("length(by_time) ", length(by_time))
      # message("by_time[1] ", by_time[1][1:10])
      # message("by_time[0] ", by_time[0][1:10])


    for (j in 1:length(rtimes)) {
      perf2 <- perf2 - as.numeric(Sys.time())
      step <- bind_rows(by_time[j])
      #step <- as.data.frame(by_time[j])
      #colnames(step) <- names(layer_subset)
      #step <- dplyr::filter(layer_subset, rtime == rtimes[j]) #len == 6
      
      #print(step)
      #print(typeof(step))
      perf2 <- perf2 + as.numeric(Sys.time()) # 0.6 sec -> 0.06 sec
      
      # Create dispersion kernel based using nearest-in-time kernel bandwidth
      perf3 <- perf3 - as.numeric(Sys.time())
      step_w <- w[find_neighbor(rtimes[j], kernel$rtime)]
      perf3 <- perf3 + as.numeric(Sys.time())

      perf4 <- perf4 - as.numeric(Sys.time()) # 1.3 sec -> 0 sec
      #message("kernel ", j, " ", xyres, " ", step_w, " ", projection)
      #k <- make_gauss_kernel(xyres, step_w, projection)
      perf4 <- perf4 + as.numeric(Sys.time())
      
      perf5 <- perf5 - as.numeric(Sys.time())
      # Array dimensions
      len <- nrow(step)
      # print(len)
      #nkx <- ncol(k)
      #nky <- nrow(k)
      #message("kernel dims ", nkx, " x ", nky)
      perf5 <- perf5 + as.numeric(Sys.time())
      
      perf6 <- perf6 - as.numeric(Sys.time()) # 0.6 sec -> 0.043 sec
      # Call permute fortran subroutine to build and aggregate kernels
      # out <- .Fortran('permute', ans = grd, nax = nx, nay = ny, k = k, 
      #                 nkx = nkx, nky = nky, len = len, lai = step$lai, 
      #                 loi = step$loi, foot = step$foot)


      d <- 3 * step_w # 3 * sigma
      stopifnot(xyres[1] == xyres[2])
      nkx <- as.integer(1 + 2 * as.integer(d/xyres[1]))
      nky <- as.integer(1 + 2 * as.integer(d/xyres[2]))
      stopifnot(nkx == nky)
      stopifnot(is.integer(nkx))
      # message("nkx ", nkx)

      # Compute sigma in units of matrix coords
      sigma_elements = step_w/xyres[1]


      .C('permute', 
          sigma=sigma_elements, nkx=nkx, nky=nky, # Kernel and its dimensions
          len = len, # Number of locations to add kernel
          lai = step$lai, loi = step$loi, # xs (lat!) and ys (long!) to add kernel
          foot = step$foot) # weights to add kernel

      perf6 <- perf6 + as.numeric(Sys.time())

      perf7 <- perf7 - as.numeric(Sys.time()) # 3.2 sec -> 0.002 sec
      perf7 <- perf7 + as.numeric(Sys.time())
    }
    ans = .Call('read_footprint')
    foot[ , , i] <- matrix(ans, nrow=ny)
  }

  message("<<< perf0: ", perf0)
  message("<<< perf1: ", perf1)
  message("<<< perf2: ", perf2)
  message("<<< perf3: ", perf3)
  message("<<< perf4: ", perf4)
  message("<<< perf5: ", perf5)
  message("<<< perf6: ", perf6)
  message("<<< perf7: ", perf7)
  message("<<<  Allocate and fill footprint output array at time: ", Sys.time() - begin_time)
  
  # Remove spatial buffer around domain used in kernel aggregation
  size <- dim(foot)
  foot <- foot[(xbuf+1):(size[1]-xbuf), (ybuf+1):(size[2]-ybuf), ] / np

  message("<<<  Remove spatial buffer around domain used in kernel aggregation at time: ", Sys.time() - begin_time)
  
  # Determine time to use in output files
  if (time_integrate) {
    time_out <- as.numeric(r_run_time) 
  } else {
    time_out <- as.numeric(r_run_time + layers * interval)
  }

  message("<<<  Determine time to use in output files at time: ", Sys.time() - begin_time)
  
  # Set footprint metadata and write to file
  write_footprint(foot, output = output, glong = glong, glati = glati,
                  projection = projection, xres = xres, yres = yres,
                  time_out = time_out)
  message("<<< Wrote to: ", output)
  message("<<<  Set footprint metadata and write to file at time: ", Sys.time() - begin_time)

}
