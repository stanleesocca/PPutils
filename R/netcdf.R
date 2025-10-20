var_descriptions <- list(chl = list(varname = "chl", 
                                    unit = "mg/m3", 
                                    longname = "Surface chlorophyll a", 
                                    other_name = c("chl_xyv"),
                                    attributes = NULL), 
                         ppPel = list(varname = "ppPel", 
                                      unit = "mg/m3/h", 
                                      longname = "Integrated Pelagic primary production", 
                                      other_name = c(""),
                                      attributes = NULL),
                         ppBen = list(varname = "ppBen", 
                                      unit = "mg/m3/h", 
                                      longname = "Integrated Benthic primary production", 
                                      other_name = c(""),
                                      attributes = NULL),
                         ndvi = list(varname = "ndvi", 
                                      unit = "-", 
                                      longname = "Normalized difference vegetation index", 
                                     other_name = c("ndvi_xyv"),
                                      attributes = NULL), 
                         WHeight = list(varname = "WHeight", 
                                     unit = "m", 
                                     longname = "Water height", 
                                     other_name = c("WHeight", "Wheight"),
                                     attributes = NULL), 
                         Irrad = list(varname = "Irrad", 
                                        unit = "umol photon/m2/s", 
                                        longname = "Photosynthetically Active Radiation", 
                                        other_name = c("PAR"),
                                        attributes = NULL), 
                         WKd = list(varname = "WKd", 
                                      unit = "/m", 
                                      longname = "Light extinction coefficient", 
                                      other_name = c("Kd"),
                                      attributes = NULL), 
                         WEopt = list(varname = "WEopt", 
                                    unit = "uEinst/m2/s", 
                                    longname = "Optimal light intensity", 
                                    other_name = c("Iopt"),
                                    attributes = NULL), 
                         WPs = list(varname = "WPs", 
                                      unit = "mg C/(mg chla)/h", 
                                      longname = "Maximum photosynthesis rate", 
                                      other_name = c("Pmax", "WPs"),
                                      attributes = NULL), 
                         WAlpha = list(varname = "WAlpha", 
                                    unit = "mg C/(mg chl)/h/(uE/m2/s)", 
                                    longname = "light-limited initial slope of the photosynthesis-irradiance (P-I) curve", 
                                    other_name = c("WAlpha"),
                                    attributes = NULL), 
                         Rad = list(varname = "Rad", 
                                       unit = "uE/m2/s", 
                                       longname = "Radiation at the bottom", 
                                       other_name = c("Rad"),
                                       attributes = NULL), 
                         BEopt = list(varname = "BEopt", 
                                    unit = "uEinst/m2/s", 
                                    longname = "Optimal light intensity in the sediment", 
                                    other_name = c("BEopt"),
                                    attributes = NULL), 
                         BPs = list(varname = "BPs", 
                                      unit = "mg C/(mg chla)/h", 
                                      longname = "Maximum benthic photosynthesis rate", 
                                      other_name = c("BPs"),
                                      attributes = NULL), 
                         BAlpha = list(varname = "BAlpha", 
                                    unit = "mg C/(mg chl)/h/(uE/m2/s)", 
                                    longname = "light-limited initial slope of the benthic photosynthesis-irradiance (P-I) curve", 
                                    other_name = c("BPs"),
                                    attributes = NULL)
                         )

rda_to_netcdf <- function(rda_path, start_depth = 1, filename){
  
  if(is.list(rda_path)){
    rda <- rda_path
    file_type <- "list"
  } else if (is.character(rda_path) & file.exists(rda_path)){
    res_env <- new.env()
    load(rda_path, envir = res_env)
    rda <- res_env$res
    file_type <- "rda"
  } else {
    stop(rda_path, "should be either a saved .rda from previous computation\n
         or a list containing the result of a computation")
  }
  

  
  # required name by convention 
  req_names <- c("output_xy", "output_t")
  
  prim_vars <- c("chl_xyv", "ndvi_xyv", "ppPel", "WHeight", "WTemperature", 
                 "Irrad", "WPs", "WAlpha", "WKd", 
                 "ppBen", "BEopt", "BPs", "BAlpha", "Rad")
  aux_vars <- c("bbox", "input_t", "Bat_xyv")
  
  
  names_accessor <- purrr::attr_getter("names")
  
  # all names in pelagic 
  all_pnms <- purrr::pluck(rda, "pel", names_accessor)
  all_pnms <- all_pnms[grep("out", all_pnms, invert = TRUE)]
  names(all_pnms) <- all_pnms

  all_bnms <- purrr::pluck(rda, "ben", names_accessor)
  names(all_bnms) <- all_bnms
  
  out_nms <- purrr::pluck(rda, start_depth, "out", names_accessor)
  pel_nms <- purrr::pluck(rda, start_depth, "pel", names_accessor)
  names(out_nms) <- out_nms
  
  
  if(!any(out_nms %in% req_names)){
    stop("Missing neccessary data 'output_xy' and 'output_t' to construct the file.")
  }
  
  output_xy <- purrr::pluck(rda, start_depth, "out", out_nms[['output_xy']])
  output_t  <- purrr::pluck(rda, start_depth, "out", out_nms[['output_t']])
  bat       <- purrr::pluck(rda, start_depth, "out", out_nms[['Bat_xyv']])
  
  mv <- NA #missing value to use
  
  # lon    <- ncdim_def("longitude", "degrees_east", output_xy[, 1])
  # lat    <- ncdim_def("latitude", "degrees_north", output_xy[, 2])
  lonlat <- ncdim_def("ij", "lon-lat-indices", output_xy[, 2])
  timedim <- ncdim_def("timestamp", "timestamp", 1)
  onedim <- ncdim_def("1D", "1D", 1)
  twodim <- ncdim_def("2D", "", 1:2)
  time   <- ncdim_def("time","day", as.double(output_t), unlim=TRUE)
  
  vardepth    <- ncvar_def("Depth", units = "m", list(lonlat, onedim), mv )
  lon     <- ncvar_def("longitude", units = "degrees_east", list(lonlat, onedim), mv)
  lat     <- ncvar_def("latitude", units = "degrees_west", list(lonlat, onedim), mv)
  XY     <- ncvar_def("lonlat", units = "degree", dim = list(lonlat, twodim), prec = typeof(output_xy[, 1]), missval =  mv)
  
  ncnew <- nc_create(filename, list(lon, lat, XY, vardepth), force_v4=TRUE)
  
  ncvar_put(ncnew, lon, output_xy[, 1], start = c(1, 1), count = c(-1, -1))
  ncvar_put(ncnew, lat, output_xy[, 2], start = c(1, 1), count = c(-1, -1))
  ncvar_put(ncnew, XY, as.matrix(output_xy[, c(1, 2)], nrow = nrow(output_xy), ncol = 2))
  ncvar_put(ncnew, vardepth, bat[, 3], start = c(1, 1), count = c(-1, -1))

  nc_close(ncnew)


  
  out_nms <- out_nms[out_nms %in% prim_vars]
  all_pnms <- all_pnms[all_pnms %in% prim_vars]
  all_bnms <- all_bnms[all_bnms %in% prim_vars]
  all_nms <- c(all_pnms, all_bnms)
  
  nms_in_var_desc <- names(var_descriptions)


  # THINK !!!!! 
  
  inms <- which(nms_in_var_desc %in% out_nms)


  if(length(inms) == 0){
    nms <- sapply(var_descriptions, function(x) x[["other_name"]])
    nms <- nms[which(nms != "")]
    inms <- which(nms %in% out_nms) 
    if(length(inms) > 0) nms <- nms[inms]
    # if(length(inms) > 0) nms <- sapply(inms, function(x) nms[x])

    out_nms <- out_nms[which(nms != "")]
    nms_in_var_desc <- names(nms)
  } else {
    nms <- var_descriptions[[inms]][["varname"]]
    out_nms <- out_nms[[nms]]
    nms_in_var_desc <- names(out_nms)
  }

    
  
  ncid <- nc_open(filename, write = TRUE)
  
  for(i in 1:length(nms_in_var_desc)){
    var_select <- nms_in_var_desc[[i]]
    var_desc <- var_descriptions[[var_select]]
  
    var_rs <- purrr::pluck(rda, start_depth, "out", out_nms[i])[, -c(1, 2)]
    
    vardef    <- ncvar_def(name  = var_desc[["varname"]],
                           units = var_desc[["unit"]],
                           dim = list(lonlat, timedim, time),
                           longname = var_desc[["longname"]],
                           prec = typeof(var_rs[, 1]), 
                           missval =  mv)

    ncid <- ncvar_add(ncid, vardef)

    for(it in 1:length(output_t)){
      ncvar_put(ncid, vardef, var_rs[, it], start = c(1, 1, it), count = c(-1, -1, 1))
    }
      
  }
  
  ## Other pelagic variables 
  nms_in_var_desc <- names(var_descriptions)
  inms <- which(nms_in_var_desc %in% all_pnms)


  # return(inms)
  if(length(inms) == 0){
    nms <- sapply(var_descriptions, function(x) x[["other_name"]])
    nms <- nms[which(nms != "")]
    all_pnms <- all_pnms[which(nms != "")]
    nms_in_var_desc <- names(nms)
  } else {
    nms <- sapply(inms, function(x){
      nms <- var_descriptions[[x]][["varname"]]
    })
    # nms <- var_descriptions[[inms]][["varname"]]
    all_pnms <- nms_in_var_desc <- all_pnms[nms]
  }

  
  for(i in 1:length(nms_in_var_desc)){
    var_select <- nms_in_var_desc[[i]]
    var_desc <- var_descriptions[[var_select]]
    
    cat(var_desc[["varname"]], "\n")
    
    vardef    <- ncvar_def(name  = var_desc[["varname"]],
                           units = var_desc[["unit"]],
                           dim = list(lonlat, timedim, time),
                           longname = var_desc[["longname"]],
                           missval =  mv)
    
    ncid <- ncvar_add(ncid, vardef)
    

    for(it in 1:length(output_t)){
      if(var_select %in% c("ppPel", "ppBen", "WEopt", "WAlpha", "WPs")){
        if(var_select %in% c("WEopt", "WAlpha", "WPs","BEopt", "BPs", "BAlpha")) {
          # variables with time as first column 
          var_rs <- t(purrr::pluck(rda, "pel", all_pnms[i])[, -1])
        } else {
          # variables without time as first column 
          var_rs <- t(purrr::pluck(rda, "pel", all_pnms[i]))
        }
      } else {
        var_rs <- purrr::pluck(rda, "pel", all_pnms[i])[, -c(1, 2)]
      }

      ncvar_put(ncid, vardef, var_rs[, it], start = c(1, 1, it), count = c(-1, -1, 1))
    }
    
  }

  
  ## Other benthic variables 
  nms_in_var_desc <- names(var_descriptions)
  inms <- which(nms_in_var_desc %in% all_bnms)
  
  
  # return(inms)
  if(length(inms) == 0){
    nms <- sapply(var_descriptions, function(x) x[["other_name"]])
    nms <- nms[which(nms != "")]
    all_bnms <- all_bnms[which(nms != "")]
    nms_in_var_desc <- names(nms)
  } else {
    nms <- sapply(inms, function(x){
      nms <- var_descriptions[[x]][["varname"]]
    })
    # nms <- var_descriptions[[inms]][["varname"]]
    all_bnms <- nms_in_var_desc <- all_bnms[nms]
  }
  
  
  for(i in 1:length(nms_in_var_desc)){
    var_select <- nms_in_var_desc[[i]]
    var_desc <- var_descriptions[[var_select]]
    
    cat(var_desc[["varname"]], "\n")
    
    vardef    <- ncvar_def(name  = var_desc[["varname"]],
                           units = var_desc[["unit"]],
                           dim = list(lonlat, timedim, time),
                           longname = var_desc[["longname"]],
                           missval =  mv)
    
    ncid <- ncvar_add(ncid, vardef)
    
    
    for(it in 1:length(output_t)){
      if(var_select %in% c("ppPel", "ppBen", "WEopt", "WAlpha", "WPs", 
                           "Rad", "BEopt", "BPs", "BAlpha")){
        if(var_select %in% c("WEopt", "WAlpha", "WPs","BEopt", "BPs", "BAlpha")) {
          # variables with time as first column 
          var_rs <- t(purrr::pluck(rda, "ben", all_bnms[i])[, -1])
        } else {
          # variables without time as first column 
          var_rs <- t(purrr::pluck(rda, "ben", all_bnms[i]))
        }

      } else {
        var_rs <- purrr::pluck(rda, "ben", all_bnms[i])[, -c(1, 2)]
      }
      
      ncvar_put(ncid, vardef, var_rs[, it], start = c(1, 1, it), count = c(-1, -1, 1))
    }
    
  }
  
  
  # # lastly gridded temperature
  vardef    <- ncvar_def(name  = "WTemp",
                         units = "celsius",
                         dim = list(lonlat, timedim, time),
                         longname = "Water temperature",
                         missval =  mv)
  
  ncid <- ncvar_add(ncid, vardef)
  for(it in 1:length(output_t)){
    var_rs <- purrr::pluck(rda, "pel", "WTemperature", "WTemp")[, -c(1, 2)]
    ncvar_put(ncid, vardef, var_rs[, it], start = c(1, 1, it), count = c(-1, -1, 1))
  }
  
  # add global attributes
  ncatt_put(ncid,0,"title", "Primary production results from dtRprimprod")
  ncatt_put(ncid,0,"institution", "NIOZ")
  ncatt_put(ncid,0,"source", "Remote sensing chlorophylla processed from sentinel-2, other data from Rijkswaterstraat and KNMI")
  history <- paste("Stanley Nmor", date(), sep=", ")
  ncatt_put(ncid,0,"history", history)
  
  # Minimum metadata based on Parinaz recommendation 
  ncatt_put(ncid, 0, "Responsible party", "EDS, NIOZ")
  ncatt_put(ncid, 0, "electronic mail address", "")
  ncatt_put(ncid, 0, "creation date", date())
  data_description <- "Data-driven estimate of remote sensed primary production and chlorophyll-a distribution. The dataset also contain gridded data used for intermediate computation"
  ncatt_put(ncid, 0, "Description", data_description)
  ncatt_put(ncid, 0, "Data creator", "Stanley Nmor, Karline Soetaert")
  ncatt_put(ncid, 0, "Data contact point", "")
  ncatt_put(ncid, 0, "Data publisher", "Stanley Nmor")
  
  if(file_type == "rda"){
    ncatt_put(ncid, 0, "Spatial coverage", res_env$res$pel$out$bbox)
    sp_res <- c(max(diff(res_env$res$pel$out$output_xy[, 1])), max(diff(res_env$res$pel$out$output_xy[, 2])))
    ncatt_put(ncid, 0, "Temporal coverage", range(res_env$res$pel$out$output_t))
  } else if(file_type == "list"){
    ncatt_put(ncid, 0, "Spatial coverage", rda$pel$out$bbox)
    sp_res <- c(max(diff(res$pel$out$output_xy[, 1])), max(diff(rda$pel$out$output_xy[, 2])))
    ncatt_put(ncid, 0, "Temporal coverage", range(rda$pel$out$output_t))
  }
  
  ncatt_put(ncid, 0, "Spatial resolution", paste0(sp_res, " deg"))
  ncatt_put(ncid, 0, "Spatial reference system", "re-gridded from a irregular curvillinear grid with WGS 64")
  ncatt_put(ncid, 0, "Temporal resolution", "Daily available snapshot")

  ncatt_put(ncid, 0, "Resource type", "Empirical model data output")
  ncatt_put(ncid, 0, "Keywords", "Pelagic and benthic primary production, Chlorophyll-a, Interpolated Water height")
  ncatt_put(ncid, 0, "License", "")
  ncatt_put(ncid, 0, "Unique Identifier", "")
  ncatt_put(ncid, 0, "Access rights", "Open access")
  ncatt_put(ncid, 0, "Distribution access URL", "")
  ncatt_put(ncid, 0, "Distribution format", "netcf4")
  ncatt_put(ncid, 0, "Distribution byte size", "")
  ncatt_put(ncid, 0, "Landing page", "")

    
  nc_close(ncid)
  cat("nc file written to :", filename)

}
