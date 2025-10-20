get_acolite_vars <- function(files, ...){

  vars <- NULL
  for(i in 1:length(files)){
    nc <- nc_open(files[i])
    vars <- c(vars, names(nc$var))
  }
  nc_close(nc)
  return(unique(vars))
}

get_datatake_date <- function(date_str){
  res <- str_extract(date_str,
                     "[0-9]{8}T[0-9]{6}")

  res_date <- paste0(str_sub(res, start = 1, end = 4), "-",
                     str_sub(res, start = 5, end = 6), "-",
                     str_sub(res, start = 7, end = 8), " ",
                     str_sub(res, start = 10, end = 11), ":",
                     str_sub(res, start = 12, end = 13), ":",
                     str_sub(res, start = 14, end = 15))

  res_date <- as.POSIXct(res_date, tz = "CET")
  res_date
}



read_acolite_output <- function(files,
                                varnames = NULL,
                                ...)
{

  xyv_chl <- input_t <- datatake_t <- NULL

  for (ff in files){
    RS.nc <- nc_open(ff)


    vars_names <- names(RS.nc$var)

    req_names <- c("chl_re_mishra", "lon", "lat", "rhow_560", "rhow_833",
                   "rhow_1614", "rhos_665", "rhos_665", "rhos_833")

    # if(!varnames %in% vars_names)
    #   stop("One of the variable requested not in the ncdf4 file")

    chl <- ncvar_get(RS.nc, "chl_re_mishra")
    lon <- ncvar_get(RS.nc, "lon")
    lat <- ncvar_get(RS.nc, "lat")

    dst <- ncatt_get(RS.nc, 0, "isodate")$value
    dd <- strsplit(dst, "T")[[1]][1]
    tt <- stringr::str_extract_all(strsplit(dst, "T")[[1]][2],
                                   pattern = "[0-9]*:[0-9]*:[0-9]{2}")
    time <- as.POSIXct(paste(dd, tt))

    datatake_date <- get_datatake_date(basename(dirname(ff)))

    # xyv <- rbind(xyv,
    #              na.omit(cbind(as.vector(lon),
    #                            as.vector(lat),
    #                            as.vector(chl)))
    # )
    # B3   <- ncvar_get(RS.nc, varid =  "rhow_560")
    # B8a  <- ncvar_get(RS.nc, varid =  "rhow_833")
    # B11  <- ncvar_get(RS.nc, varid =  "rhow_1614")

    B4   <- ncvar_get(RS.nc, varid =  "rhos_665")
    B8   <- ncvar_get(RS.nc, varid =  "rhos_833")

    NDVI <- (B8 - B4)/(B8 + B4)


    if (is.null(xyv_chl)){
      xyv_chl <- cbind(as.vector(lon),
                       as.vector(lat),
                       as.vector(chl))

      xyv_ndvi  <- cbind(as.vector(lon),
                         as.vector(lat),
                         as.vector(NDVI))


    } else {
      xyv_chl <- cbind(xyv_chl, as.vector(chl))
      xyv_ndvi <- cbind(xyv_ndvi, as.vector(NDVI))
    }
    input_t <- c(input_t, time)
    datatake_t <- c(datatake_t, datatake_date)
  }

  input_t <- as.POSIXct(input_t)
  datatake_t <- as.POSIXct(datatake_t)
  nc_close(RS.nc)

  return(list(xyv_chl = xyv_chl, xyv_ndvi = xyv_ndvi,
              input_t = input_t,
              datatake_t = datatake_t))
}


create_spatio_temporal <- function(files,
                                   bathymetric_file,
                                   RS_data,
                                   method = c("average", "interpolate"),
                                   use_datatake_time = TRUE,
                                   ...
){

  method <- match.arg(method)


  xyv_chl <- RS_data[["xyv_chl"]]
  xyv_ndvi <- RS_data[["xyv_ndvi"]]

  if(use_datatake_time) {
    input_t <-  RS_data[["datatake_t"]]
  } else {
    input_t <-  RS_data[["input_t"]]
  }


  lon_out <- c(min(bathymetric_file[["longitude"]]),
               max(bathymetric_file[["longitude"]]),
               range(diff(bathymetric_file[["longitude"]]))[1])
  lat_out <- c(min(bathymetric_file[["latitude"]]),
               max(bathymetric_file[["latitude"]]),
               range(diff(bathymetric_file[["latitude"]]))[1])

  min_x = min(xyv_chl[, 1]); max_x = max(xyv_chl[, 1])
  min_y = min(xyv_chl[, 2]); max_y = max(xyv_chl[, 2])

  bbox <- c(min_x, min_y, max_x, max_y)

  by <- 1
  ix <- seq(1, length(bathymetric_file[["longitude"]]), by = by)
  iy <- seq(1, length(bathymetric_file[["latitude"]]) , by = by)

  Bat_xyv   <- with (bathymetric_file,
                     data.frame(expand.grid(longitude = longitude[ix],
                                            latitude  = latitude [iy]),
                                depth = as.vector(depth[ix, iy])))


  chl_xyv <- expand.grid(longitude = bathymetric_file[["longitude"]][ix],
                         latitude  = bathymetric_file[["latitude"]][iy])

  ndvi_xyv <- expand.grid(longitude = bathymetric_file[["longitude"]][ix],
                          latitude  = bathymetric_file[["latitude"]][iy])


  if(method == "average"){
    for(i in 3:(ncol(xyv_chl))){
      intChl <- average_xy(xyv_chl[, c(1, 2, i)], output_x = unique(Bat_xyv[, 1]), output_y = unique(Bat_xyv[, 2]))[[3]]
      intChl[intChl < 0] <- NA
      chl_xyv <- cbind(chl_xyv, as.vector(intChl[ix, iy]))

      intNDVI <- average_xy(xyv_ndvi[, c(1, 2, i)], output_x = unique(Bat_xyv[, 1]), output_y = unique(Bat_xyv[, 2]))[[3]]
      # intNDVI[intNDVI < 0] <- NA
      ndvi_xyv <- cbind(ndvi_xyv, as.vector(intNDVI[ix, iy]))
    }
  } else if(method == "interpolate"){
    chl_xyv <- average_xyt(input_xyv = xyv_chl,
                           input_t = input_t,
                           output_xy = Bat_xyv[, -3],
                           output_t = input_t)
    chl_xyv[chl_xyv < 0] <- NA

    ndvi_xyv <- average_xyt(input_xyv = xyv_ndvi,
                            input_t = input_t,
                            output_xy = Bat_xyv[, -3],
                            output_t = input_t)

    colnames(chl_xyv) <- colnames(ndvi_xyv)  <- c("longitude", "latitude", as.character(input_t))

    chl_xyv <- as.data.frame(chl_xyv)
    ndvi_xyv <- as.data.frame(ndvi_xyv)
  }



  chl_xyv <- subset(chl_xyv, subset = longitude > min_x & longitude <= max_x & latitude > min_y & latitude <= max_y)
  colnames(chl_xyv) <- c("longitude", "latitude", as.character(input_t))

  ndvi_xyv <- subset(ndvi_xyv, subset = longitude > min_x & longitude <= max_x & latitude > min_y & latitude <= max_y)
  colnames(ndvi_xyv) <- c("longitude", "latitude", as.character(input_t))



  ## Get marsdiep area
  Bat_xyv <- subset(Bat_xyv, subset = longitude > min_x & longitude <= max_x & latitude > min_y & latitude <= max_y)


  ii <- !is.na(Bat_xyv[["depth"]])
  chl_xyv <- chl_xyv[ii, ]
  ndvi_xyv <- ndvi_xyv[ii, ]
  Bat_xyv <- Bat_xyv[ii, ]

  output_xy <- Bat_xyv[, 1:2]
  row.names(output_xy) <- NULL

  output_t <- input_t

  # browser()
  output_xy <- mask_shape(output_xy, shape=Shape)
  non_mask_points <- which(!is.na(output_xy$mask))
  output_xy <- output_xy[non_mask_points, ]
  chl_xyv <- chl_xyv[non_mask_points, ]
  ndvi_xyv <- ndvi_xyv[non_mask_points, ]
  Bat_xyv <- Bat_xyv[non_mask_points, ]


  return(list(Bat_xyv    =  Bat_xyv,
              chl_xyv    =  chl_xyv,
              ndvi_xyv   =  ndvi_xyv,
              input_t    =  input_t,
              output_t   =  output_t,
              output_xy  =  output_xy,
              bbox       =  bbox
  ))

}



transpose_sparse_data <- function(df, ...){
  df_mt <- t(df)
  df_mt <- df_mt[-c(1, 2), ]
  return(as.matrix(df_mt))
}

spatial_operation <- function(X, fun, ...){
  n <- dim(X)[1]; m <- dim(X)[2]
  if(m <= 1){
    out <- as.vector(apply(X[, drop = FALSE], MARGIN = 1, FUN = fun,...))
    # out <- matrix(out, ncol = 1)
  } else {
    out <- apply(X, MARGIN = 2, FUN = fun, ...)
  }
  return(out)
}


load_irradiance <- function(dir, file, bbox,
                            output_xy, output_t,...){
  Irradiance <- read_KNMI(dir  = dir,
                         file = file)[ ,c("station", "datetime", "radiation")]


  stat_idx <- subset(meta(Irradiance)$stations, subset = longitude >= bbox[1] & longitude <= bbox[3] & latitude >= bbox[2] & latitude <= bbox[4]+0.2)
  # print(stat_idx)
  # subset two station corresponding to grid
  Irr <- subset(Irradiance, subset =  station %in% stat_idx$station)

  Weather.stations <- meta(Irr)$station
  # Weather.stations <- Weather.stations[Weather.stations$Wadden, ]

  Weather.stations <- meta(Irr)$station
  # knitr::kable(Weather.stations)

  Watt.m2_to_uE.m2.s <- 4.57*0.45
  Irr$par <- Irr$radiation * Watt.m2_to_uE.m2.s


  Irr <- merge(Irr, Weather.stations[, c("station", "longitude", "latitude")])
  Irr <- Irr[, c("longitude", "latitude", "datetime", "par")]

  Irrad <- interpolate_xyt(input_xytv = na.omit(Irr),
                           output_xy = output_xy,
                           output_t = output_t)


  ## Need to attach metadata like weather station, etc as attributes
  attributes(Irrad)$Weather.stations <- Weather.stations

  return(Irrad)

}


reshape_water_height_wide <- function(RWS_file,  ...){

  stats <- attributes(RWS_file)$stations

  HeightWad <- subset(RWS_file, subset = variable == "Height" &
                        value < 1e5          &
                        !is.na(datetime))


  HeightWad <- HeightWad[order(HeightWad$station, HeightWad$datetime),
                         c("station", "datetime", "value")]


  HH <- average_timeseries(HeightWad,
                           avgOver = "min",
                           avgTime = 30,
                           by = "station")
  colnames(HH)[3] <- "height"

  WadHeightHR <- reshape(data = HH, direction = "wide", timevar = "station", idvar = "datetime")
  colnames(WadHeightHR) <- sub("height.", "", colnames(WadHeightHR))

  attributes(WadHeightHR) <- c(attributes(WadHeightHR),
                               attributes(HeightWad)[c("variables", "stations", "datasource", "EPSG",
                                                       "file", "processing", "fun")])
  attributes(WadHeightHR)$format <- "wide"

  return(WadHeightHR)

}



load_water_height <- function(RWS_file,
                              output_xy,
                              output_t, ...){


  WadHeightHR <- reshape_water_height_wide(RWS_file)


  WadHeightHR_Long <- dt_tolong(WadHeightHR, vname = "Height")
  WadHeightHR_Long_station <- meta(WadHeightHR_Long)$stations


  WadHeightHR_Long_merged_data <- merge(WadHeightHR_Long, y = WadHeightHR_Long_station[, c("station", "longitude", "latitude")])


  WadHeightHR_Long_merged_data <- WadHeightHR_Long_merged_data[, c("station", "longitude", "latitude", "datetime", "Height")]
  head(WadHeightHR_Long_merged_data)

  WadHeightHR_Long_merged_data$Height[WadHeightHR_Long_merged_data$Height > 1000] <- NA
  # meta(xx) <- meta(xx_orig)

  WHeight_RWS <- average_timeseries(input = WadHeightHR_Long_merged_data, avgOver = "hour", avgTime = 1,
                                    by = c("longitude", "latitude"),
                                    value = "Height")


  WHeight <- interpolate_xyt(input_xytv = WHeight_RWS,
                             output_xy = output_xy,
                             output_t = output_t,
                             nmean = nrow(WadHeightHR_Long_station),
                             ...)


  return(WHeight)

}


load_water_temperature <- function(RWS_file, output_xy,
                                   output_t,...){

  stats <- attributes(RWS_file)$stations

  TempWad <- subset(RWS_file, subset = variable == "T" &
                      value < 1e5     &
                      !is.na(datetime))


  TempWad <- TempWad[, c('station', 'datetime', 'value')]

  Tnum <- table(TempWad$station)

  HRstation <- names(Tnum)[Tnum > 10000]

  TempWadHR <- subset(TempWad,
                      subset = station %in% HRstation)

  statHR <- subset(stats,
                   subset = station %in% HRstation)

  THR <- average_timeseries(TempWadHR,
                            avgOver = "hour",
                            avgTime = 1,
                            by = "station")
  colnames(THR)[3] <- "T"


  WadTempHR <- reshape(data = THR,
                       direction = "wide",
                       timevar = "station",
                       idvar = "datetime")

  attributes(WadTempHR)$stations <- statHR
  cn <- colnames(WadTempHR)[-1]
  cn <- sub("T.", "", cn)
  colnames(WadTempHR)[-1] <- cn
  attributes(WadTempHR)$format <- "wide"

  WadTempHR <- WadTempHR[order(WadTempHR$datetime), ]


  WT <- dt_tolong(WadTempHR, vname = "Temperature")
  Temp_stations <- meta(TempWad)$stations


  WT <- merge(Temp_stations[, c("station", "longitude", "latitude")],
              WT[, c("station", "datetime", "Temperature")])

  WT <- WT[order(WT$station, WT$datetime), ]

  WT$datetime <- as.POSIXct(WT$datetime)


  WTemp    <- interpolate_xyt(input_xytv = WT[, -1],
                              output_xy = output_xy,
                              output_t = output_t,
                              nmean = 15)

  return(list(WTemp = WTemp, WT = WT))
}


load_photosynthesis <- function(PPS, WT, biogeo_data,
                                WChl,
                                WTemp,
                                water_mask,
                                output_xy,
                                output_t,
                                ...){

  PPS$date <- as.POSIXct(PPS$date)

  WPI_stations <- meta(PPS)$stations

  variables <- subset(meta(PPS)$variables,
                      subset = variable %in% c( "alpha", "eopt", "ps"))


  colnames(WT)[4] <- "date"

  select <- c("station", "date", "alpha", "eopt", "ps", "Chla_ug.L", "Kd")
  #select <- c("station", "datetime", "alpha", "eopt", "pmax", "Chla_ug.L", "Kd")

  PS_data <- merge(WT[, c("station", "date", "Temperature")],
                   PPS [, select], # PS parameters
                   all = TRUE,
                   by = c("station", "date"))

  PS_data <- merge(WPI_stations[, c("station", "longitude", "latitude")],
                   PS_data)
  PS_data$date <- as.POSIXct(PS_data$date)


  # Eopt
  #select <- c("station", "date", "alpha", "eopt", "ps", "Chla_ug.L", "Kd")
  PS_data <- merge(biogeo_data[, c("station", "date", "Temperature")],
                   PPS    [, select],    # PS parameters
                   all = TRUE,
                   by  = c("station", "date"))
  PS_data <- merge(WPI_stations[, c("station", "longitude", "latitude")],
                   PS_data)
  PS_data$date <- as.POSIXct(PS_data$date)

  LM <- with(PS_data,  lm(eopt~Temperature))
  # abline()

  cc.eopt <- coef(LM)

  # optimal light intensity for all output stations
  WEopt  <- data.frame(date = output_t,
                       cc.eopt[1] + cc.eopt[2]*transpose_sparse_data(WTemp))


  # Ps
  Q10fun <- function(p,   # Q10 parameters
                     T)   # Temperature
    p[1] * p[2]^((10-T)/10)

  # data to fit the Q10 function to
  PPS   <- na.omit(PS_data[,c("ps", "Temperature")])

  # function that returns the sum of squared residuals of model and data
  SRfun <- function(p)
    sum( (PPS$ps-Q10fun(p, PPS$Temperature))^2)

  # nonlinear model fitting
  ZZ <- nlm(f = SRfun,
            p = c(mean(PPS$ps), 2) )

  # best-fit parameters
  cc.ps <- ZZ$estimate



  WPs  <- cbind(date = output_t,
                t(cc.ps[[1]] * cc.ps[[2]]^((10- WTemp[,-c(1, 2)])/10) *
                    WChl[, -c(1, 2)]))


  # Alpha
  PS_data$alpha[PS_data$alpha > 0.2] <- NA

  WAlpha  <- cbind(date = output_t,
                   t(mean(PS_data$alpha, na.rm=TRUE) * water_mask * WChl[, -c(1, 2)]))


  WKd <- interpolate_xyt(
    input_xytv = PS_data[, c("longitude", "latitude", "date", "Kd")],
    output_xy  = output_xy,
    output_t   = output_t,
    nmean      = 6)


  return(list(PS_data    = PS_data,
              WPs        = WPs,
              WEopt      = WEopt,
              WAlpha     = WAlpha,
              WKd        = WKd))
}


create_water_mask <- function(bat, WHeight){
  mask_array <- bat

  mask_array$depth <- mask_array$depth * -1

  immersion_mask <- apply(WHeight[, -c(1:2)], MARGIN = 2,
                          FUN = function(x) ifelse(mask_array$depth < x, 1, 0))

  immersion_mask_na <- immersion_mask
  # immersion_mask_na[immersion_mask_na == 0] <- NA

  return(immersion_mask_na)
}

calculate_pelagic_primary_production <- function(bathymetric_file,
                                                 WKd, Irrad, WAlpha,
                                                 WEopt, WPs, WHeight){

  tKd     <- t(as.matrix(WKd[, -c(1, 2)]))
  tIrrad  <- t(as.matrix(Irrad[, -c(1, 2)]))
  tAlpha  <- as.matrix(WAlpha[, -1])
  tEopt   <- as.matrix(WEopt[, -1])
  tPs     <- as.matrix(WPs[, -1])
  tHeight <- t(as.matrix(WHeight[, -c(1, 2)]))


  # system.time(
  #   ppPel <- intPP_mixed(bathymetric_file$depth,
  #                        tKd,
  #                        tIrrad,
  #                        tAlpha,
  #                        tEopt,
  #                        tPs,
  #                        tHeight
  #   )
  # )

  system.time(
    ppPel <- rs_mixed_PP(zmax = bathymetric_file$depth,
                         kz = tKd,
                         par = tIrrad,
                         alfa = tAlpha,
                         eopt = tEopt,
                         pmax = tPs,
                         height = tHeight
    )
  )

  return(ppPel)
}


load_sediment <- function(sediment_data,
                          poro,
                          output_xy,
                          output_t
                          ) {

  Sediment_grid <- interpolate_xy(
    input_xyv = sediment_data[ ,c("longitude", "latitude", "mdGrain")],
    output_xy = output_xy)

  Sediment_grid$sand <- interpolate_xy(
    input_xyv = sediment_data[ ,c("longitude", "latitude", "sand")],
    output_xy = output_xy)$sand

  Sediment_grid$silt <- interpolate_xy(
    input_xyv = sediment_data[ ,c("longitude", "latitude", "silt")],
    output_xy = output_xy)$silt

  ## Porosity data
  # linear relationship
  LM <- lm(poro$porosity~I(poro$sandperc/100))

  coeff.por <- coef(LM)

  # Sediment porosity based on sand fraction
  Sediment_grid$poro <- coeff.por[[1]] + coeff.por[[2]]*Sediment_grid$sand/100

  # g solid sediment per cm3 of bulk
  rho  <- 2.5             # sediment dry density g/cm3
  Sediment_grid$dens <- (1 - Sediment_grid$poro)*rho

  # kd 4000/m per g/cm3 for sand and 45000/m for silt
  Sediment_grid$Kd <- Sediment_grid$dens*Sediment_grid$sand/100 * 4000 +
    Sediment_grid$dens*Sediment_grid$silt/100 * 45000

  return(Sediment_grid)
}


load_photosynthesis_benthic <- function(BPS, porosity, biogeo_data,
                                        BChl, WTemp,
                                        immersion_mask_na,
                                        output_xy, output_t){
  # the Ems benthicstations
  BPI_stations  <- meta(BPS)$stations
  BPI_stations

  B_PS <- BPS[ , c("station", "date", "alpha", "eopt", "ps", "Chl_sed_ug.g")]

  B_PS <- na.omit(B_PS)


  # Show information about the variables
  variables <- subset(meta(B_PS)$variables,
                      subset = variable %in% c( "alpha", "eopt", "ps"))


  B_PS <- merge(B_PS, porosity, by.x="station", by.y="nr")
  B_PS$Chl_sed_mg.m3 <- B_PS$Chl_sed_ug.g * B_PS$ug.g2mg.m3
  B_PS$date <- as.POSIXct(B_PS$date)


  WT <- biogeo_data[, c("station", "date", "Temperature")]
  WT <- merge(meta(biogeo_data)$stations, WT)
  WT <- na.omit(WT)

  Tmean <-
    with (WT, aggregate(
      x   = list(Temperature = Temperature),
      by  = list(date = date),
      FUN = mean, na.rm=TRUE))

  B_PS$Temperature <- approx(x    = julian(Tmean$date),
                             y    = Tmean$Temperature,
                             xout = julian(B_PS$date))$y


  B_PS$eopt[B_PS$eopt > 1e4] <- NA

  LM <- with(B_PS,  lm(eopt~Temperature))
  # summary(LM)

  cc.eopt <- coef(LM)

  # cat(length(output_t), "\n",
  #     nrow(transpose_sparse_data(WTemp)), "\n",
  #     nrow(t(immersion_mask_na_benth)))

  # optimal light intensity for all output stations
  BEopt <- data.frame(date = output_t,
                      cc.eopt[1] + cc.eopt[2]*transpose_sparse_data(WTemp) * t(immersion_mask_na))

  ## fitting function for Q10
  Q10fun <- function(p, T)
    p[1]*p[2]^((10-T)/10)

  PP <- na.omit(B_PS[,c("ps", "Temperature")])

  fun <- function(p)
    sum( (B_PS$ps-Q10fun(p, B_PS$Temperature))^2)

  ZZ <- nlm(f=fun, p=c(mean(B_PS$ps), 2) )

  ## best-fit parameters
  cc.ps <- ZZ$estimate

  BPs  <- cbind(date = output_t,
                t(cc.ps[[1]] * cc.ps[[2]]^((10- WTemp[,-c(1, 2)])/10) *
                    immersion_mask_na * BChl))

  BAlpha  <- cbind(date = output_t,
                   t(mean(B_PS$alpha, na.rm=TRUE) * immersion_mask_na * BChl))

  return(list(B_PS       = B_PS,
              BPs        = BPs,
              BEopt      = BEopt,
              BAlpha     = BAlpha))
}

estimate_chla_from_ndvi <- function(WNDVI,
                                    immersion_mask_na,
                                    a = 532, b = 48){

  WNDVI[, -c(1, 2)] <- apply(WNDVI[, -c(1, 2)], MARGIN = 2,
                             function(x){
                               ifelse(x > 1, 1, # overboard => forest
                                      ifelse(x < -1, -1, # clear water
                                             x)) # water + vegetation
                             })

  # WNDVI[WNDVI[, -c(1, 2)] >  1,  ] <- 1
  # WNDVI[WNDVI[, -c(1, 2)] < -1, ]  <- -1

  BChl <- a * (immersion_mask_na * WNDVI[,-c(1, 2)]) + b

  # use NDVI top of atmosphere as the surface reflectance seem to not work well
  # BChl <- a * (immersion_mask_na * WNDVI_rhot_Ems[,-c(1, 2)]) + b

  BChl[BChl < 0] <- NA

  return(BChl)
}

create_water_mask_benthic <- function(bathymetric_file, WHeight){
  mask_array <- bathymetric_file

  mask_array$depth <- mask_array$depth * -1

  immersion_mask <- apply(WHeight[, -c(1:2)], MARGIN = 2,
                          FUN = function(x) ifelse(mask_array$depth < x, 1, 0))

  immersion_mask_na <- immersion_mask_na_benth <- immersion_mask
  # immersion_mask_na[immersion_mask_na == 0] <- NA
  immersion_mask_na_benth[immersion_mask_na_benth == 1] <- NA
  immersion_mask_na_benth[immersion_mask_na_benth == 0] <- 1

  return(immersion_mask_na_benth)
}


calculate_benthic_primary_production <- function(bathymetric_file,
                                                 Sediment_grid,
                                                 WKd, Irrad, BAlpha,
                                                 BEopt, BPs, WHeight,
                                                 zn = 0.002  # depth of chlorophyll layer
){
  Rad <- rad_bot(bathymetric_file$depth,
                 t(as.matrix(WKd    [, -c(1, 2)])),
                 t(as.matrix(Irrad  [, -c(1, 2)])),
                 t(as.matrix(WHeight[, -c(1, 2)])))


  # Benthic primary production,
  # exponentially declining chlorophyll concentration, a function of silt fraction
  # ppBen <- intPP_exp(as.vector(rep(zn, times = nrow(bathymetric_file))),
  #                    as.vector(Sediment_grid$Kd),
  #                    as.vector(Sediment_grid$silt/100),
  #                    as.matrix(Rad),
  #                    as.matrix(BAlpha      [, -1]),
  #                    as.matrix(BEopt       [, -1]),
  #                    as.matrix(BPs         [, -1]))
  
  ppBen <- rs_exp_PP(zmax = as.vector(rep(zn, times = nrow(bathymetric_file))),
                     kz = as.vector(Sediment_grid$Kd),
                     pmud = as.vector(Sediment_grid$silt/100),
                     par = as.matrix(Rad),
                     alfa = as.matrix(BAlpha      [, -1]),
                     eopt = as.matrix(BEopt       [, -1]),
                     pmax = as.matrix(BPs         [, -1]))

  return(list(Rad = Rad, ppBen = ppBen))

}

rs_files <- function(rs_dir, year, region){
  basepath <- paste0(rs_dir, year, "/", year, "/")
  files <- switch(region,
                  "Marsdiep" = list.files(basepath, pattern = glob2rx(pattern = "*T31UFU*L2W*.nc"), recursive = TRUE, full.names = TRUE),
                  "Vlieland" = list.files(basepath, pattern = glob2rx(pattern = "*T31UFV*L2W*.nc"), recursive = TRUE, full.names = TRUE),
                  "Ems" = list.files(basepath, pattern = glob2rx(pattern = "*T32ULE*L2W*.nc"), recursive = TRUE, full.names = TRUE))
  return(files)
}

rs_files_lisemijn <- function(rs_dir, year, region){
  basepath <- paste0(rs_dir, year, "/", "/")
  files <- switch(region,
                  "Marsdiep" = list.files(basepath, pattern = glob2rx(pattern = "*T31UFU*L2W*.nc"), recursive = TRUE, full.names = TRUE),
                  "Vlieland" = list.files(basepath, pattern = glob2rx(pattern = "*T31UFV*L2W*.nc"), recursive = TRUE, full.names = TRUE),
                  "Ems" = list.files(basepath, pattern = glob2rx(pattern = "*T32ULE*L2W*.nc"), recursive = TRUE, full.names = TRUE),
                  "Ems_2" = list.files(basepath, pattern = glob2rx(pattern = "*T31UGV*L2W*.nc"), recursive = TRUE, full.names = TRUE))
  return(files)
}


run_pelagic <- function(files,
                        year, region, RWS_parse, PPS,
                        biogeo, knmi_dir = "D:/KMNI"){

  # Sentinel-2 data
  RS_data <- read_acolite_output(files)

  out <- create_spatio_temporal(files = files,
                                bathymetric_file = Wad_depth_HR,
                                RS_data = RS_data)

  output_xy <- out$output_xy
  output_t  <- out$output_t

  if(is.null(knmi_dir)) stop("Please provide directory where KNMI data are stored.")

  Irrad <- load_irradiance(dir       = knmi_dir,
                           file      = paste0("KNMI_Wad_", year, ".txt"),
                           bbox      = out$bbox,
                           output_xy = output_xy,
                           output_t  =  output_t)
  meanPAR <- apply(transpose_sparse_data(Irrad), MARGIN = 2, FUN = mean, na.rm = TRUE)

  WHeight <- load_water_height(RWS_parse, output_xy = output_xy, output_t = output_t)

  rangeH <- apply(transpose_sparse_data(WHeight), MARGIN = 2, FUN = function(x) diff(range(x, na.rm = TRUE)))

  WTemperature <- load_water_temperature(RWS_parse,
                                         output_xy = output_xy,
                                         output_t  = output_t)

  WTemp <- WTemperature$WTemp

  meanT <- apply(transpose_sparse_data(WTemp), MARGIN = 2, FUN = mean)

  immersion_mask_na <- create_water_mask(out$Bat_xyv, WHeight)

  PS_data <- load_photosynthesis(PPS, WTemperature$WT,
                                 biogeo, out$chl_xyv,
                                 WTemp,
                                 immersion_mask_na,
                                 output_xy = output_xy,
                                 output_t = output_t)
  WAlpha <- PS_data$WAlpha
  WPs    <- PS_data$WPs
  WEopt  <- PS_data$WEopt
  WKd    <- PS_data$WKd

  meanAlpha <- apply(WAlpha[,-1],
                     MARGIN = 2,
                     FUN    = mean,  na.rm = TRUE)
  meanPs <- apply(WPs[,-1],
                  MARGIN = 2,
                  FUN    = mean, na.rm = TRUE)
  meanEopt <- apply(WEopt[,-1],
                    MARGIN = 2,
                    FUN    = mean, na.rm = TRUE)

  ppPel <- calculate_pelagic_primary_production(out$Bat_xyv,
                                                Irrad   = Irrad,
                                                WKd     = WKd,
                                                WEopt   = WEopt,
                                                WPs     = WPs,
                                                WAlpha  = WAlpha,
                                                WHeight = WHeight)

  return(list(ppPel = ppPel, out = out, Irrad = Irrad, WKd = WKd, WHeight = WHeight,
              WEopt = WEopt, WPs = WPs, WAlpha = WAlpha, WTemperature = WTemperature))
}

run_benthic <- function(year, region, out, WHeight, Wad_sediment, BPS, poro,
                        biogeo, WTemp, Irrad, WKd){

  output_xy <- out$output_xy
  output_t <- out$output_t

  immersion_mask_na_benth <- create_water_mask_benthic(out$Bat_xyv, WHeight)
  BChl <- estimate_chla_from_ndvi(out$ndvi_xyv, immersion_mask_na_benth)

  sediment <- subset(Wad_sediment,
                     subset = longitude >= out$bbox[1] & longitude <= out$bbox[3] &
                       latitude >= out$bbox[2]  & latitude <= out$bbox[4])

  Sediment_grid <- load_sediment(sediment, poro, output_xy, output_t)


  BS_data <- load_photosynthesis_benthic(BPS, poro, biogeo,
                                         BChl, WTemp,
                                         immersion_mask_na_benth,
                                         output_xy, output_t)

  BEopt  <- BS_data$BEopt
  BPs    <- BS_data$BPs
  BAlpha <- BS_data$BAlpha

  meanEopt <- apply(BEopt[,-1],
                    MARGIN = 2,
                    FUN    = mean, na.rm = TRUE)
  meanPs <- apply(BPs[,-1],
                  MARGIN = 2,
                  FUN    = mean, na.rm = TRUE)

  meanAlpha <- apply(BAlpha[,-1],
                     MARGIN = 2,
                     FUN    = mean, na.rm = TRUE)

  ppBenthic <- calculate_benthic_primary_production(out$Bat_xyv,
                                                    Sediment_grid = Sediment_grid,
                                                    Irrad   = Irrad,
                                                    WKd     = WKd,
                                                    BEopt   = BEopt,
                                                    BPs     = BPs,
                                                    BAlpha  = BAlpha,
                                                    WHeight = WHeight)
  ppBen <- ppBenthic$ppBen
  Rad   <- ppBenthic$Rad

  return(list(ppBen = ppBen, Rad = Rad, Sediment_grid = Sediment_grid,
              BEopt = BEopt, BPs = BPs, BAlpha = BAlpha))

}

run_primary_production <- function(files = NULL,
                                   year = NULL,
                                   region = NULL,
                                   bathymetric_file = NULL,
                                   RWS_parse = NULL,
                                   knmi_dir = NULL,
                                   Wad_sediment = NULL,
                                   PPS = NULL,
                                   BPS = NULL,
                                   poro = NULL,
                                   biogeo = NULL,
                                   ...){

  if(is.null(bathymetric_file)) stop("Please provide bathymetry file...")

  load(file = bathymetric_file)

  if(is.null(RWS_parse) ){
    stop("Please provide Rijkwaterstraat data.")
  }

  pel <- run_pelagic(files = files,year = year, region = region, RWS_parse,
                     PPS, biogeo, knmi_dir = knmi_dir)
  ppPel    <- pel[["ppPel"]]
  out      <- pel[["out"]]
  output_xy <- out$output_xy
  # output_t <- out$output_t
  Irrad    <- pel[["Irrad"]]
  WKd      <- pel[["WKd"]]
  WHeight  <- pel[["WHeight"]]
  WEopt    <- pel[["WEopt"]]
  WPs      <- pel[["WPs"]]
  WAlpha   <- pel[["WAlpha"]]
  WTemperature <- pel[["WTemperature"]]
  WTemp <- WTemperature$WTemp

  ben <- run_benthic(year = year, region = region,
                     out, WHeight, Wad_sediment, BPS, poro, biogeo,
                     WTemp, Irrad, WKd)

  return(list(pel = pel, ben = ben))
}

year_average_wadden <- function(year, dir = "D:/RS_analysis_data/PP_results/", regions = NULL){
  res_env <- new.env()
  Pelagic_xy_wad <- Benthic_xy_wad <- Pelagic_t_wad <- Benthic_t_wad <- NULL

  if(is.null(regions)){
    regions <- c("Marsdiep", "Vlieland", "Ems")
  }


  for(rr in regions){

    load(paste0(dir, rr, "_", year, ".rda"),
         envir = res_env)

    res_env$Pelagic_t <- data.frame(date = res_env$res$pel$out$output_t, pprod = apply(res_env$res$pel$ppPel, MARGIN = 1, FUN = mean, na.rm = TRUE))
    res_env$Benthic_t <- data.frame(date = res_env$res$pel$out$output_t, pprod = apply(res_env$res$ben$ppBen, MARGIN = 1, FUN = mean, na.rm = TRUE))

    res_env$Pelagic_xy <- data.frame(res_env$res$pel$out$Bat_xyv, # longitude, latitude, depth
                                     ppPel = apply(res_env$res$pel$ppPel, MARGIN = 2, FUN = mean, na.rm = TRUE))

    res_env$Benthic_xy <- data.frame(res_env$res$pel$out$Bat_xyv,
                                     ppBen = apply(res_env$res$ben$ppBen, MARGIN = 2, FUN = mean, na.rm = TRUE))
    Pelagic_xy_wad <- rbind(Pelagic_xy_wad, res_env$Pelagic_xy)
    Benthic_xy_wad <- rbind(Benthic_xy_wad, res_env$Benthic_xy)

    Pelagic_t_wad <- rbind(Pelagic_t_wad, res_env$Pelagic_t)
    Benthic_t_wad <- rbind(Benthic_t_wad, res_env$Benthic_t)
  }

  rm(res_env)

  Pelagic_t_wad <- Pelagic_t_wad[order(Pelagic_t_wad$date), ]
  Benthic_t_wad <- Benthic_t_wad[order(Benthic_t_wad$date), ]

  return(list(Pelagic_xy_wad = Pelagic_xy_wad,
              Benthic_xy_wad = Benthic_xy_wad,
              Pelagic_t_wad  = Pelagic_t_wad,
              Benthic_t_wad  = Benthic_t_wad))

}

