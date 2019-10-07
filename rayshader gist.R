```{r rayshader_prep, eval = F, echo = F}
library(rayshader)
library(raster)

dim_out = 1000


blank_raster <- raster(shp1, ncols = dim_out, nrows = dim_out, vals = 0)
each_adm_value <- sapply(1:length(shp1$predictions), function(i){
  #print(i)
  intersect_country <- rasterize(shp1[i, ], crop(blank_raster, extent(shp1[i, ])), mask = TRUE)
  intersect_country[!is.na(intersect_country)] <- shp1$predictions[i]
  intersect_country
  
}, simplify = FALSE)

each_adm_value$fun <- mean

full_shape <- do.call(mosaic, each_adm_value)

plot(full_shape)
```

```{r rayshader_out, eval = F, echo = F}

bath = focal(flip(full_shape, "x"), w=matrix(1,3,3), mean)

tempfilename = tempfile()
png(tempfilename,width = 401,height=401)
par(mar = c(0,0,0,0))
raster::image(bath,axes = FALSE, col = snapalette::snapalette(snapal, 1000, type = "continuous") )
dev.off()
tempmap = png::readPNG(tempfilename)

bath[is.na(bath)] = -1
depth_exaggeration = 100

bath %>%
  as.matrix() %>%
  sphere_shade(texture="desert",
               zscale = 1/depth_exaggeration) %>%
  add_shadow(ray_shade(as.matrix(bath), 
                       zscale = 1/depth_exaggeration), 
             max_darken = 1) %>%
  add_shadow(ambient_shade(as.matrix(bath), 
                           zscale = 1/depth_exaggeration)) %>%
  
  add_overlay(tempmap, alphacolor = 0.5) %>%
  
  plot_3d(as.matrix(bath),
          solid=TRUE,
          soliddepth="auto",
          solidcolor = 'white',
          solidlinecolor = "grey50",
          water = TRUE,
          waterdepth = 0,
          wateralpha = 0.9,
          theta = 90,
          shadow = TRUE,
          background='white',
          zscale = 1/depth_exaggeration,
          windowsize = c(1200, 700)) 

render_snapshot(clear = TRUE)

```

