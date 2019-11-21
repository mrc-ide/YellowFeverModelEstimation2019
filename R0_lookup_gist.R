df = R0_lookup
df = as.data.frame(df)


df %<>% mutate(adm1 = rownames(df))
df %<>% mutate(adm0 = substr(adm1, 1,3))

df %<>% gather(R0, infections, -c(adm1, adm0))
View(df %>% group_by(adm0) %>% summarise(inf = sum(infections)))
View(df %>% group_by(R0) %>% summarise(inf = sum(infections)) %>% mutate(R0 = as.numeric(R0)) %>% arrange(-R0)) 


#-------------------------------------------------------------------------------
tmp = df %>% group_by(adm0) %>% summarise(inf = sum(infections)))

colours = snapalette("Pop", 100, type = "continuous")

shp0$inf = NA
shp0$inf = tmp$inf[match(shp0$GID_0, tmp$adm0)]

mybreaks= seq(min(log10(shp0$inf), na.rm = TRUE), 
              max(log10(shp0$inf), na.rm = TRUE)+0.01, 
              length.out=101)

vcols = findInterval(log10(shp0$inf),mybreaks)

plot(shp0)
plot(shp0,col=colours[vcols], lty=2, add=T)

image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
           legend.mar = 3.5)

#-------------------------------------------------------------------------------
tmp = df %>% filter(R0 == 10) %>% group_by(adm1) %>% summarise(inf = sum(infections))

colours = snapalette("Pop", 100, type = "continuous")

shp1$inf = NA
shp1$inf = tmp$inf[match(shp1$GID_1, tmp$adm1)]

mybreaks= seq(0, 
              max((shp1$inf), na.rm = TRUE)+0.01, 
              length.out=101)

vcols = findInterval((shp1$inf),mybreaks)

plot(shp1)
plot(shp1,col=colours[vcols], lty=2, add=T)

image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
           legend.mar = 3.5)
