# burden gist

library(YFburden)
library(dplyr)
library(magrittr)
library(tidyr)
library(snapalette)
library(KsetupR)
library(R.utils)
library(YFestimation)
library(ggplot2)
library(snapalette)
library(rgdal)

sourceDirectory("Functions", modifiedOnly = FALSE)

# read shapefiles in #
shp0 = rgdal::readOGR("../../shapefiles/Africa_SAmerica_0_smooth.shp",
                      stringsAsFactors = FALSE,
                      encoding = "UTF-8",
                      use_iconv = TRUE)
shp1 = rgdal::readOGR("../../shapefiles/Africa_SAmerica_1_smooth.shp",
                      stringsAsFactors = FALSE,
                      encoding = "UTF-8",
                      use_iconv = TRUE)

# ---------------------------------------------------------------------------
### POPULATION DATA ###

all_res_pop_3d = get_pop_data_3d() 

pop1 = all_res_pop_3d$pop1                                      
P_tot = all_res_pop_3d$P_tot                                   #total populations for each adm and year
p_prop = all_res_pop_3d$p_prop                                   #proportions of population

pop1 %<>% rename(adm0_adm1 = adm1)

pop1 %<>% arrange(year) %>% select(-c(country_code, country))

names(pop1)[3:ncol(pop1)] = paste0("a", 0:100)

pop1 %<>% filter(year<2051)
p_prop %<>% filter(year<2051)
P_tot %<>% filter(year<2051)

# ---------------------------------------------------------------------------
# VACCINATION DATA #
vc2d = readRDS("vacc_1940_1950.RDS")
vc2d = as.data.frame(vc2d)

# expand vc2d to match pop1
vc2d %<>% mutate(year = as.numeric(as.character(year)))
vc_tmp = pop1 %>% filter(!adm0_adm1 %in% vc2d$adm0_adm1)
vc_tmp[, grep("adm0_adm1|year", names(vc_tmp), invert = TRUE)] = 0

vc2d %<>% bind_rows(vc_tmp)
vc2d %<>% unique()
vc2d %<>% filter(!is.na(adm0_adm1))


# ---------------------------------------------------------------------------
# FOI #
foi = read.csv( "FOI_samples.csv", stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
fun_calcBurden = function(fois, year , adm1 = NA, vc2d, pop1) { # Tini's burden function
  
  if(anyNA(adm1)){
    adm1 = names(fois)
  }
  
  foi.mat = matrix(fois, nrow=length(fois), ncol=101, byrow=FALSE)
  
  class(foi.mat) = "numeric"
  
  a.mat = matrix(0:100, nrow=length(fois), ncol=101, byrow=TRUE)
  
  burden = fois * rowSums(exp(-foi.mat*a.mat) * 
                            (1-vc2d[vc2d$year==year & vc2d$adm0_adm1 %in% adm1,-(1:3)]) *
                            pop1[pop1$year==year & pop1$adm0_adm1 %in% adm1, -(1:2)], na.rm=TRUE)
  return(burden)
}

# ---------------------------------------------------------------------------

infections_out = fun_calcBurden(foi[1,], 2018, adm1 = names(foi), vc2d, pop1)

infections_df = data.frame(adm0_adm1 = names(infections_out),
                           year = 2018,
                           infections = as.numeric(infections_out))

infections_df %<>% mutate(adm0 = substr(adm0_adm1, 1, 3))

infections_df %<>% mutate(severe_infections = 0.12*infections,
                          deaths = 0.47*severe_infections)

pop_df = data.frame(adm0_adm1 = pop1$adm0_adm1,
                    year = pop1$year,
                    pop = rowSums(pop1[, -c(1:2)], na.rm = TRUE))

infections_df %<>% left_join(pop_df)

infections_df %<>% mutate(severe_inci = severe_infections/pop,
                          deaths_inci = deaths/pop)
# ---------------------------------------------------------------------------

ggplot(infections_df) + geom_col( aes(x = adm0, y = infections))


# ---------------------------------------------------------------------------
colours = snapalette("Pop", 1000, type = "continuous")

shp1$infections = NA
shp1$infections = infections_df$infections[match(shp1$GID_1, infections_df$adm0_adm1)]

mybreaks= seq(min(log10(shp1$infections), na.rm = TRUE), max(log10(shp1$infections), na.rm = TRUE)+0.01, length.out=1001)

vcols = findInterval(log10(shp1$infections),mybreaks)

plot(shp0, main = "Infections in 2018")
plot(shp1, col=colours[vcols], lty=2, add=T)
plot(shp0, add = TRUE, lwd = 2)
image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
           legend.mar = 3.5)

# ---------------------------------------------------------------------------
shp1$deaths = NA
shp1$deaths = infections_df$deaths[match(shp1$GID_1, infections_df$adm0_adm1)]

mybreaks= seq(min(log10(shp1$deaths), na.rm = TRUE), max(log10(shp1$deaths), na.rm = TRUE)+0.01, length.out=1001)

vcols = findInterval(log10(shp1$deaths),mybreaks)

plot(shp0, main = "Deaths in 2018")
plot(shp1, col=colours[vcols], lty=2, add=T)
plot(shp0, add = TRUE, lwd = 2)
image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
           legend.mar = 3.5)

#sum(shp1$deaths, na.rm = TRUE)

# ---------------------------------------------------------------------------

tmp =  infections_df %>%group_by(adm0) %>% summarise(sum_deaths = sum(deaths, na.rm = TRUE))

shp0$deaths = NA
shp0$deaths = tmp$sum_deaths[match(shp0$GID_0, tmp$adm0)]

mybreaks= seq(min(log10(shp0$deaths), na.rm = TRUE), max(log10(shp0$deaths), na.rm = TRUE)+0.01, length.out=1001)

vcols = findInterval(log10(shp0$deaths),mybreaks)

plot(shp0, main = "Deaths in 2018")
plot(shp0, col=colours[vcols], lty=2, add=T)
plot(shp0, add = TRUE, lwd = 2)
image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
           legend.mar = 3.5)


# ---------------------------------------------------------------------------
shp1$deaths_inci = NA
shp1$deaths_inci = infections_df$deaths_inci[match(shp1$GID_1, infections_df$adm0_adm1)]

mybreaks= seq(min(log10(shp1$deaths_inci), na.rm = TRUE), max(log10(shp1$deaths_inci), na.rm = TRUE)+0.01, length.out=1001)

vcols = findInterval(log10(shp1$deaths_inci),mybreaks)

plot(shp0, main = "Deaths per capita in 2018")
plot(shp1, col=colours[vcols], lty=2, add=T)
plot(shp0, add = TRUE, lwd = 2)
image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
           legend.mar = 3.5)
