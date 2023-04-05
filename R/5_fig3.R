#### Make figures for manuscript
library(tidyverse)
library(lubridate)
library(cowplot)
library(sf)
library(rmapzen) # needed for mz_ functions to get roads
theme_set(theme_cowplot())

source('R/infection-estimation-fxns-jags.R')
source('R/data-helper-fxns.R')
source('R/load-results.R')

# Fig 3 zip cumulative map results ----------------------------------------
load('raw-data/spatial-files/zip_geometries.rda')

mz_set_tile_host_nextzen(key="hxNDKuWbRgetjkLAf_7MUQ")
## Get roads for maps
get_vector_tiles <- function(bbox){
  mz_box=mz_rect(bbox$xmin,bbox$ymin,bbox$xmax,bbox$ymax)
  mz_vector_tiles(mz_box)
}
zcta_geom <- st_union(zip_cov_df %>% 
                        left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
                        pull(geometry))
zcta_bbox <- st_bbox(zcta_geom)
zcta_vector_tiles <- get_vector_tiles(zcta_bbox)
zcta_roads <- as_sf(zcta_vector_tiles$roads) 

I35<- zcta_roads %>%
  mutate(st_transform(geometry, st_crs(zip_geom))) %>% 
  filter(ref== "I 35" | ref=="I 35;US 190" | ref=="I 35;US 290" | ref=="I 35;US 77") %>% # | ref=="I 35-H Business") %>%
  pull(geometry) %>%
  st_union() %>%
  st_transform(st_crs(zcta_geom)) %>%
  st_intersection(zcta_geom)

US183 = zcta_roads %>%
  mutate(st_transform(geometry, st_crs(zip_geom))) %>% 
  # filter(grepl('183', ref)) %>% 
  filter(ref=="US 183" | ref=="US 183;FM 20" | ref=="US 183;TX 80" | ref=="US 90;US 183" | ref== "US 183;US 190;US 281"
         | ref=="US 183;US 190;US 281;FM 580" | ref=="US 183;US 190;US 281;Truck" | ref=="US 183;US 190" | id=="2c2066c3b6fafe2574463154f5cc877d" ) %>%
  pull(geometry) %>%
  st_union() %>%
  st_transform(st_crs(zcta_geom)) %>%
  st_intersection(zcta_geom)

plot_zip_maps <- function(sf_obj, fill_col, plot_title, scales_func){
  fill_col <- enquo(fill_col)
  # browser()
  sf_obj %>% 
    ggplot() +
    geom_sf(aes(geometry=geometry, fill=!!fill_col), 
            size = 0.05, color="black")+ # , show.legend = FALSE
    scale_fill_gradient(low="lightyellow", high="slateblue4", labels = scales_func)+
    guides(fill = guide_colourbar(frame.colour = "black", ticks.colour="black"))+
    geom_sf(data = I35, col = "black", size=0.3)+
    geom_sf(data = US183, col = "black", size=0.3)+
    annotate(geom="text", x=-97.68, y=30.50, label="I-35", size=3)+
    annotate(geom="text", x=-97.80, y=30.50, label="US 183", size=3)+
    labs(title = plot_title, fill = NULL)+
    theme_map() +
    theme(plot.title = element_text(hjust = 0.5,size = 14), plot.title.position = 'plot') +
    NULL
}


## Case rate
zip_cov_df %>% 
  mutate(case_rate = cases/population*100000) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = case_rate, 
                plot_title = 'Cases per 100k', 
                scales::label_comma()) -> zip_case_rate_map
zip_case_rate_map

## Hosp rate
age_zip_df %>% 
  group_by(zip) %>% 
  summarize(n_admits = sum(n_admits)) %>% 
  left_join(zip_cov_df) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  mutate(admit_rate = n_admits/population*100000) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = admit_rate, 
                plot_title = 'Hospitalizations per 100k', 
                scales::label_comma()) -> zip_hosp_rate_map
zip_hosp_rate_map

## SVI
zip_cov_df %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = svi, 
                plot_title = 'Social vulnerability index (SVI)', 
                scales::label_comma(accuracy = .1)) -> zip_svi_map
zip_svi_map

## Infection hospitalization rate
age_zip_df %>% 
  group_by(zip) %>% 
  summarize(ihr = sum(ihrmean*population)/sum(population)) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = ihr, 
                plot_title = 'Infection hospitalization rate (IHR)', 
                scales::label_percent(accuracy = 1)) -> zip_ihr_map
zip_ihr_map


## Zip cumulative infections
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(zip, samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  summarize(meaninf = mean(inf),
            lowerinf = quantile(inf, probs = 0.025),
            upperinf = quantile(inf, probs = 0.975)) %>% 
  left_join(zip_cov_df) %>% 
  mutate_at(vars(meaninf:upperinf), ~./population) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) -> zip_cuminf_df 

zip_cuminf_df %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = meaninf, 
                plot_title = 'Cumulative infections (%)', 
                scales::label_percent(accuracy = 1)) -> zip_cum_inf_map
zip_cum_inf_map


## Zip reporting rate
age_zip_inf_df %>% 
  unnest(est_infections) %>% 
  group_by(zip, age_group_paper) %>% 
  mutate(samp = seq_along(est_infections)) %>% 
  group_by(zip, samp) %>% 
  summarize(inf = sum(est_infections)) %>% 
  left_join(zip_cov_df) %>% 
  mutate(cases = ifelse(inf < cases, inf, cases)) %>% 
  mutate(rr = samp_conditional_reporting_rate(detected_cases = cases, 
                                              estimated_infections = inf, 
                                              num_samps = n())) %>% 
  summarize(meanrr = mean(rr),
            lowerrr = quantile(rr, probs = 0.025),
            upperrr = quantile(rr, probs = 0.975)) %>% 
  left_join(zip_geom, by = c('zip' = 'ZIP')) -> zip_rr_df

zip_rr_df %>% 
  plot_zip_maps(sf_obj = .,
                fill_col = meanrr, 
                plot_title = 'Reporting rate (%)', 
                scales::label_percent(accuracy = 1)) -> zip_rr_map
zip_rr_map

##Combine into grid
plot_grid(zip_case_rate_map, zip_hosp_rate_map, 
          zip_svi_map, zip_ihr_map,
          zip_cum_inf_map, zip_rr_map, 
          nrow = 3, labels = 'AUTO', align = 'hv') -> f3_zip_maps
save_plot('figs/fig3_zip_maps.png', f3_zip_maps, base_height =8, base_asp = 1.1,bg = "white")
save_plot('figs/fig3.tiff', f3_zip_maps, base_height =8, base_asp = 1.1,bg = "white")

# Zip code level statistics -----------------------------------------------

## Case rates
zip_cov_df %>% 
  mutate(case_rate = cases/population*100000) %>% 
  filter(case_rate == min(case_rate) | case_rate == max(case_rate))

## Hospitalization rates
age_zip_df %>% 
  group_by(zip) %>% 
  summarize(n_admits = sum(n_admits)) %>% 
  left_join(zip_cov_df) %>% 
  mutate(admit_rate = n_admits/population*100000) %>% 
  filter(admit_rate == min(admit_rate) | admit_rate == max(admit_rate))

## IHR
age_zip_df %>% 
  group_by(zip) %>% 
  summarize(ihr = sum(ihrmean*population)/sum(population),
            ihr_lo = sum(ihrlower*population)/sum(population),
            ihr_hi = sum(ihrupper*population)/sum(population)) %>% 
  filter(ihr == min(ihr) | ihr == max(ihr)) %>% 
  mutate(infections = 1/ihr)

## Cumulative infection estimates
zip_cuminf_df %>% 
  filter(meaninf == min(meaninf) | meaninf == max(meaninf))

zip_rr_df %>% 
  filter(meanrr == min(meanrr) | meanrr == max(meanrr))

save.image(file='processed-data/fig3-image.rda')

