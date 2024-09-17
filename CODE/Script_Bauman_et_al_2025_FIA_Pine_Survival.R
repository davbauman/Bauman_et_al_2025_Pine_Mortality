#
# Supplementary material of Bauman et al. 2025, Title, Journal.
# ____________________________________________________________

# David Bauman, 17/09/2024.
# ________________________

# See R code outline to easily navigate this long R script (initially made of several R scripts).
# Species codes: 111, 121, and 131 are Pinus elliottii, P. palustris, and P. taeda, respectively.

# Note that we interchangeably refer to the logit of the survival plateau in a given spatial cell and 
# year as "K_L" and "K_tmp". In the manuscript, we only use K_L (latent K, on the logit scale; see 
# Methods). Because we initially called this variable K_tmp, K_tmp often appears in object and variable
# names. It must therefore be read as synonym of K_L.

# Load packages ----
# **************
pacman::p_load(tidyverse, data.table, lubridate, sf, rstan, brms, tidybayes, marginaleffects,
               ggridges, viridis, rFIA, rnaturalearth, rnaturalearthdata, patchwork,
               gridExtra, doParallel, cowplot, grid,
               adespatial, spdep, adegraphics)

# helper functions:
`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# Extract FIA data for survival modeling ----
# **************************************
# Daniel J. Johnson 2023-02-28

# To download the raw data (this may take a long time)
getFIA(states = c('AL', 'AR', 'FL', 'GA', 'LA', 'MS', 'NC',
                  'OK','SC','TN','TX', 'VA'), 
       dir = "DATA/southeast 2023/")

# read in southeastern US data:
se <- readFIA(dir = "DATA/southeast 2023/") 

# Note: The data downloaded by the getFIA() function above is quite heavy (~5Gb), and
# is therefore not on the github repository (but can easily be generated running the two
# functions above). Instead the DATA subfolder of the repository simply contains what all this
# first section of code generates, i.e. data2_50km_2023.csv and data2_100km_2023.csv.

se.plots <- subset(se$PLOT, se$PLOT$KINDCD > 0 & se$PLOT$KINDCD <4 # only Annual plots
                   & se$PLOT$PLOT_STATUS_CD == 1 & se$PLOT$INVYR < 9999) # only forested plots with proper years

se.cond <- subset(se$COND, se$COND$PLT_CN %in% se.plots$CN) # conditions for plots
se.tree <- subset(se$TREE, se$TREE$PLT_CN %in% se.plots$CN) # tree data in plots
colnames(se.tree)[1] <- "TRE_CN"

# get only potentially needed columns
tr <- subset(se.tree, select = c("TRE_CN","PLT_CN","PREV_TRE_CN","SUBP",'TREE',"CONDID","SPCD","AGENTCD",
                                 "MORTYR", "STATUSCD", "DIA", "PREVDIA","STANDING_DEAD_CD","PREV_STATUS_CD"))
colnames(tr)[10] <- "TRE_STATUSCD"

# grab needed plot columns
pl <- subset(se.plots , select = c("CN","SRV_CN","CTY_CN","PREV_PLT_CN","INVYR","STATECD","UNITCD",
                                   "COUNTYCD","PLOT","PLOT_STATUS_CD",
                                   "MEASYEAR","MEASMON","MEASDAY","REMPER","KINDCD",
                                   "RDDISTCD","WATERCD","LAT","LON","ELEV","ECOSUBCD","MANUAL"))
colnames(pl)[1] <- "PLT_CN"

#grab condition columns
co <- subset(se.cond, select = c("CN","PLT_CN","CONDID","COND_STATUS_CD","COND_NONSAMPLE_REASN_CD","RESERVCD",
                                 "OWNCD","OWNGRPCD","FORINDCD","ADFORCD","FORTYPCD","FLDTYPCD",
                                 "STDAGE","STDSZCD","FLDSZCD","STDORGCD","STDORGSP","SLOPE","ASPECT",
                                 "PHYSCLCD","DSTRBCD1","DSTRBYR1","DSTRBCD2","DSTRBYR2",
                                 "DSTRBCD3","DSTRBYR3","TRTCD1","TRTYR1","TRTCD2","TRTYR2",
                                 "TRTCD3","TRTYR3","BALIVE","FLDAGE"))

colnames(co)[1] <- "COND_CN"

# remove old plots as they are often "kind of strange"
sr <- subset(pl, pl$INVYR >= 1999)

# # # join tree and plot level data
sr2 <- merge(sr, tr, by = "PLT_CN")
# # # join condition, tree, and plot level data
sr3 <- merge(sr2, co, by = c("PLT_CN", "CONDID"))

# remove plantations and focus on natural forests
sr3 <- subset(sr3, sr3$STDORGCD == 0)

# species selection of pine species
sr3 <- subset(sr3, sr3$SPCD == 111 | sr3$SPCD == 121 | sr3$SPCD == 131)
# make unique id for each stem
sr3$tr.tag <- paste0(sr3$STATECD,"_",sr3$COUNTYCD,"_",sr3$PLOT,"_",sr3$SUBP, "_",sr3$TREE)
# check for human caused deaths
hums <- ifelse(is.na(sr3$AGENTCD), FALSE, ifelse(sr3$AGENTCD == 80 | sr3$TRE_STATUSCD == 3, TRUE, FALSE))
hum.rm <- sr3$tr.tag[hums]
#removed stems eventually killed by humans
sr4 <- subset(sr3, sr3$tr.tag %not in% hum.rm)
sr4 <- subset(sr4, sr4$TRE_STATUSCD != 3)

# # plot level unique codes
sr4$plotid <- paste0(sr4$STATECD,"_",sr4$UNITCD,"_",sr4$COUNTYCD,"_",sr4$PLOT)
# make measurement date for time calculations
sr4$measdate <- as.Date(ymd(paste(sr4$MEASYEAR,sr4$MEASMON,sr4$MEASDAY, sep="-")))

# sort by plot and stem for lagging variables
setorderv(sr4, c('tr.tag', 'measdate'))

# make flag for remeasured tree
sr4$remp <- ifelse(dplyr::lag(sr4$tr.tag,n=1)==sr4$tr.tag,TRUE,FALSE)
sr4$remp[1] <- FALSE

# add census number  
cns <- rep(1, length=nrow(sr4))

for(j in 2:nrow(sr4)){
  cns[j] <- ifelse(sr4$remp[j], cns[j-1]+1, cns[j]) 
}

sr4$census <- cns

# subset out just columns needed
sr5 <- sr4[ , c("plotid", "tr.tag", "census", 
                "SPCD","TRE_STATUSCD", "DIA", 
                "MEASYEAR", "PREVDIA", 
                "PREV_STATUS_CD", "measdate", "AGENTCD",
                "MORTYR", "ECOSUBCD",
                "LAT", "LON")]

# unit conversion
sr5$dbh_0 <- round(sr5$DIA*25.4, 0)

# set initial measurement
sr5$year_0 <- sr5$MEASYEAR
sr5 <- data.table(sr5)
# lead values to get interval data
sr5[, year_1:= shift(year_0, -1), by=tr.tag]
sr5[, dbh_1:= shift(dbh_0, -1), by = tr.tag]
sr5[, stat_1:= shift(TRE_STATUSCD, -1), by = tr.tag]
sr5[, date_1:= shift(measdate, -1), by = tr.tag]
sr5[, agent1:= shift(AGENTCD, -1), by = tr.tag]
sr5[, mortyr:= shift(MORTYR, -1), by = tr.tag]
sr5$nbdays <- as.numeric(sr5$date_1-sr5$measdate)

# lots of standing dead remeasured and need removed
sr5$surv <- ifelse(sr5$TRE_STATUSCD==1 & sr5$stat_1==1, 1,
                   ifelse(sr5$TRE_STATUSCD==1 & sr5$stat_1==2, 0,
                          ifelse(sr5$TRE_STATUSCD==2 & sr5$stat_1==2, NA, NA)))

sr5 <- subset(sr5, !is.na(sr5$surv))

sr5$sp_name <- as.character(sr5$SPCD)

# clean up final data
df.stan <- data.table(plot=sr5$plotid, long=sr5$LON, lat=sr5$LAT, dbh=sr5$dbh_0, 
                      stem=sr5$tr.tag, Sp=sr5$sp_name, agent = sr5$agent1,
                      mortyr = sr5$mortyr,
                      year_0= sr5$year_0,
                      year_1=sr5$year_1, nbdays=sr5$nbdays, surv=sr5$surv)

# remove years with small sample sizes
df.stan <- subset(df.stan, df.stan$year_1 < 2022 & df.stan$year_0 > 2002)

# look at number of stems by species
#df.stan[,.(count = length(unique(stem))), by = Sp]

# add grouping cells (50 km or 100 km on a side)
dat.u <-  unique(df.stan[,c('long','lat')]) # get unique plots not trees
dat_sf <- st_as_sf(dat.u, coords = c('long', 'lat'), crs = 4326) # define projection
dat_tr <- st_transform(dat_sf, crs  = 32616)  # transform to UTM

# set cell size in meters (50 x 50 km grain, for the manuscript main analyses):
# ____________________________________________________________________________
cellsz <- c(50000, 50000)

#make grid
area_hex_grid_50 <- st_make_grid(dat_tr, cellsize = cellsz,  what = 'polygons', square = TRUE)

#add cell id
honey_50 <- st_sf(area_hex_grid_50) %>%
  mutate(cell_50 = 1:length(lengths(area_hex_grid_50)))
#saveRDS(honey_50, file = "./Output/surv_grid_50.RDS")

honey_50$n_plot <- lengths(st_intersects(honey_50, dat_tr))

# retrieve center coordinates of each grid cell:
honey_50 <- honey_50 %>% 
  bind_cols(st_coordinates(st_centroid(area_hex_grid_50))) %>% 
  rename(long_50 = X, lat_50 = Y)

# get rid of empty cells
honey_50 <- subset(honey_50, honey_50$n_plot > 0)
#summary(honey_50$n_plot)

#remove plot count
honey_50$n_plot<-NULL

# get original data in correct format for spatial join
dat_trees <- st_as_sf(df.stan, coords = c('long', 'lat'), crs = 4326) # define projection
dat_trees <- st_transform(dat_trees, crs  = 32616)  # transform to UTM

#join cell id and tree data
dat_cell_50 <- st_join(dat_trees, honey_50, left = TRUE)

# do it again for 100 x 100 km grain   #
# ______________________
cellsz <- c(100000, 100000)

area_hex_grid_100 <- st_make_grid(dat_tr, cellsize = cellsz,  what = 'polygons', square = TRUE)
#saveRDS(area_hex_grid_100, file = "./Output/surv_grid_100.RDS")

#add cell id
honey_100 <- st_sf(area_hex_grid_100) %>%
  mutate(cell_100 = 1:length(lengths(area_hex_grid_100)))
#saveRDS(honey_100, file = "./Output/surv_grid_100.RDS")

honey_100$n_plot <- lengths(st_intersects(honey_100, dat_tr))
# retrieve center coordinates of each grid cell:
honey_100 <- honey_100 %>% 
  bind_cols(st_coordinates(st_centroid(area_hex_grid_100))) %>% 
  rename(long_100 = X, lat_100 = Y)

# get rid of empty cells
honey_100 <- subset(honey_100, honey_100$n_plot > 0)
#summary(honey_100$n_plot)

honey_100$n_plot<-NULL

#join tree data with cells
all_dat <- st_join(dat_cell_50, honey_100, left = TRUE)

#remove sf geometry
all_dat2 <- st_drop_geometry(all_dat)
all_dat2$yrstart <- all_dat2$year_0 - min(all_dat2$year_0) + 1
all_dat2$yrend <- all_dat2$year_1 - min(all_dat2$year_0) + 1


# Add column for mortality agent as character string (see metadata .docx file):
all_dat2$m.agent <- NA
all_dat2$m.agent[which(all_dat2$agent == 10)] <- "insect"
all_dat2$m.agent[which(all_dat2$agent == 20)] <- "disease"
all_dat2$m.agent[which(all_dat2$agent == 30)] <- "fire"
all_dat2$m.agent[which(all_dat2$agent == 40)] <- "animal"
all_dat2$m.agent[which(all_dat2$agent == 50)] <- "weather"
all_dat2$m.agent[which(all_dat2$agent == 60)] <- "vegetation"
all_dat2$m.agent[which(all_dat2$agent == 70)] <- "unknown_or_mult"
all_dat2$m.agent[which(all_dat2$agent == 80)] <- "human"


# cells of 50 km:
data2 <- all_dat2 %>% 
  select(plot:surv, 
         plot2 = cell_50, long2 = long_50, lat2 = lat_50,
         yrstart, yrend, agent, m.agent) %>% 
  mutate(plot2 = paste0("c_", plot2), 
         grain = "50km")
# Save:
# ____
fwrite(data2, "DATA/data2_50km_2023.csv")
# data2 <- fread("DATA/data2_50km.csv") 


# cells of 100 km:
data100 <- all_dat2 %>% 
  select(plot:surv, 
         plot2 = cell_100, long2 = long_100, lat2 = lat_100,
         yrstart, yrend, agent, m.agent) %>% 
  mutate(plot2 = paste0("c_", plot2), 
         grain = "100km")
# Save:
# ____
fwrite(data100, "DATA/data2_100km_2023.csv")


# Bayesian size-dependent multilevel model of survival ----
# ****************************************************

## Load data ----
# __________
data2 <- fread("DATA/data2_50km.csv") 

## Write the STAN survival model ----
# ______________________________

# Non-centered parametrization of the model, for faster sampling.

stan_surv_model <- "
data {
  int<lower=0> n;                   // no. of observations
  int<lower=0> n_plot2;             // no. of plots (rather cells, see Methods in manuscript)
  int<lower=0> T;                   // no. of Years
  real dbh[n];                      // diameter in mm
  int surv[n];                      // survival of each tree, 1 for alive, 0 for dead
  real thresh;                      // dbh (in mm) threshold separating first from second curve
  int cstart[n];                    // index for census first year
  int cend[n];                      // index for census last year
  int plot2[n];                     //  plot (cell) identifier
  int Tind[T];                      //  year identifier
}

parameters {
  // Main model parameters
  real K_mu;

  // Non-centered priors (see transformed parameters block):
  vector[T] K_T_z;
  vector[n_plot2] K_P_z;

  // varying effect parameters
  real<lower=0> sigma_K_T;
  real<lower=0> sigma_K_P;

  // These are stricter priors that better resolve the identity of the p and r  
  // If dbh is untransformed:
  real<lower=-1*thresh, upper=thresh> p1;  
  real<lower=0.0001, upper=0.25> r1; 
  real<lower=thresh, upper=3*max(dbh)> p2;  // 1.5 usually, but we test less constrain on p2
  real<lower=-0.45, upper=-0.001> r2;
}

transformed parameters {
  vector[T] K_T;
  vector[n_plot2] K_P;

  // Non-centered priors (smuggle the hyperparameters out of the adaptive priors):
  K_T = K_T_z * sigma_K_T;
  K_P = K_P_z * sigma_K_P;
}

model {
  real theta;
  real p;
  real K;
  real K_tmp;

  // Non-centered priors, to avoid having hyperparameters in the adaptive priors (more efficient MCMC sampling):
  K_T_z ~ normal(0, 1);
  K_P_z ~ normal(0, 1);
  K_mu ~ normal(2, 2);   // boot::inv.logit(rnorm(1000, 2, 2)) to explore implications

  // Priors on hyperparameters for varying effects are here
  sigma_K_T ~ gamma(3, 4);
  sigma_K_P ~ gamma(3, 4);

  for(i in 1:n) {
    if(dbh[i] < thresh) {
      for(t in cstart[i]:cend[i]) {
        K_tmp = K_mu + K_P[plot2[i]] + K_T[Tind[t]];  //  + eps[i]
        K = inv_logit(K_tmp);
        p = (K) / (1 + exp(-r1 * (dbh[i] - p1)));
        theta = p;
      }
    }
    else {
      for(t in cstart[i]:cend[i]) {
        K_tmp = K_mu + K_P[plot2[i]] + K_T[Tind[t]];  
        K = inv_logit(K_tmp);
        p = (K) / (1 + exp(-r2 * (dbh[i] - p2)));
        theta = p;
      }
    }
    surv[i] ~ bernoulli(theta);
  }
}
"

## Run the survival model ----
# _______________________

n_iter <- 2000
n_warmup <- 500
n_chains <- 4    # defines both the number of chains and cores used
seed <- 2022
n_samples <- 200  # nb of samples from the posterior draws for each observation (to generate K_tmp)

sp_name <- "111" # code of the focal species

sp.df <- subset(data2, Sp == sp_name)

# Run the following chunk once only to build and save the stan model.
# Build the survival model on one species, instead of building it every time in run.stan():
# ___________________________________________________________________
tmp <- sp.df
# tmp$dbh <- as.numeric(scale(tmp$dbh))

thresh <- mean(tmp$dbh) # 300, for 131 only (mean dbh for 111 and 121)
tot.years <- max(tmp$yrend) - min(tmp$yrstart) + 1
surv.data <- tidybayes::compose_data(tmp)

surv.data <- list_modify(.x = surv.data,
                         T = tot.years, # not great to use T as a variable name
                         thresh = thresh,
                         cstart = tmp$yrstart - min(tmp$yrstart) + 1,
                         cend = tmp$yrend - min(tmp$yrstart) + 1,
                         Tind = seq(tot.years))

surv.model <- stan(model_code = stan_surv_model,
                   data = surv.data, chains = 0, save_dso = TRUE, verbose = FALSE)
# surv.model <- stan(file = "CODE/MORTALITY_Year_Plot_nc_fishnet_no_eps_varying_K.stan",
#                    data = surv.data, chains = 0, save_dso = TRUE, verbose = FALSE)
# save(surv.model, file = "CODE/surv_model_Year_Plot_nc_fishnet_no_eps_varying_K.RData")

surv.fit <- rstan::stan(fit = surv.model, data = surv.data, 
                        iter = n_iter, warmup = n_warmup, chains = n_chains, cores = n_chains, 
                        verbose = TRUE, #init = 0, 
                        control = list(adapt_delta = 0.98, max_treedepth = 15),
                        seed = seed,
                        # save_warmup(FALSE), 
                        include = FALSE, pars = c("K_T_z", "K_P_z")) 

## Extract a parameters' posteriors ----
# _________________________________
# tidybayes::spread_draws() generates for each level of K_P for example (i.e. each plot), 
# a number of rows equal to the number of iterations (after removing the warmups), 
# i.e. n_iter/2 by default, multiplied by the n_chains. 
# tidybayes::recover_types() is a handy function that back-transforms the index variables
# used in Stan (i.e. plots and years) into their original values.

# K:
# __
extract_Kmu <- tidybayes::spread_draws(model = surv.fit,
                                       K_mu,  seed = seed, ndraws = n_samples)
extract_KP <- tidybayes::spread_draws(model = recover_types(surv.fit, sp.df),
                                        K_P[plot2], seed = seed, ndraws = n_samples) 
extract_KT <- tidybayes::spread_draws(model = surv.fit,
                                      K_T[T], seed = seed, ndraws = n_samples) 
extract_KT$T <- extract_KT$T + min(sp.df$year_0) - 1 # Convert 'T' back in corresponding years
names(extract_KT) <- c('year','K_T', ".chain", ".iteration", ".draw")

# r1 and r2:
# __________
extract_r1 <- tidybayes::spread_draws(model = recover_types(surv.fit, sp.df),
                                      r1, seed = seed, ndraws = n_samples) 
extract_r2 <- tidybayes::spread_draws(model = recover_types(surv.fit, sp.df),
                                      r2, seed = seed, ndraws = n_samples) 

# p1 and p2:
# __________
extract_p1 <- tidybayes::spread_draws(model = surv.fit, 
                                      p1,  seed = seed, ndraws = n_samples) 

extract_p2 <- tidybayes::spread_draws(model = surv.fit,
                                      p2,  seed = seed, ndraws = n_samples) 

sp.df <- data.table(sp.df)

cyear <- sp.df[, list(plot = plot,
                      long = long,
                      lat = lat,
                      plot2 = plot2,
                      long2 = long2,
                      lat2 = lat2,
                      grain = grain,
                      surv = surv,
                      dbh = dbh,
                      stem = as.character(stem),
                      year = seq(year_0, year_1)),
               by = 1:nrow(sp.df) ]

cyear <- left_join(cyear, extract_KT, by = 'year', relationship = "many-to-many")
cyear <- left_join(cyear, extract_Kmu, by = c('.chain', '.iteration', '.draw'))
cyear <- left_join(cyear, extract_KP) #  by = c('plot', '.chain', '.iteration', '.draw') # don't specify 'by' to avoid having to adjust 'plot' or 'plot2'

cyear <- left_join(cyear, extract_r1)
cyear <- left_join(cyear, extract_r2)

cyear <- left_join(cyear, extract_p1)
cyear <- left_join(cyear, extract_p2)

df_Ktmp <- data.table(species = rep(sp_name, sum(surv.data$cend - surv.data$cstart + 1) * n_samples), 
                      obs = cyear$nrow, 
                      stem = cyear$stem,
                      plot = cyear$plot, 
                      long = cyear$long,
                      lat = cyear$lat,
                      plot2 = cyear$plot2, 
                      long2 = cyear$long2,
                      lat2 = cyear$lat2,
                      year = cyear$year, 
                      dbh = cyear$dbh, 
                      surv = cyear$surv,
                      K_tmp = with(cyear, K_mu + K_T + K_P), 
                      K_mu = cyear$K_mu,
                      K_P = cyear$K_P,
                      K_T = cyear$K_T,
                      p1 = cyear$p1,
                      p2 = cyear$p2)

df_Ktmp$r1 <- cyear$r1
df_Ktmp$r2 <- cyear$r2

df_Ktmp[, theta := ifelse(dbh <= thresh, 
                          boot::inv.logit(K_tmp) / 
                            (1 + exp(-(r1 * ((dbh - p1))))^(1)),
                          boot::inv.logit(K_tmp) / 
                            (1 + exp(-(r2 * ((dbh - p2))))^(1)))]

df_Ktmp$seed <- seed

# Save:
suffix <- paste0("varying_K_", sp_name, "_", sp.df$grain[1], "_nsamples", n_samples)

# optional:    ### silence or leave
suffix <- paste0(suffix, "_threshmean")  # _thresh300 (for species 131) ### adapt here

### Save df_Ktmp ----
save(df_Ktmp, file = sprintf("OUTPUT/df_Ktmp_plot_year_%s.RData", suffix))  

# Create lighter version of df_Ktmp (keeping columns necessary to generate figures):
df_Ktmp <- dplyr::select(df_Ktmp,
                         -c(theta, seed, obs))
save(df_Ktmp, file = sprintf("OUTPUT/df_Ktmp_plot_year_%s_lesscolumns.RData", suffix))  

# All future analyses and figures will generally start from object 'df_Ktmp'.

## Fig. 1 and Fig. S4 - Visualize fitted survival curves ----
# *****************************************************

# Customized figure theme:
mytheme <- theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 15),
        plot.caption = element_text(size = 9))

# Helper functions:
# ________________
predict.surv <- function(x, params, time = 1) {
  K <- params[1]
  p <- params[2]
  r <- params[3]
  # Linear predictor parameters
  pred.y <-  K / (1 + exp(-(r * ((x - p) )))) ^ (time)
  return(pred.y)
}

get.curve.mat <- function(df, subsample = TRUE, s.no = 500, thresh = 0, 
                          min_dbh = -2, max_dbh = 2, suppl_dbh = 0, length_dbh = 1000,
                          dbh_vec = NULL,
                          eps = FALSE,
                          var_r = FALSE
) {
  # get.curve() from 'df_Ktmp'
  
  # Function to generate survival probabilities on raw scale (theta, [0, 1]) for a sequence of 1000 dbh values ranging from
  # min_dbh to max_dbh + suppl_dbh. The function generates survival predictions based on either all posterior draws (subsample = FALSE)
  # or a subset of 's.no' posterior draws (subsample = TRUE).
  
  # - df: df_Ktmp or a subset of df_Ktmp (posterior draws for all survival parameters (see output from run.stan function))
  # - subsample, s.no, min_dbh, max-dbh, suppl_dbh: see function explanation above.
  # - eps: logical; whether epsilon was used to define K_tmp or not, in stan model. 
  # - var_r: whether the survival model was fitted with varying intervals per spatial unit (e.g. plot) for r1 and r2 (r_mu + r_P).
  # - dbh_vec: vector of dbhs of interest; to use instead of min_dbh, max_dbh, suppl_dbh and length_dbh if one is interested in particular
  # values of dbh, instead of regularly spaced vector of dbh between a min and a max.
  # - length_dbh: length of the simulated dbh vector going from min_dbh to max_dbh; default is 1000.
  
  if (subsample == TRUE) {
    samp.id <- sample(seq(nrow(df)), s.no)
  } else samp.id <- seq(nrow(df)) # all rows
  
  if (eps == FALSE) {
    K  <- with(df, boot::inv.logit(K_mu[samp.id] + K_P[samp.id] + K_T[samp.id]))
  } else K  <- with(df, boot::inv.logit(K_mu[samp.id] + K_P[samp.id] + K_T[samp.id] + eps[samp.id]))
  p1 <- df$p1[samp.id]
  p2 <- df$p2[samp.id]
  if (var_r == FALSE) {
    r1 <- df$r1[samp.id]
    r2 <- df$r2[samp.id]
  } else {
    r1 <- with(df, r1_mu[samp.id] + r1_P[samp.id])
    r2 <- with(df, r2_mu[samp.id] + r2_P[samp.id])
  }
  
  params.lo <- split(cbind(K, p1, r1), seq(length(K))) # combinations of 3 param.
  params.hi <- split(cbind(K, p2, r2), seq(length(K))) # combinations of 3 param.
  if(is.null(dbh_vec)) {
    dbh <- seq(from = min_dbh, to = max_dbh + suppl_dbh, length = length_dbh)
  } else dbh <- dbh_vec
  res.mat <- as.data.frame(matrix(0, ncol = length(dbh), nrow = length(K)))
  colnames(res.mat) <- as.character(round(dbh, 0))
  
  for(i in 1:nrow(res.mat)) {
    lo.prob <- predict.surv(x = dbh[dbh < thresh], params = params.lo[[i]])
    hi.prob <- predict.surv(x = dbh[dbh >= thresh], params = params.hi[[i]])
    res.mat[i, ] <- c(lo.prob, hi.prob)
  }
  return(res.mat)
}

make.surv.fig <- function(res.mat, min_dbh = -2, max_dbh = 2, 
                          scale = FALSE, mean_dbh, sd_dbh,
                          suppl_dbh = 1, species = "", 
                          probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                          y_lower = 0,
                          theme = NULL) {
  
  # - res.mat: matrix of surv. prob. predicted values (rows) per dbh (columns), for 1000 dbh values (columns) ranging from min_dbh to max_dbh
  # + suppl_dbh (as generated by function get.curv.mat()).
  # - scale: TRUE if survival model was run on scaled dbh (i.e. centered and scaled); default is FALSE (raw dbh). If true, mean_dbh and sd_dbh
  # are used to back-transform dbhs to raw scale.
  # - probs: vector defining the lower and upper bounds of the credibility intervals to represent on figure; first and last values
  # together, second and second to last together. 0.5 is there to obtain the median -> line defining the most likely curve.
  # - y_lower: allows setting a lower bound to the y-axis of survival probabilities; default is 0.
  # - theme: define a personalised theme for the figure; default is theme_classic()
  
  if (is.null(theme)) theme <- theme_classic()
  
  dbh <- seq(from = min_dbh, to = max_dbh + suppl_dbh, length = 1000)
  quants <- apply(res.mat, 2, quantile, probs = probs)
  rownames(quants) <- paste("q", str_replace(rownames(quants), "%", ""), sep = "")
  quants_t <- as.data.frame(cbind(t(quants), dbh))
  
  if (scale == TRUE) { 
    # mean_dbh and sd_dbh are used to back-transform the dbh to original scale:
    quants_t$dbh <- (quants_t$dbh * sd_dbh) + mean_dbh
    max <- (max_dbh * sd_dbh) + mean_dbh
  } else max <- max_dbh
  
  ggplot(data = quants_t, aes(x = dbh, y = q50)) +
    geom_line(colour = "red3", linewidth = 0.8) +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.3) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.5) +
    geom_vline(xintercept = max, linetype = 2) +
    coord_cartesian(ylim = c(y_lower, 1)) +
    labs(x = "DBH (mm)",
         y = "Survival probability",
         subtitle = paste0("Species: ", sp_name),
         caption = "Median, 50%- and 90%-credibility intervals of posterior surv. probabilities") +
    theme
}

# Functions to visualise fitted survival model as a function of DBH:
source("CODE/surv_plotting_functions.R")

### Load data at the spatial grain used for the analysis of interest ----
# _________________________________________________________________

# For Fig. 1:
data <- fread("DATA/data2_50km_2023.csv") # what sp.df was based on

# Or for supplementary figure, Fig. S4:
data <- fread("DATA/data2_100km_2023.csv") # what sp.df was based on

### Generate posterior predictions for the figure ----

#### Species 111 ----
# *************
# *************

# For Fig. 1:
load("OUTPUT/df_Ktmp_plot_year_varying_K_111_50km_threshmean_lesscolumns.RData") # sp 111

# Or for supplementary figure Fig. S4 (grain of 100 km):
load("OUTPUT/df_Ktmp_plot_year_varying_K_111_100km_threshmean_lesscolumns.RData") # sp 111

sp_name <- df_Ktmp$species[1]

dfsub <- subset(data, Sp == sp_name)

# Extract posterior draws:
# _______________________
# samps <- rstan::extract(surv.fit)

min_dbh <- min(dfsub$dbh)
max_dbh <- max(dfsub$dbh) # species 111: 780; species 121: 714; species 131: 1024
mean_dbh <- mean(dfsub$dbh)

# We define an extra dbh beyond the observed max, until which we want to observe predictions of surv. prob.:
suppl_dbh <- 0 # 400 for sp 121, and 600 for sp 111

# or set to fixed value, to compare different species on same scale:
x_upper_limit <- max_dbh # 1024 # 1024 is the max_dbh of species 131. Choosing this value allows to have the same x-axis for all species.
suppl_dbh <- x_upper_limit - max_dbh 

ndraws <- 1000

res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = mean_dbh,  # 300     ### adapt here
                         min_dbh = min_dbh, # (min_dbh - mean_dbh) / sd_dbh
                         max_dbh = max_dbh, # ((max_dbh - mean_dbh) / sd_dbh)
                         suppl_dbh = suppl_dbh) # in mm

fig111_50km <-      ### adapt name here (from the species code of the considered species) 
  make.surv.fig(res.mat,
                min_dbh = min_dbh,
                max_dbh = max_dbh,
                scale = FALSE,
                # mean.dbh = mean_dbh,
                # sd.dbh = 1, # sd_dbh
                suppl_dbh = suppl_dbh, # sd(dfsub$dbh)
                species = sp_name,
                y_lower = 0,
                theme = mytheme)

# Save:
saveRDS(fig111_50km, file = "OUTPUT/figures/figures_manuscript/Fig. 1/fig111_50km_surv_curve.rds")

saveRDS(fig111_100km, file = "OUTPUT/figures/figures_supp_manuscript/Fig_S4/fig111_100km_surv_curve.rds")

# Survival posterior at a given size:
sizes <- c(25, 100, 500) # dbh in mm 
res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = mean_dbh,  # 300     ### adapt here
                         dbh_vec = sizes) # in mm

apply(res.mat, 2, median)
apply(res.mat, 2, function(x) qi(x, .width = 0.9))
# apply(res.mat, 2, function(x) hdi(x, credMass = c(0.8, 0.95)))


#### Species 121 ----
# *************

load("OUTPUT/df_Ktmp_plot_year_varying_K_121_50km_threshmean_lesscolumns.RData") # sp 121

load("OUTPUT/df_Ktmp_plot_year_varying_K_121_100km_threshmean_lesscolumns.RData") # sp 121

sp_name <- df_Ktmp$species[1]

dfsub <- subset(data, Sp == sp_name)

# Extract posterior draws:
# ________________________
# samps <- rstan::extract(surv.fit)

min_dbh <- min(dfsub$dbh)
max_dbh <- max(dfsub$dbh) # species 111: 780; species 121: 714; species 131: 1024
mean_dbh <- mean(dfsub$dbh)

# We define an extra dbh beyond the observed max, until which we want to observe predictions of surv. prob.:
suppl_dbh <- 0 # 400 for sp 121, and 600 for sp 111

# or set to fixed value, to compare different species on same scale:
x_upper_limit <- max_dbh # 1024 # 1024 is the max_dbh of species 131. Choosing this value allows to have the same x-axis for all species.
suppl_dbh <- x_upper_limit - max_dbh 

ndraws <- 1000

res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = mean_dbh,     ### adapt here
                         min_dbh = min_dbh, # (min_dbh - mean_dbh) / sd_dbh
                         max_dbh = max_dbh, # ((max_dbh - mean_dbh) / sd_dbh)
                         suppl_dbh = suppl_dbh) # in mm

fig121_50km <-      ### adapt name here
  make.surv.fig(res.mat,
                min_dbh = min_dbh,
                max_dbh = max_dbh,
                scale = FALSE,
                # mean.dbh = mean_dbh,
                # sd.dbh = 1, # sd_dbh
                suppl_dbh = suppl_dbh, # sd(dfsub$dbh)
                species = sp_name,
                y_lower = 0,
                theme = mytheme)

# Save:
saveRDS(fig121_50km, file = "OUTPUT/figures/figures_manuscript/Fig. 1/fig121_50km_surv_curve.rds")

saveRDS(fig121_100km, file = "OUTPUT/figures/figures_supp_manuscript/Fig_S4/fig121_100km_surv_curve.rds")

# Survival posterior at a given size:
sizes <- c(25, 100, 500) # dbh in mm
res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = mean_dbh,  # 300     ### adapt here
                         dbh_vec = sizes) # in mm

apply(res.mat, 2, median)
apply(res.mat, 2, function(x) qi(x, .width = 0.9))
# apply(res.mat, 2, function(x) hdi(x, credMass = c(0.8, 0.95)))


#### Species 131 ----
# *************

load("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData")  # sp 131

load("OUTPUT/df_Ktmp_plot_year_varying_K_131_100km_thresh300_lesscolumns.RData")  # sp 131

sp_name <- df_Ktmp$species[1]

dfsub <- subset(data, Sp == sp_name)

# Extract posterior draws:
# ________________________
# samps <- rstan::extract(surv.fit)

min_dbh <- min(dfsub$dbh)
max_dbh <- max(dfsub$dbh) # species 111: 780; species 121: 714; species 131: 1024
mean_dbh <- mean(dfsub$dbh)

# We define an extra dbh beyond the observed max, until which we want to observe predictions of surv. prob.:
suppl_dbh <- 0 # 400 for sp 121, and 600 for sp 111

# or set to fixed value, to compare different species on same scale:
x_upper_limit <- max_dbh # 1024 # 1024 is the max_dbh of species 131. Choosing this value allows to have the same x-axis for all species.
suppl_dbh <- x_upper_limit - max_dbh 

ndraws <- 1000

res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = 300,  # mean_dbh     ### adapt here
                         min_dbh = min_dbh, # (min_dbh - mean_dbh) / sd_dbh
                         max_dbh = max_dbh, # ((max_dbh - mean_dbh) / sd_dbh)
                         suppl_dbh = suppl_dbh) # in mm

fig131_100km <-      ### adapt name here
  make.surv.fig(res.mat,
                min_dbh = min_dbh,
                max_dbh = max_dbh,
                scale = FALSE,
                # mean.dbh = mean_dbh,
                # sd.dbh = 1, # sd_dbh
                suppl_dbh = suppl_dbh, # sd(dfsub$dbh)
                species = sp_name,
                y_lower = 0,
                theme = mytheme)

saveRDS(fig131_50km, file = "OUTPUT/figures/figures_manuscript/Fig. 1/fig131_50km_surv_curve.rds")

saveRDS(fig131_100km, file = "OUTPUT/figures/figures_supp_manuscript/Fig_S4/fig131_100km_surv_curve.rds")

# Survival posterior at a given size:
sizes <- c(25, 100, 500) 
res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = 300,  # 300     ### adapt here
                         dbh_vec = sizes) # in mm

apply(res.mat, 2, median)
apply(res.mat, 2, function(x) qi(x, .width = 0.9))
# apply(res.mat, 2, function(x) hdi(x, credMass = c(0.8, 0.95)))


### Figure of fitted size-dependent survival curve ----
# ***********************************************
# ***********************************************

#### Fig. 1 (grain of 50 km) ----
# _________
fig111_50km <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 1/fig111_50km_surv_curve.rds")
fig121_50km <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 1/fig121_50km_surv_curve.rds")
fig131_50km <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 1/fig131_50km_surv_curve.rds")

fig1 <- (fig111_50km + 
           labs(subtitle = expression(italic("Pinus elliottii"))) +
           coord_cartesian(ylim = c(0.45, 1)) + 
           scale_y_continuous(breaks = seq(0.5, 1, 0.1)) +
           scale_x_continuous(limits = c(0, 1024)) +
           theme(plot.caption = element_blank())) / (fig121_50km + 
                                                       labs(subtitle = expression(italic("Pinus palustris"))) +
                                                       coord_cartesian(ylim = c(0.45, 1)) + 
                                                       scale_y_continuous(breaks = seq(0.5, 1, 0.1)) +
                                                       scale_x_continuous(limits = c(0, 1024)) +
                                                       theme(plot.caption = element_blank())) / (fig131_50km + 
                                                                                                   labs(subtitle = expression(italic("Pinus taeda"))) +
                                                                                                   coord_cartesian(ylim = c(0.45, 1)) + 
                                                                                                   scale_y_continuous(breaks = seq(0.5, 1, 0.1)) +
                                                                                                   scale_x_continuous(limits = c(0, 1024)) + 
                                                                                                   theme(plot.caption = element_blank())) +
  
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 15, face = "bold"))

fig1

# Save Fig. 1:
# ___________
ggsave(plot = fig1, 
       filename = "OUTPUT/figures/figures_manuscript/Fig. 1/Fig1.pdf", width = 8, height = 8, units = "in")

ggsave(plot = fig1, 
       filename = "OUTPUT/figures/figures_manuscript/Fig. 1/Fig1.png", width = 8, height = 8, units = "in", dpi = 400)

#### Fig. S4 (grain of 100 km) ----
# _________
fig111_100km <- readRDS("OUTPUT/figures/figures_supp_manuscript/Fig_S4/fig111_100km_surv_curve.rds")
fig121_100km <- readRDS("OUTPUT/figures/figures_supp_manuscript/Fig_S4/fig121_100km_surv_curve.rds")
fig131_100km <- readRDS("OUTPUT/figures/figures_supp_manuscript/Fig_S4/fig131_100km_surv_curve.rds")

figs4 <- (fig111_100km + 
            labs(subtitle = expression(italic("Pinus elliottii"))) +
            coord_cartesian(ylim = c(0.45, 1)) + 
            scale_y_continuous(breaks = seq(0.5, 1, 0.1)) +
            scale_x_continuous(limits = c(0, 1024)) +
            theme(plot.caption = element_blank())) / (fig121_100km + 
                                                        labs(subtitle = expression(italic("Pinus palustris"))) +
                                                        coord_cartesian(ylim = c(0.45, 1)) + 
                                                        scale_y_continuous(breaks = seq(0.5, 1, 0.1)) +
                                                        scale_x_continuous(limits = c(0, 1024)) +
                                                        theme(plot.caption = element_blank())) / (fig131_100km + 
                                                                                                    labs(subtitle = expression(italic("Pinus taeda"))) +
                                                                                                    coord_cartesian(ylim = c(0.45, 1)) + 
                                                                                                    scale_y_continuous(breaks = seq(0.5, 1, 0.1)) +
                                                                                                    scale_x_continuous(limits = c(0, 1024)) + 
                                                                                                    theme(plot.caption = element_blank())) +
  
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 15, face = "bold"))

figs4

ggsave(plot = figs4, 
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S4/FigS4.pdf", width = 8, height = 8, units = "in")

ggsave(plot = figs4, 
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S4/FigS4.png", width = 8, height = 8, units = "in", dpi = 400)

## Fig. S2 - Interspecific comparisons of survival curves (contrasts) ----
# *******************************************************************

# If not already loaded/created, create functions predict.surv, get.curv.mat, make.surv.fig,
# introduced above for Fig. 1.

# Theme for figure:
mytheme <- theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        plot.title = element_text(size = 16, colour = "black"),
        plot.subtitle = element_text(size = 15, colour = "black"),
        plot.caption = element_text(size = 11, colour = "black"))

# Load data used to run the model and stan model output:
# ______________________________________________________

data <- fread("DATA/data2_50km_2023.csv") # what sp.df was based on

# Define the range of DBHs for which we'll want predicted survival to compare the species two by two:
min_dbh <- 25
max_dbh <- 700

ndraws <- 1000 # nb of posterior draws to use

### Species 111: P. elliottii ----
# ___________________________

load("OUTPUT/df_Ktmp_plot_year_varying_K_111_50km_threshmean_lesscolumns.RData") # sp 111

sp_name <- df_Ktmp$species[1]

dfsub <- subset(data, Sp == sp_name)
mean_dbh <- mean(dfsub$dbh)

res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = mean_dbh,  # 300     ### adapt here
                         min_dbh = min_dbh, # (min_dbh - mean_dbh) / sd_dbh
                         max_dbh = max_dbh,
                         length_dbh = 1000) # in mm

# res.mat generates 1000 equally distant dbhs for the predictions
dbh_vec <- seq(min_dbh, max_dbh, length = 1000)

# transform in long format:
posterior_111 <- res.mat %>% 
  pivot_longer(cols = everything(), names_to = "class", values_to = "surv") %>% 
  mutate(dbh = rep(dbh_vec, times = ndraws),
         sp = "P. elliottii") %>% 
  select(sp, dbh, surv)

# make.surv.fig(res.mat,
#                 min_dbh = min_dbh,
#                 max_dbh = max_dbh,
#                 scale = FALSE,
#                 species = sp_name,
#                 y_lower = 0,
#                 theme = mytheme)

### Species 121: P. palustris ----
# __________________________

load("OUTPUT/df_Ktmp_plot_year_varying_K_121_50km_threshmean_lesscolumns.RData") # sp 121

sp_name <- df_Ktmp$species[1]

dfsub <- subset(data, Sp == sp_name)
mean_dbh <- mean(dfsub$dbh)

res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = mean_dbh,  # 300     ### adapt here
                         min_dbh = min_dbh, # (min_dbh - mean_dbh) / sd_dbh
                         max_dbh = max_dbh,
                         length_dbh = 1000) # in mm
# res.mat generates 1000 equally distant dbhs for the predictions
dbh_vec <- seq(min_dbh, max_dbh, length = 1000)

# transform in long format:
posterior_121 <- res.mat %>% 
  pivot_longer(cols = everything(), names_to = "class", values_to = "surv") %>% 
  mutate(dbh = rep(dbh_vec, times = ndraws),
         sp = "P. palustris") %>% 
  select(sp, dbh, surv)

# make.surv.fig(res.mat,
#                 min_dbh = min_dbh,
#                 max_dbh = max_dbh,
#                 scale = FALSE,
#                 species = sp_name,
#                 y_lower = 0,
#                 theme = mytheme)

save(posterior_111, posterior_121, file = "OUTPUT/posteriors_sp_111_121_for_contrasts_size_surv_traj.RData")

### Species 131: P. taeda ----
# _______________________

load("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData") # sp 121

sp_name <- df_Ktmp$species[1]

dfsub <- subset(data, Sp == sp_name)
mean_dbh <- mean(dfsub$dbh)

res.mat <- get.curve.mat(df = df_Ktmp, 
                         subsample = TRUE,
                         eps = FALSE,  ### adapt here
                         var_r = FALSE, ### adapt here
                         s.no = ndraws, # nb of posterior samples to take for the figure (used if subsample = TRUE)
                         thresh = mean_dbh,  # 300     ### adapt here
                         min_dbh = min_dbh, # (min_dbh - mean_dbh) / sd_dbh
                         max_dbh = max_dbh,
                         length_dbh = 1000) # in mm
# res.mat generates 1000 equally distant dbhs for the predictions
dbh_vec <- seq(min_dbh, max_dbh, length = 1000)

# transform in long format:
posterior_131 <- res.mat %>% 
  pivot_longer(cols = everything(), names_to = "class", values_to = "surv") %>% 
  mutate(dbh = rep(dbh_vec, times = ndraws),
         sp = "P. taeda") %>% 
  select(sp, dbh, surv)

# make.surv.fig(res.mat,
#                 min_dbh = min_dbh,
#                 max_dbh = max_dbh,
#                 scale = FALSE,
#                 species = sp_name,
#                 y_lower = 0,
#                 theme = mytheme)

rm(df_Ktmp)

load("OUTPUT/posteriors_sp_111_121_for_contrasts_size_surv_traj.RData")

### Generate the posteriors of the contrasts ----
# ******************************************
# posterior_all <- bind_rows(posterior_111, posterior_121, posterior_131) # To overlay curves of the 3 species (aes(colour), for example)

contrasts <- posterior_131 %>% 
  select(-c(sp, surv)) %>% 
  mutate(cont_131_111 = posterior_131$surv - posterior_111$surv,
         cont_131_121 = posterior_131$surv - posterior_121$surv,
         cont_111_121 = posterior_111$surv - posterior_121$surv)

# Save:
save(contrasts, file = "OUTPUT/contrasts_for_figure_comparison_species.RData")
load("OUTPUT/contrasts_for_figure_comparison_species.RData")

# Summarise contrasts for the figure:

# contrasts_summ <- contrasts %>% 
#   group_by(dbh) %>% 
#   point_interval(.point = median, .interval = qi, .width = 0.9)

contrasts_summ <- contrasts %>% 
  pivot_longer(cols = starts_with("cont"), names_to = "sp_pair", values_to = "contrast") %>% 
  group_by(dbh, sp_pair) %>% 
  point_interval(.point = median, .interval = qi, .width = 0.9)

### Create Fig. S2 ----
# ________________

mypalette <- c("#f499ff", "#006431", "#c40067")

contrasts_summ %>% 
  ggplot(aes(dbh, contrast)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = sp_pair), alpha = 0.4) +
  geom_line(aes(colour = sp_pair), size = 1) +
  scale_fill_manual(values = mypalette, labels = c(expression(italic(P.~elliottii)~-~italic(P.~palustris)),
                                                   expression(italic(P.~taeda)~-~italic(P.~elliottii)~~~~~~""),
                                                   expression(italic(P.~taeda)~-~italic(P.~palustris)~~~"")),
                    name = "Contrasts") +
  scale_colour_manual(values = mypalette, labels = c(expression(italic(P.~elliottii)~-~italic(P.~palustris)),
                                                     expression(italic(P.~taeda)~-~italic(P.~elliottii)~~~~~~""),
                                                     expression(italic(P.~taeda)~-~italic(P.~palustris)~~~"")),
                      name = "Contrasts") +
  geom_hline(yintercept = 0, lty = 2, colour = "black", size = 1) +
  scale_x_continuous(breaks = seq(0, 700, 100)) +
  labs(x = "DBH (mm)", y = "Interspecific survival prob. difference") +
  mytheme +
  theme(legend.position = c(0.85, 0.12))

# ggsave(...)


## Fig. S1 - Schematic illustration of size-dependent survival model and fitted curve ----
# ***********************************************************************************

# Theme for figures:
mytheme <- theme_classic() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"), # 12
        axis.title = element_text(size = 12, colour = "black"), # 14
        legend.text = element_text(size = 10, colour = "black"), # 12
        legend.title = element_text(size = 12, colour = "black"), # 14
        plot.title = element_text(size = 16, colour = "black"),
        plot.subtitle = element_text(size = 15, colour = "black"),
        plot.caption = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"))

# Helper functions:
predict.surv <- function(x, params, time = 1) {
  K <- params[1]
  p <- params[2]
  r <- params[3]
  # Linear predictor parameters
  pred.y <-  K / (1 + exp(-(r * ((x - p) )))) ^ (time)
  return(pred.y)
}

make.pred.sim <- function(K = 0.99, p1 = 20, r1 = 0.05, p2 = 1800, r2 = -0.02, 
                          thresh = 100, dbh.min = 10, dbh.max = 2000) {
  params.lo <- c(K, p1, r1)
  params.hi <- c(K, p2, r2)
  dbh <- seq(from = dbh.min, to = dbh.max, length = 1000)
  lo.prob <- predict.surv(x = dbh[dbh < thresh], params = params.lo)
  hi.prob <- predict.surv(x = dbh[dbh >= thresh], params = params.hi)
  pred <- c(lo.prob, hi.prob)
  df_fig <- as.data.frame(cbind(as.matrix(pred), dbh))
  colnames(df_fig) <- c("pred", "dbh")
  return(df_fig)
}

# Figure with one plot and one year, and one K_L for this plot and year:
# _____________________________________________________________________

# Plot 1: 
K <- 0.98
p1 <- 20
r1 <- 0.05
p2 <- 1800
r2 <- -0.02
thresh <- 270

df_fig1 <- make.pred.sim(K, p1, r1, p2, r2, thresh)
df_fig1 <- cbind(df_fig1, K = rep(K, 1000))

# Year 1: 
K <- 0.86
df_fig2 <- make.pred.sim(K, p1, r1, p2, r2, thresh)
df_fig2 <- cbind(df_fig2, K = rep(K, 1000))

# Average (inv.logit(K_mu)): 
# Doesn't need to be the mean(plots and years above) (supposes unrepresented groups)
K <- 0.95  
df_fig_u <- make.pred.sim(K, p1, r1, p2, r2, thresh)

# K_L: # mean of the K chosen above
(K <- mean(c(0.95, 0.98, 0.86)))
boot::inv.logit(boot::logit(0.95) + boot::logit(0.03) - boot::logit(0.09))
df_fig_KL <- make.pred.sim(K, p1, r1, p2, r2, thresh)
df_fig_KL <- cbind(df_fig_KL, K = rep(K, 1000))

df_fig_groups <- cbind(rbind(df_fig1, 
                             df_fig2,
                             df_fig_KL), 
                       group = rep(c("plot1", "year1", "KL"), each = 1000),
                       type = c(rep("plot", 1000), rep("year", 1000), rep("KL", 1000)))

# Main figure:
(fig <- ggplot(data = df_fig_groups, aes(x = dbh, y = pred)) +
    geom_line(aes(lty = type, group = group, colour = K), size = 0.7) +
    scale_linetype_manual(values = c(6, 2, 3)) +
    scale_colour_gradient(expression(paste(italic(K))), low = "darkred", high = "blue",
                          limits = c(0.86, 0.98), breaks = seq(0.86, 0.98, 0.04)) +
    guides(linetype = "none") +
    geom_line(data = df_fig_u, aes(dbh, pred), size = 1) +
    geom_vline(xintercept = thresh, colour = "grey60") +
    lims(y = c(0, 1),
         x = c(0, 2000)) +
    labs(x = "DBH (mm)",
         y = expression(paste(Annual~survival~probability~(italic(theta))))) +
    theme_classic() +
    mytheme)

# Main figure:
(fig_1 <- ggplot(data = df_fig_groups, aes(x = dbh, y = pred)) +
    geom_line(aes(lty = type, group = group, colour = K), size = 0.7) +
    scale_linetype_manual(values = c(6, 2, 3)) +
    scale_colour_gradient(expression(paste(italic(K))), low = "darkred", high = "blue",
                          limits = c(0.86, 0.98), breaks = seq(0.86, 0.98, 0.04)) +
    guides(linetype = "none") +
    geom_line(data = df_fig_u, aes(dbh, pred), size = 1) +
    # DBH threshold vertical line:
    geom_vline(xintercept = thresh, colour = "grey60") +
    annotate(geom = "text", x = 350, y = 0.03, label = expression(italic(threshold))) +
    # Highlight p1 and p2:
    geom_point(x = p1, y = 0.5, size = 2.5) +
    annotate(geom = "text", x = p1 + 50, y = 0.5, label = expression(italic(p1))) +
    geom_point(x = p2, y = 0.5, size = 2.5) +
    annotate(geom = "text", x = p2 + 50, y = 0.5, label = expression(italic(p2))) +
    # Highlight r1 and r2:
    annotate(geom = "text", x = 80, y = 0.64, label = expression(italic(r1))) +
    geom_segment(x = 100, y = 0.61, xend = 140, yend = 0.68, linewidth = 0.3) +
    geom_segment(x = 100, y = 0.61, xend = 140, yend = 0.61, linewidth = 0.3) +
    geom_segment(x = 140, y = 0.61, xend = 140, yend = 0.68, linewidth = 0.3) +
    annotate(geom = "text", x = p2, y = 0.64, label = expression(italic(r2))) +
    geom_segment(x = 1820, y = 0.61, xend = 1860, yend = 0.68, linewidth = 0.3) +
    geom_segment(x = 1820, y = 0.61, xend = 1860, yend = 0.61, linewidth = 0.3) +
    geom_segment(x = 1860, y = 0.61, xend = 1860, yend = 0.68, linewidth = 0.3) +
    # Highlight K:
    annotate(geom = "text", x = 750, y = 1, label = expression(italic(K)~"="~"inv.logit ("*italic(K[L])*")")) +
    
    lims(y = c(0, 1),
         x = c(0, 2000)) +
    labs(x = "DBH (mm)",
         y = expression(paste(Annual~survival~probability~(italic(theta))))) +
    theme_classic() +
    mytheme)

# Zoom in portion of fig_1:
(fig_2 <- fig_1 +
    coord_cartesian(xlim = c(1000, 1400), ylim = c(0.85, 1)) +
    theme_void() +
    guides(colour = "none") #+
  # theme(aspect.ratio = 1)
)

(fig_merge <- fig_1 + 
    # Add contour for inset:
    geom_segment(x = 650, xend = 1350, y = 0.07, yend = 0.07, linewidth = 0.2) +
    geom_segment(x = 650, xend = 1350, y = 0.73, yend = 0.73, linewidth = 0.2) +
    geom_segment(x = 650, xend = 650, y = 0.07, yend = 0.73, linewidth = 0.2) +
    geom_segment(x = 1350, xend = 1350, y = 0.07, yend = 0.73, linewidth = 0.2) +
    # Add contour for the magnified portion of fig_1:
    geom_segment(x = 1000, xend = 1170, y = 0.84, yend = 0.84, linewidth = 0.2) +
    geom_segment(x = 1000, xend = 1170, y = 1, yend = 1, linewidth = 0.2) +
    geom_segment(x = 1000, xend = 1000, y = 0.84, yend = 1, linewidth = 0.2) +
    geom_segment(x = 1170, xend = 1170, y = 0.84, yend = 1, linewidth = 0.2) +
    # Link the two rectangles:
    geom_segment(x = 1000, xend = 650, y = 0.84, yend = 0.73, lty = 2, linewidth = 0.2) +
    geom_segment(x = 1170, xend = 1350, y = 0.84, yend = 0.73, lty = 2, linewidth = 0.2) +
    # Add annotations in the inset space:
    annotate(geom = "text", x = 1000, y = 0.6, label = expression(italic(K*"_"*mu))) +
    annotate(geom = "text", x = 1000, y = 0.5, label = expression(italic(K[L~k*","~t])~"="~italic(K*"_"*mu)~+~italic(K_P[k])~+~italic(K_T[t])),
             colour = "purple3") +
    annotate(geom = "text", x = 1250, y = 0.67, label = expression(site[italic(k)])) +
    annotate(geom = "text", x = 1250, y = 0.11, label = expression(year[italic(t)])) +
    annotate(geom = "text", x = 1210, y = 0.43, label = expression(site[italic(k)]*","~year[italic(t)])) +
    # Add arrows for K_P and K_T:
    # K_T:
    geom_segment(x = 680, xend = 680, y = 0.16, yend = 0.54, linewidth = 0.4,
                 arrow = arrow(ends = "both", length = unit(0.1, "inches")),
                 colour = "darkred") +
    annotate(geom = "text", x = 750, y = 0.32, label = expression(italic(K_T[t])), colour = "darkred") +
    # K_P:
    geom_segment(x = 680, xend = 680, y = 0.58, yend = 0.69, linewidth = 0.4,
                 arrow = arrow(ends = "both", length = unit(0.1, "inches")),
                 colour = "blue") +
    annotate(geom = "text", x = 750, y = 0.64, label = expression(italic(K_P[k])), colour = "blue") +
    # Add the inset:
    inset_element(fig_2, 0.35, 0.1, 0.65, 0.8) 
)

# The above locations for the text and segments (and inset) are thought for a ration width to height of 8.8 to 4.8.
ggsave(plot = fig_merge,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S1_Surv_Model/Fig. S1.pdf",
       width = 8.8, height = 4.8)

ggsave(plot = fig_merge,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S1_Surv_Model/Fig. S1.png",
       width = 8.8, height = 4.8, dpi = 400)

# The png figure is then opened in PPT to add the two equations of the logistic functions on top, and save the final figure from the PPT file.


# Spatial analyses of K_tmp (Moran's eigenvector maps, MEM) ----
# *********************************************************
# *********************************************************

# Function to ease the visualization of connectivity matrices:
nb2ggplot <- function(nb, coord) {
  # 'coord' must be a matrix/dataframe with two columns (called "long" and "lat")
  # 'nb' is an object of class nb
  
  # take out the connections from the nb object
  # and assign them the lat and long in a dataframe
  n <- length(attributes(nb$neighbours)$region.id)
  DA <- data.frame(
    from = rep(1:n, sapply(nb$neighbours, length)),
    to = unlist(nb$neighbours),
    weight = unlist(nb$weights)
  )
  DA <- cbind(DA, coord[DA$from, 1:2], coord[DA$to, 1:2])
  colnames(DA)[4:7] = c("long2","lat2","long_to","lat_to")
  return(DA)
}

## Load data and model output ----
# ___________________________
# Adapt object name, depending on whether we want to analize the data at the grain of 50 or 100 km:
load("OUTPUT/df_Ktmp_plot_year_varying_K_111_50km_threshmean_lesscolumns.RData")
load("OUTPUT/df_Ktmp_plot_year_varying_K_121_50km_threshmean_lesscolumns.RData")
load("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData")

## Use MEM to test for and visualize spatial structures in cells' average K_tmp ----
# *****************************************************************************

# Spatial coordinates:
# ___________________
# Data for spatial analyses (MEMs) -> keep K_tmp on logit scale.
# For the spatial analysis, we summarize the survival plateau, K_tmp (K_L in manuscript), of
# each cell (50 x 50 km, or 100 x 100 km) across the multiple years available for each cell.

df_Ktmp2 <- df_Ktmp %>%
  rename(sp = species) %>%
  group_by(sp, plot2, long2, lat2) %>%
  summarise(K_tmp = mean(K_tmp),
            K_P = mean(K_mu + K_P)) %>%
  ungroup # %>% 
# # inverse logit transformation
# mutate(across(contains("K_"), ~boot::inv.logit(.)))

# Filter out an odd value (error):
df_Ktmp2 <- df_Ktmp2 %>% filter(K_tmp > 0) # for species 111
df_Ktmp2 <- df_Ktmp2 %>% filter(lat2 < 4600000) # for species 131

# Object of cells' coordinates:
xy <- df_Ktmp2 %>% 
  select(plot2, long2, lat2) %>% 
  column_to_rownames("plot2")

# Visualize plot network:
xy %>% 
  ggplot(aes(long2, lat2)) +
  geom_point(shape = 21, size = 1, alpha = 0.6, fill = "blue") + 
  coord_fixed() +
  theme_classic()

# Visualise logit-survival (K_tmp) on the map:
df_Ktmp2 %>% 
  ggplot(aes(long2, lat2)) +
  geom_point(shape = 21, size = 3, alpha = 0.6, aes(fill = K_tmp)) + 
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_classic()

### Generate candidate spatial weighting matrices (SWMs) ----
# ____________________________________________________
# Irregular sampling design --> favours graph-based connectivity matrices over
# distance-based MEMs (Bauman et al. 2018, Ecology). 
# Not too many candidates, to balance gain of accuracy from optimization procedure
# with the loss of stat. power (due to p-value adjustment for multiple tests).

swm <- adespatial::listw.candidates(xy,
                                    nb = c("gab", "mst"),
                                    weights = c("flin", "fdown"), y_fdown = 5)

# Visualise some graph-based connectivity schemes:
# _______________________________________________
# nbtri <- tri2nb(xy) # Delaunay triangulation
nbgab <- graph2nb(gabrielneigh(as.matrix(xy), nnmult = 4), sym = TRUE) # Gabriel graph
# nbrel <- graph2nb(relativeneigh(xy), sym = TRUE) # Relative neighbourhood graph
nbmst <- mst.nb(dist(xy)) # minimum spanning tree

# Other options visualisable through:
# other <- chooseCN(xy)
# listw.explore()

# The spdep package is designed to create object plotted with base R.
# This is not very handy, so that I wrote a function to be able to plot the links in ggplot2.
# We first need to transform the nb object into listw objetcts:

# nbtri_listw <- nb2listw(nbtri) # 
nbgab_listw <- nb2listw(nbgab)
# nbrel_listw <- nb2listw(nbrel)
nbmst_listw <- nb2listw(nbmst)

# Visualise connectivity matrix:

DA_gab <- nb2ggplot(nbgab_listw, xy)
# Figure:
xy %>% 
  as.data.frame() %>% 
  ggplot(aes(long2, lat2)) +
  geom_point(size = 1) +
  coord_fixed() +
  geom_segment(data = DA_gab, 
               aes(xend = long_to, yend = lat_to), 
               size = 0.3, alpha = 0.5, colour = "darkred") +
  labs(title = "Gabriel graph") +
  theme_bw()

DA_mst <- nb2ggplot(nbmst_listw, xy)

# Figure:
xy %>% 
  as.data.frame() %>% 
  ggplot(aes(long2, lat2)) +
  geom_point() +
  coord_fixed() +
  geom_segment(data = DA_mst, 
               aes(xend = long_to, yend = lat_to), size = 0.3, alpha = 0.5, colour = "blue") +
  labs(title = "MST") +
  theme_bw()

# If deemed necessary, some links can be edited (i.e. removed/added) in
# the above connectivity matrices, before building MEM variables (e.g. they cross
# physical barriers).

### Test for spatial structure and select SWM and MEM variables ----
# ___________________________________________________________

# Is there a clear spatial structure in the response variable(s)?
# If so, what does it look like?

# Test for spatial structure in the response data, and select an optimal subset
# of MEM variables in the best SWM using a selection criterion (here, the 
# forward selection with double stopping criterion, Blanchet et al. 2008,
# Bauman et al. 2018, Ecography):

spa_select <- adespatial::listw.select(df_Ktmp2$K_tmp, 
                                       swm,
                                       MEM.autocor = "positive",
                                       method = "FWD",
                                       p.adjust = T)

# Save:
save(spa_select, file = "OUTPUT/spa_select_sp111_50km.RData")
save(spa_select, file = "OUTPUT/spa_select_sp121_50km.RData")
save(spa_select, file = "OUTPUT/spa_select_sp131_50km.RData")

# Selected SWM:
spa_select$best.id
# Summarises the candidate SWMs and their characterisitcs:
spa_select$candidates
# Summarises MEM var. selected: 
spa_select$best$summary

# Selected MEM variables within 'select$best.id':
MEM.sel <- spa_select$best$MEM.select # to be used as explanatory variables 

### Visualization of the selected MEM variables ----
# ____________________________________________

# Prepare data for figure:
data_fig <- MEM.sel %>% 
  as_tibble %>% 
  mutate(plot = row.names(xy)) %>% 
  select(plot, everything()) %>% 
  # Transform into long format:
  pivot_longer(cols = 2:(ncol(MEM.sel)+1), # all MEM columns
               names_to = "MEM", values_to = "value") %>% 
  left_join(xy %>% 
              rownames_to_column("plot")) %>% 
  mutate(colour = ifelse(value > 0, "black", "white"))

mem_names <- colnames(MEM.sel)

# To generate single figure of multiple selected MEMs:
gg_list <- list()
for(i in 1:ncol(MEM.sel)) {
  gg_list[[i]] <- data_fig %>% 
    filter(MEM == mem_names[i]) %>% 
    ggplot(aes(long2, lat2)) +
    geom_point(shape = 21, alpha = 0.3,
               aes(size = abs(value), fill = colour)) +
    # scale_fill_viridis_d(option = "viridis") +
    scale_fill_manual(values = c("mediumblue", "red")) +
    scale_size_continuous(range = c(0.2, 3)) +
    coord_equal() +
    labs(subtitle = mem_names[i],) +
    guides(fill = "none", size = "none") +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.subtitle = element_text(size = 9))
}

# Explore:
do.call(gridExtra::grid.arrange, c(gg_list[1:10], # MEMs we want to visualise 
                                   ncol = 4)) # number of columns

# Create figure to be saved:
# 131, 50 km: 1:16
# 131, 100 km: 1:10
fig <- do.call(gridExtra::grid.arrange, c(gg_list[1:10], # MEMs we want to visualise 
                                          ncol = 4)) # number of columns

# Save:

# Note that these are not Fig. S6, eventhough they are saved in the same folder.
# These supplementary figures are not in the published manuscript, even though visualizing
# them is handy to better understand the way Moran's eigenvector maps are selected and 
# combined through a fitted linear combination to highlight complex spatial patterns 
# (resulting from combinations of these individual spatial eigenvectors).

ggsave(fig,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S6_MEM/Fig_selected_MEMs_131_50km.pdf", 
       width = 9, height = 7.5, units = "in")
ggsave(fig,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S6_MEM/Fig_selected_MEMs_131_50km.png", 
       width = 9, height = 7.5, units = "in", dpi = 400)

## Fig. S6 - Visualization of the resulting spatial structure in K_tmp ----
# (fitted linear combination of the selected MEM variables)
# _________________________________________________________

# Only generated for P. taeda (131), since no significant spatial structure is detected is the
# other species.

scores <- fitted(lm(df_Ktmp2$K_tmp ~ ., data = MEM.sel))

data_fig2 <- scores %>% 
  as_tibble() %>% 
  mutate(plot = row.names(xy)) %>% 
  select(plot, everything()) %>% 
  left_join(xy %>% 
              rownames_to_column("plot")) 

# Make nice figure of the fitted model (spatial pattern in K_tmp) on map of the USA:
# ____________________________________

# Base map for pine survival analysis:
world <- ne_countries(scale = "medium", returnclass = "sf") # load world map
NoAm <- world[world$continent=='North America',] # subset out only North America

d <- st_transform(st_as_sf(data_fig2, coords = c("long2", "lat2"), 
                           crs = 32616), # currently, Projected CRS: WGS 84 / UTM zone 16N
                  crs = 4326) # transform into WGS84, to match the coord. system used for the map (NoAm)

# No spatial pattern detected in P. elliottii or P. palustris.

# Sp 131:
# _______
R2_adj <- round(spa_select$candidates$R2Adj.select[spa_select$best.id], 3) # adjusted R2
p_adj <- round(spa_select$best$global.test$pvalue, 4) # p-value adjusted for multiple tests (candidate SWMs)

fig_fitted_131_50km <- ggplot(data = NoAm) +
  geom_sf(fill = 'white') +
  geom_sf(data = d, shape = 21, size = 4, alpha = 0.8, # size = 4 for 50 km, size = 8 for 100 km
          aes(fill = value, colour = value)) +
  # geom_sf(data = g50dat, fill = NA, colour = "grey80") +  # could join on the cell id to the mean values for fill
  coord_sf(xlim = c(-100, -73), ylim = c(24, 40), expand = FALSE) +
  scale_fill_viridis(option = "inferno", name = expression(Fitted~italic(K[L]))) +
  scale_colour_viridis(option = "inferno", name = expression(Fitted~italic(K[L]))) +
  
  labs(subtitle = expression(Spatial~pattern~detected~"in"~italic("P. taeda")*"'s survival"~"("*italic(K[L])*")"~across~"2003-2022")) +
  
  # annotate(geom = 'text', x = -97, y = 38.5, size = 14/.pt, 
  #          label = expression(italic("P. taeda"))) +
  annotate(geom = 'text', x = -94, y = 39, size = 12/.pt,
           label = substitute(
             italic(R)*"²"[adj]~"="~nn*","~italic(p)*"-value"[adj]~"="~oo,
             list(nn = R2_adj, oo = format(p_adj, scientific = FALSE)))) +
  
  mytheme +
  theme(panel.background = element_rect(fill = 'lightblue'),
        legend.background = element_rect(fill = "lightblue"),
        legend.position = c(.9, .35),
        axis.title = element_blank()) 

fig_fitted_131_50km

saveRDS(fig_fitted_131_50km, 
        file = "OUTPUT/figures/figures_supp_manuscript/Fig_S6_MEM/fig_fitted_131_50km.rds")

ggsave(fig_fitted_131_50km,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S6_MEM/Fig_fitted_MEMs_131_50km.pdf", 
       width = 9, height = 7.5, units = "in")
ggsave(fig_fitted_131_50km,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S6_MEM/Fig_fitted_MEMs_131_50km.png", 
       width = 9, height = 7.5, units = "in", dpi = 400)

# Rerun the spatial analysis above from the 100 km grain, and create corresponding figure.

# Fig. S6, of MEMs fitted to K_tmp, at 50 and 100 km grain:
# _________________________________________________________
# First, create fig_fitted_131_50km and _100km (or load them).

library(patchwork)

fig_s <- fig_fitted_131_50km + theme(plot.subtitle = element_blank()) + 
  fig_fitted_131_100km + theme(plot.subtitle = element_blank()) +
  plot_annotation(tag_levels = "a",
                  title = expression(Spatial~pattern~detected~"in"~italic("P. taeda")*"'s survival"~"("*italic(K[L])*")"~across~"2003-2022"),
                  theme = theme(plot.title = element_text(hjust = 0.5)))

fig_s

# Save Fig. S6:
ggsave(fig_s,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S6_MEM/Fig_S6_fitted_MEMs_131_50km_and_100km.pdf", 
       width = 14, height = 5.5, units = "in")
ggsave(fig_s,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S6_MEM/Fig_S6_fitted_MEMs_131_50_and_100km.png", 
       width = 14, height = 5.5, units = "in", dpi = 400)


# Temporal variation and trends in species' survival (K_tmp) - Fig. 3 and Fig. S8 ----
# *******************************************************************************
# *******************************************************************************

# Run the code below on the grain of 50 km and 100 km to generate the elements of Fig. 3 and Fig. S8,
# respectively.

## Species 111 ----
# ************
# Load model output:
load("OUTPUT/df_Ktmp_plot_year_varying_K_111_50km_threshmean_lesscolumns.RData")

fig111time_50km <- df_Ktmp %>% 
  filter(year >= 2004) %>% 
  # mutate(K_tmp = boot::inv.logit(K_tmp)) %>% 
  mutate(K_T = boot::inv.logit(K_mu + K_T)) %>% 
  # Figure:
  ggplot(aes(year, K_T)) +
  # stat_halfeye(.width = c(0.75, 0.9)) +
  stat_pointinterval(.width = c(0.75, 0.9)) +
  coord_cartesian(ylim = c(0.85, 1)) +
  labs(x = "Years", y = "Survival probability",
       subtitle = expression(italic(P*"."~elliottii))) +
  mytheme

fig111time_50km + labs(x = "") # to remove axis title without changing the width to height ratio

# Save:
saveRDS(fig111time_50km, file = "OUTPUT/figures/figures_manuscript/Fig. 3/fig111time_50km.rds")

ggsave(plot = fig111time_50km,
       filename = "OUTPUT/figures/figures_manuscript/Fig. 3/Fig3_elliottii_per_year.pdf",
       width = 6.7, height = 4, units = "in")

# # Supplementary figure with 2003, and full y-axis:
# fig111time_50km_supp <- df_Ktmp %>% 
#   # filter(year >= 2004) %>% 
#   # mutate(K_tmp = boot::inv.logit(K_tmp)) %>% 
#   mutate(K_T = boot::inv.logit(K_mu + K_T)) %>% 
#   # Figure:
#   ggplot(aes(year, K_T)) +
#   # stat_halfeye(.width = c(0.75, 0.9)) +
#   stat_pointinterval(.width = c(0.75, 0.9)) +
#   coord_cartesian(ylim = c(0.5, 1)) +
#   labs(x = "Years", y = "Survival probability",
#        subtitle = expression(italic(P*"."~elliottii))) +
#   mytheme
# fig111time_50km_supp
# 
# ggsave(plot = fig111time_50km_supp,
#        filename = "OUTPUT/figures/figures_supp_manuscript/Fig3_suppmat_elliottii_per_year.pdf", width = 6.7, height = 4, units = "in")
# 
# ggsave(plot = fig111time_50km_supp + coord_cartesian(ylim = c(0, 1)),
#        filename = "OUTPUT/figures/figures_supp_manuscript/Fig3_suppmat_elliottii_per_year_0_to_1.pdf", width = 6.7, height = 4, units = "in")

## Species 121 ----
# **************
# Load model output:
load("OUTPUT/df_Ktmp_plot_year_varying_K_121_50km_threshmean_lesscolumns.RData")

fig121time_50km <- df_Ktmp %>% 
  filter(year >= 2004) %>% 
  # mutate(K_tmp = boot::inv.logit(K_tmp)) %>% 
  mutate(K_T = boot::inv.logit(K_mu + K_T)) %>% 
  # Figure:
  ggplot(aes(year, K_T)) +
  # stat_halfeye(.width = c(0.75, 0.9)) +
  stat_pointinterval(.width = c(0.75, 0.9)) +
  coord_cartesian(ylim = c(0.85, 1)) +
  labs(x = "Years", y = "Survival probability",
       subtitle = expression(italic(P*"."~palustris))) +
  mytheme

fig121time_50km + labs(x = "") # to remove axis title without changing the width to height ratio

saveRDS(fig121time_50km, file = "OUTPUT/figures/figures_manuscript/Fig. 3/fig121time_50km.rds")

ggsave(plot = fig121time_50km,
       filename = "OUTPUT/figures/figures_manuscript/Fig. 3/Fig3_palustris_per_year.pdf", width = 6.7, height = 4, units = "in")

## Species 131 ----
# ____________
# Load model output:
load("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData")  # sp 131

fig131time_50km <- df_Ktmp %>% 
  filter(year >= 2004) %>% 
  # mutate(K_tmp = boot::inv.logit(K_tmp)) %>% 
  mutate(K_T = boot::inv.logit(K_mu + K_T)) %>% 
  # Figure:
  ggplot(aes(year, K_T)) +
  # stat_halfeye(.width = c(0.75, 0.9)) +
  stat_pointinterval(.width = c(0.75, 0.9)) +
  coord_cartesian(ylim = c(0.85, 1)) +
  labs(x = "Years", y = "Survival probability",
       subtitle = expression(italic(P*"."~taeda))) +
  mytheme
fig131time_50km

ggsave(plot = fig131time_50km,
       filename = "OUTPUT/figures/figures_manuscript/Fig. 3/Fig3_taeda_per_year.pdf", width = 6.7, height = 4, units = "in")

saveRDS(fig131time_50km, file = "OUTPUT/figures/fig131time_50km.rds")

# Fig. 3, merging figures of the three species is now built in Inkscape.

# ### Single figure for the three species: (now done in Inkscape)
# # **************************************
# # **************************************
# fig111time_50km <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 3/fig111time_50km.rds")
# fig121time_50km <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 3/fig121time_50km.rds")
# fig131time_50km <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 3/fig131time_50km.rds")
# 
# (fig111time_50km + theme(axis.title.x = element_blank())) / 
#   (fig121time_50km + theme(axis.title.x = element_blank())) / 
#   (fig131time_50km)
# 
# grid.arrange(fig111time_50km + theme(axis.title.x = element_blank()), 
#              fig121time_50km + theme(axis.title.x = element_blank()),
#              fig131time_50km,
#              cols = 1)


# Multinomial (softmax) regression models ----
# ***************************************
# ***************************************

## Create summarized K_tmp posteriors for multinomial models ----
# __________________________________________________________

# Run for the three species and two spatial grains: (here only ran for one species and one grain)
load("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData")

df_Ktmp_summarised <- df_Ktmp %>%
  group_by(stem, year, plot2, long2, lat2) %>%
  summarise(K_tmp = median(K_tmp)) %>%
  ungroup

saveRDS(df_Ktmp_summarised, file = "OUTPUT/df_Ktmp_median_stemyear_131_50km.rds")

## Prepare dataset for softmax analyses ----
# _____________________________________

# Chose the spatial grain by replacing all "50" by "100" or vice versa, for the 50 x 50 km
# and 100 x 100 km grid, respectively.

# Select species:
sp.id <- 111

# Load data:
meta.tmp <- read.csv(file = "DATA/data2_50km_2023.csv", na = c("", "NA"))

# Retrieve only data with mortality agents listed:
meta.sub <- meta.tmp[which(!is.na(meta.tmp$agent)), 
                     c("dbh", "stem", "Sp", "year_0", "year_1", "surv",
                       "agent")]

# Keep only the focal species:
meta.df <- subset(meta.sub, Sp == sp.id)

# Load survival model output for the corresponding species:
# (median K_tmp per tree per plot2 (cell) and per year):
df_Ktmp <- readRDS(file = sprintf("OUTPUT/df_Ktmp_median_stemyear_%i_50km.rds", sp.id))

mort.ag <- factor(meta.df$agent) # vector of mort. agent for all the trees of sp.id that died

# Create a key to merge the mortality agent info to the survival model output (df_Ktmp):
meta.df$stemyear <- paste0(meta.df$stem, meta.df$year_0) 
df_Ktmp$stemyear <- paste0(df_Ktmp$stem, df_Ktmp$year)

# Add mortality agent info to survival model output:
df_Ktmp.a <- merge(as.data.frame(df_Ktmp), meta.df, by = "stemyear")

# Create an explicit column for the mortality agent (named "m.cause" here):
df_Ktmp.a$m.cause <- factor(df_Ktmp.a$agent, 
                            labels = c("Insects", "Disease", 
                                       "Fire", "Animal", "Weather", "Competition", "Unknown"))

# df_Ktmp.a$n.cause <- factor(as.integer(df_Ktmp.a$agent/10))

# Make Competition the level of reference for the m.cause factor, so the softmax regression
# will use that agent as point of comparison for the others:
df_Ktmp.a$m.cause <- relevel(df_Ktmp.a$m.cause, ref = "Competition")

table(data.frame(Code = seq(10, 70, 10), 
                 Cause = c("Insects", "Disease", 
                           "Fire", "Animal", "Weather", "Competition", "Unknown")))
table(df_Ktmp.a$m.cause)

# Save:
save(df_Ktmp.a, file = sprintf("OUTPUT/DF_Ktmp_%i_50km.Rdata", sp.id))

# Code    Description
# 10  Insects
# 20  Disease
# 30  Fire
# 40  Animal
# 50  Weather
# 60  Vegetation (competition/suppression)
# 70  Unknown (or multiple causes)
# 80  Human cause death # not present here, as these were removed prior to running the surv. models
# table(df_Ktmp.a$m.cause)

## Run spatial softmax analysis ----
# ****************************
library(brms)

# if run on the cluster, re create object sp.id for corresponding species.
sp.id <- 111
load(file = sprintf("OUTPUT/DF_Ktmp_%i_50km.Rdata", sp.id)) # loads df_Ktmp.a

agent.sp.cell <-
  brm(data = df_Ktmp.a,
      family = categorical(link = logit),
      m.cause ~ 1 + K_tmp + (1 + K_tmp | plot2),
      iter = 2000, cores = 4, chains = 4,
      control = list(adapt_delta = 0.97, max_treedepth = 15))

# Save model output:
saveRDS(agent.sp.cell,
        file = sprintf("OUTPUT/softmax_out_50_cell_%i_new.rds", sp.id))

## Run temporal softmax analysis ----
# ******************************
library(brms)

# if run on the cluster, re create object sp.id for corresponding species.
sp.id <- 111
load(file = sprintf("OUTPUT/DF_Ktmp_%i_50km.Rdata", sp.id))

agent.sp.cell.yr <- 
  brm(data = df_Ktmp.a, 
      family = categorical(link = logit),
      m.cause ~ 1 + K_tmp + (1 + K_tmp | year_0),
      iter = 2000, cores = 4, chains = 4, silent = 2,
      control = list(adapt_delta = 0.97, max_treedepth = 15))

# Save model output:
saveRDS(agent.sp.cell.yr,
        file = sprintf("OUTPUT/softmax_out_50_cell_year_%i_new.rds", sp.id))


## Extract softmax models' posteriors and generate conditional effects ----
# ********************************************************************

# Replace 100 by 50 or 50 by 100 below, depending on the grain desired.

### Spatial softmax regression ----
# ***************************

# Select a focal species: (only run for one species here --> repeat code below for the other two)
sp.id <- 111

agent.sp.cell <- readRDS(file = sprintf("OUTPUT/softmax_out/softmax_out_50_cell_%i_new.rds", sp.id))

post <- gather_draws(agent.sp.cell, `b_mu.*.tmp`, regex = TRUE)
post_summary <- point_interval(post, .point = median, .interval = hdi, .width = 0.90)

# Posterior predictions:
# ______________________

new.data <- df_Ktmp.a[, c("plot2", "K_tmp", "m.cause")]

e.pred.cell <- add_epred_draws(newdata = new.data,
                               object = agent.sp.cell, 
                               re_formula = NULL, 
                               ndraws = 1000)

# create a matrix summarizing the output from e.pred.cell
cell_vals <- tapply(e.pred.cell$.epred,
                    list(e.pred.cell$.category, e.pred.cell$plot2), mean)

# head(cell_vals)

# Save (will be used to generate the figure):
save(cell_vals, file = sprintf("OUTPUT/softmax_out/Pred_matrix_%i_50km.Rdata", sp.id))

# Conditional effects using K_tmp:
# ________________________________

# Get conditional estimates for cells:
load(sprintf("OUTPUT/df_Ktmp_plot_year_varying_K_%i_50km_threshmean_lesscolumns.RData", sp.id))

# For sp.id <- 131:
# load(sprintf("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData", sp.id))

# For the figure displaying all plots, but with plot2 mean K-tmp for all plots within a given plot2:
new.data <- df_Ktmp %>%
  group_by(plot2) %>%
  summarise(K_tmp = mean(K_tmp)) %>%
  ungroup

fit.data <- agent.sp.cell$data # data used in the fit

conditions <- data.frame(plot2 = unique(fit.data$plot2))

rownames(conditions) <- unique(fit.data$plot2)
head(conditions)

plot_Ktmp <- df_Ktmp$K_tmp[match(conditions$plot2, new.data$plot2)]
names(plot_Ktmp) <- conditions$plot2

# If working with sp.id 131:
# rm(df_Ktmp) # Need to clear some memory, otherwise cond.eff hits a RAM wall

# Conditional effects per cell: (will be used for the figure)
cond.eff <- conditional_effects(agent.sp.cell, 
                                categorical = TRUE,
                                conditions = conditions,
                                re_formula = NULL,
                                int_conditions = list(K_tmp = plot_Ktmp))

# # Conditional effects across cells:
# cond.eff_grand <- conditional_effects(agent.sp.cell, effects = "K_tmp",
#                                       categorical = TRUE,
#                                       conditions = NULL, re_formula = NULL,
#                                       int_conditions = list(K_tmp = plot_Ktmp),
#                                       spaghetti = FALSE, prob = 0.5)

save(cond.eff, file = sprintf("OUTPUT/softmax_out/Cond_eff_%i_50km.Rdata", sp.id))

rm(agent.sp.cell, cond.eff) # free memory for what comes next


### Temporal softmax regression ----
# ***************************

sp.id <- 131

# fitted model
agent.sp.yr <- readRDS(file = sprintf("OUTPUT/softmax_out/softmax_out_50_cell_year_%i_new.rds", sp.id))

# Get conditional estimates for cells for non-year model:
load(sprintf("OUTPUT/df_Ktmp_plot_year_varying_K_%i_50km_threshmean_lesscolumns.RData", sp.id))

# For sp.id <- 131:
# load(sprintf("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData", sp.id))

# Prepare the grid of covariate combination for which we want predictions:
dat.grid <- datagrid(model = agent.sp.yr, 
                     K_tmp = seq(2, 5, length = 100), # 100 K_tmp values between 2 and 5 (see Methods)
                     year_0 = seq(2003, 2018))

cond.eff <- marginaleffects::predictions(
  agent.sp.yr,
  newdata = dat.grid,
  # by = c("K_tmp", "year_0"),
  conf_level = 0.9,
  re_formula = NULL
)

cond.eff.df <- as.data.frame(cond.eff)
sort.ix <- order(cond.eff.df$group, cond.eff.df$year_0)
cond.eff.df2 <- cond.eff.df[sort.ix, ]

# Save, for figure:
write.csv(cond.eff.df2[, -c(1, 8)], # cond.eff.df2[, -c(1, 6)],
          file = sprintf("OUTPUT/softmax_out/Cond_effects_year_%i_100km.csv", sp.id), row.names = FALSE,
          quote = FALSE)
cond.eff.df <- read.csv(file = sprintf("OUTPUT/Cond_effects_year_%i_100km.csv", sp.id))

# Code for sp.id 131: (to avoid reaching memory limit)
# __________________
df_Ktmp <- select(df_Ktmp, -c(species:lat2, dbh:surv, K_mu:r2))

new.data <- df_Ktmp %>%
  group_by(year) %>%
  summarise(K_tmp = mean(K_tmp)) %>%
  ungroup 

fit.data <- agent.sp.yr$data # data used in the fit

conditions <- data.frame(year_0 = unique(fit.data$year))

rownames(conditions) <- unique(fit.data$year)
head(conditions)

year_Ktmp <- df_Ktmp$K_tmp[match(conditions$year, new.data$year)]

names(year_Ktmp) <- conditions$year

# If working with sp.id 131:
rm(df_Ktmp) # Need to clear some memory, otherwise cond.eff hits a RAM wall

# Conditional effects per year:
cond.eff <- conditional_effects(agent.sp.yr, categorical = TRUE, 
                                conditions = conditions, re_formula = NULL, 
                                int_conditions = list(K_tmp = year_Ktmp))

save(cond.eff, file = sprintf("OUTPUT/softmax_out/Cond_eff_%i_year_50km.Rdata", sp.id))


# Fig. 2 - Maps of species' survival (K_L) and dominant mortality agents (spatial softmax) ----
# **********************************************************************
# **********************************************************************

# The figure is made of two columns of three panels each. The left-hand side shows the
# map of the mean K_L (logit survival plateau) of the species across the multiple censuses,
# whereas the right-hand side displays the cell-level dominant mortality agent.

## Map of species' survival (K_L) ----
# *******************************

# Theme for figures:
mytheme <- theme_classic() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 15),
        plot.caption = element_text(size = 11))

# Load data:
# **********
# data <- readRDS("DATA/data_s_pines_50_and_100km_from2003.rds") # old version of the dataset -> find up-to-date version of this for proper match of cell names with surv_grid_.RDS
data <- fread("DATA/data2_50km_2023.csv") # what sp.df was based on
data$plot2 <- str_remove(data$plot2, "c_")

# Add plot coordinates (only the cell coordinates are available in 'data'):
coord <- fread("OUTPUT/Pine Plot locations 2023.csv")
data <- left_join(data, coord)

# Base map for pine survival analysis:
# ************************************
world <- ne_countries(scale = "medium", returnclass = "sf") # load world map
NoAm <- world[world$continent=='North America',] # subset out only North America

g50 <- readRDS("DATA/surv_grid_50.RDS") # read in grid 50
g50t <- st_transform(g50, crs  = st_crs(NoAm)) # align coordinate systems

g100 <- readRDS("DATA/surv_grid_100.RDS") # read in grid 100
g100t <- st_transform(g100, crs  = st_crs(NoAm)) # align coordinate systems

# Common scale of color for all three species, to map average K_tmp:
lower_col <- 0.7

### Species 111 ----
# _____________
load("OUTPUT/df_Ktmp_plot_year_varying_K_111_50km_threshmean_lesscolumns.RData")

# Correct some missing long lat:
df_Ktmp <- df_Ktmp %>% 
  select(-c(long, lat)) %>% 
  left_join(coord)

# For the figure displaying all plots, but with plot2 mean K-tmp for all plots within a given plot2:
df_Ktmp2 <- df_Ktmp %>%
  rename(sp = species) %>%
  group_by(sp, plot2, long, lat) %>%
  summarise(K_tmp = mean(K_tmp),
            K_P = mean(K_mu + K_P)) %>%
  ungroup %>%
  # inverse logit transformation
  mutate(across(contains("K_"), ~boot::inv.logit(.)))

# Only keep grids where data are present:
g50dat <- subset(g50t, g50t$cell_50 %in% subset(data, Sp == "111")$plot2) # simple subset of only cells with data

# transform coordinates into WGS84 to match US map crs:
d <- st_transform(st_as_sf(df_Ktmp2, coords = c("long", "lat"), crs = 4326),
                  crs = 4326)

# Create the figure:
figspa111_50 <-
  ggplot(data = NoAm) +
  geom_sf(fill = 'white') +
  geom_sf(data = d, shape = 21, size = 1, aes(fill = K_tmp, colour = K_tmp)) +
  geom_sf(data = g50dat, fill = NA, colour = "grey80") +  # could join on the cell id to the mean values for fill
  coord_sf(xlim = c(-100, -73), ylim = c(24, 43), expand = FALSE) +
  scale_fill_viridis(name = expression(~italic("K")), limits = c(lower_col, 1)) +
  scale_colour_viridis(name = expression(~italic("K")), limits = c(lower_col, 1)) +
  labs(x= "", y = "") +
  annotate(geom = 'text', x= -97, y = 38, size = 14/.pt, 
           label = expression(~italic("P. elliottii"))) +
  mytheme +
  theme(panel.background = element_rect(fill = 'lightblue'),
        legend.background = element_rect(fill = "lightblue"),
        legend.position = c(.9, .35)) 

# Show and save the figure zoomed in on the region of interest (for final figure):
figspa111_50 +
  coord_sf(xlim = c(-100, -73), ylim = c(24, 40), expand = FALSE)

# Save:
saveRDS(figspa111_50, file = "OUTPUT/figures/figures_manuscript/Fig. 2/figspa111_50km.rds")

# # Save the legend separately:
# legend <- cowplot::get_legend(figspa111_50)	
# ggpubr::as_ggplot(legend)
# 
# saveRDS(legend, file = "OUTPUT/figures/figures_manuscript/Fig. 2/figspa_legend_50km.rds")

### Species 121 ----
# _____________
load("OUTPUT/df_Ktmp_plot_year_varying_K_121_50km_threshmean_lesscolumns.RData")

# Correct some missing long lat:
df_Ktmp <- df_Ktmp %>% 
  select(-c(long, lat)) %>% 
  left_join(coord)

# For the figure displaying all plots, but with plot2 mean K-tmp for all plots within a given plot2:
df_Ktmp2 <- df_Ktmp %>%
  rename(sp = species) %>%
  group_by(sp, plot2, long, lat) %>%
  summarise(K_tmp = mean(K_tmp),
            K_P = mean(K_mu + K_P)) %>%
  ungroup %>%
  # inverse logit transformation
  mutate(across(contains("K_"), ~boot::inv.logit(.)))

# Only keep grids where data are present:
g50dat <- subset(g50t, g50t$cell_50 %in% subset(data, Sp == "121")$plot2) # simple subset of only cells with data

# transform coordinates into WGS84 to match US map crs:
d <- st_transform(st_as_sf(df_Ktmp2, coords = c("long", "lat"), crs = 4326),
                  crs = 4326)

figspa121_50 <-
  ggplot(data = NoAm) +
  geom_sf(fill = 'white') +
  geom_sf(data = d, shape = 21, size = 1, aes(fill = K_tmp, colour = K_tmp)) +
  geom_sf(data = g50dat, fill = NA, colour = "grey80") +  # could join on the cell id to the mean values for fill
  coord_sf(xlim = c(-100, -73), ylim = c(24, 43), expand = FALSE) +
  scale_fill_viridis(name = expression(~italic("K"["L"])), limits = c(lower_col, 1)) +
  scale_colour_viridis(name = expression(~italic("K"["L"])), limits = c(lower_col, 1)) +
  labs(x = "", y = "") +
  annotate(geom = 'text', x= -97, y = 38, size = 14/.pt, 
           label = expression(~italic("P. palustris")))+
  mytheme +
  theme(panel.background = element_rect(fill = 'lightblue'),
        legend.position = 'none')

figspa121_50 +
  coord_sf(xlim = c(-100, -73), ylim = c(24, 40), expand = FALSE) 

saveRDS(figspa121_50, file = "OUTPUT/figures/figures_manuscript/Fig. 2/figspa121_50km.rds")

### Species 131 ----
# _____________
load("OUTPUT/df_Ktmp_plot_year_varying_K_131_50km_thresh300_lesscolumns.RData")

# Correct some missing long lat:
df_Ktmp <- df_Ktmp %>% 
  select(-c(long, lat))

df_Ktmp <- df_Ktmp %>% 
  select(-c(stem, dbh, p1, p2, r1, r2, K_T, surv, year, species)) 

df_Ktmp <- df_Ktmp %>% 
  left_join(coord)

# For the figure displaying all plots, but with plot2 mean K-tmp for all plots within a given plot2:
df_Ktmp2 <- df_Ktmp %>%
  mutate(sp = "131") %>%
  group_by(sp, plot2, long, lat) %>%
  summarise(K_tmp = mean(K_tmp),
            K_P = mean(K_mu + K_P)) %>%
  ungroup %>%
  # inverse logit transformation
  mutate(across(contains("K_"), ~boot::inv.logit(.)))

# Dan spotted one mistake, one plot too high up in the country to not be a wrong coordinate. We remove it here:
df_Ktmp2 <- df_Ktmp2[-which.max(df_Ktmp2$lat), ] # removes plot2 "c_18549" (plot2 from previous cell IDs, when they were first generated from all species, not only pine species); corresponding lat: 4695385

# Only keep grids where data are present:
g50dat <- subset(g50t, g50t$cell_50 %in% subset(data, Sp == "131")$plot2) # simple subset of only cells with data
# g50dat <- subset(g50t, g50t$cell_50 %in% subset(data, Sp == "131")$plot2 & g50t$lat_50 != max(g50t$lat_50)) # simple subset of only cells with data

d <- st_transform(st_as_sf(df_Ktmp2, coords = c("long", "lat"), crs = 4326),
                  crs = 4326)

figspa131_50 <-
  ggplot(data = NoAm) +
  geom_sf(fill = 'white') +
  geom_sf(data = d, shape = 21, size = 1, aes(fill = K_tmp, colour = K_tmp)) +
  geom_sf(data = g50dat, fill = NA, colour = "grey80") +  # could join on the cell id to the mean values for fill
  coord_sf(xlim = c(-100, -73), ylim = c(24, 43), expand = FALSE) +
  scale_fill_viridis(name = expression(~italic("K"["L"])), limits = c(lower_col, 1)) +
  scale_colour_viridis(name = expression(~italic("K"["L"])), limits = c(lower_col, 1)) +
  annotate(geom = 'text', x= -97, y = 38, size = 14/.pt, 
           label = expression(~italic("P. taeda")))+
  labs(x= "", y = "") +
  mytheme +
  theme(panel.background = element_rect(fill = 'lightblue'),
        legend.position = 'none')

figspa131_50 +
  coord_sf(xlim = c(-100, -73), ylim = c(24, 40), expand = FALSE) 

saveRDS(figspa131_50, file = "OUTPUT/figures/figures_manuscript/Fig. 2/figspa131_50km.rds")

# Final figure is generated in further down, after the spatial softmax figures are generated.


## Map of dominant mortality agents ----
# *********************************

# Theme for figures:
mytheme <- theme_classic() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 15),
        plot.caption = element_text(size = 11))

mypalette <- c("#5cd86f", "#ff9e5d", "#ffbcfb", "#947800", "#c00062", "#2b0b70")

### Species 111 ----
# _____________
sp.id <- 111

load(file = sprintf("OUTPUT/softmax_out/Pred_matrix_%i_50km.Rdata", sp.id))
c_vals <- as.data.frame(t(cell_vals))

c_vals$cell <- str_remove(row.names(c_vals), "c_")

c_vals_long <-c_vals %>% 
  pivot_longer(!cell, names_to = 'cause', values_to = 'prob')

c_vals_main <- c_vals_long %>% 
  group_by(cell) %>% 
  filter(prob == max(prob))

# Base map for pine survival analysis:
# ************************************
world <- ne_countries(scale = "medium", returnclass = "sf") # load world map
NoAm <- world[world$continent=='North America',] # subset out only North America

g50 <- readRDS("DATA/surv_grid_50.RDS") # read in grid 50
g50_probs <- merge(g50, c_vals_main, by.x = 'cell_50', by.y = 'cell')

# align coordinate systems
g50t <- st_transform(g50_probs, crs  = st_crs(NoAm)) 

# Create the figure:
fig_softmax_50km_111 <-
  ggplot(data = NoAm) +
  geom_sf(fill = 'white', color='grey10') +
  geom_sf(data = g50t,
          aes(fill = cause), colour = 'grey80') +  # replace "ESTIMATE" by the name of the column based on which cells should be coloured (cause-related)
  coord_sf(xlim = c(-100, -73), ylim = c(24, 43), expand = FALSE) +
  
  # scale_fill_viridis_d(name = "Agent of\nmortality") +
  # scale_colour_viridis_d(name = "Agent of\nmortality") + 
  
  scale_colour_manual(name = "Agent of\nmortality",
                      values = mypalette) +
  scale_fill_manual(name = "Agent of\nmortality",
                    values = mypalette) +  
  
  labs(x=NULL, y = NULL)+
  annotate(geom = 'text', x= -97, y = 38, size = 14/.pt, 
           label = expression(~italic("P. elliottii")))+
  mytheme +
  theme(panel.background = element_rect(fill = 'lightblue'),
        legend.background = element_rect(fill = "lightblue"),
        legend.position = c(.87, .24))

fig_softmax_50km_111 + coord_sf(xlim = c(-100, -73), ylim = c(24, 40), expand = FALSE)

# Save the figure object and the figure:
ggsave(fig_softmax_50km_111, 
       file = "OUTPUT/figures/figures_manuscript/Fig. 2/fig111_softmax_50km.pdf",
       width = 8, height = 7)

saveRDS(fig_softmax_50km_111, 
        file = "OUTPUT/figures/figures_manuscript/Fig. 2/fig_softmax_50km_111.rds")

# # Save the legend separately:
# legend <- cowplot::get_legend(fig_softmax_50km_111)	
# ggpubr::as_ggplot(legend)
# 
# saveRDS(legend, file = "OUTPUT/figures/figures_manuscript/fig_softmax_50km_111_legend.rds")

### Species 121 ----
# _____________
sp.id <- 121

load( file = sprintf("OUTPUT/softmax_out/Pred_matrix_%i_50km.Rdata", sp.id))
c_vals <- as.data.frame(t(cell_vals))

c_vals$cell <- str_remove(row.names(c_vals), "c_")

c_vals_long <-c_vals %>% 
  pivot_longer(!cell, names_to = 'cause', values_to = 'prob')

c_vals_main <- c_vals_long %>% 
  group_by(cell) %>% 
  filter(prob == max(prob))

g50 <- readRDS("DATA/surv_grid_50.RDS") # read in grid 50
g50_probs <- merge(g50, c_vals_main, by.x = 'cell_50', by.y = 'cell')

# align coordinate systems
g50t <- st_transform(g50_probs, crs  = st_crs(NoAm)) 

# Create the figure:
fig_softmax_50km_121 <-
  ggplot(data = NoAm) +
  geom_sf(fill = 'white') +
  geom_sf(data = g50t,
          aes(fill = cause), colour = 'grey80') +  # replace "ESTIMATE" by the name of the column based on which cells should be coloured (cause-related)
  coord_sf(xlim = c(-100, -73), ylim = c(24, 38), expand = FALSE) +
  
  # scale_fill_viridis_d(name = "Agent of\nmortality") +
  # scale_colour_viridis_d(name = "Agent of\nmortality") + 
  
  scale_colour_manual(name = "Agent of\nmortality",
                      values = mypalette) +
  scale_fill_manual(name = "Agent of\nmortality",
                    values = mypalette) +  
  labs(x=NULL, y = NULL)+
  annotate(geom = 'text', x= -96, y = 38, size = 14/.pt, 
           label = expression(~italic("P. palustris")))+ 
  mytheme +
  theme(panel.background = element_rect(fill = 'lightblue')) 

fig_softmax_50km_121 + coord_sf(xlim = c(-100, -73), ylim = c(24, 40), expand = FALSE)

# Save the figure object and the figure:
ggsave(fig_softmax_50km_121, 
       file = "OUTPUT/figures/figures_manuscript/Fig. 2/fig121_softmax_50km.pdf",
       width = 8, height = 7)

saveRDS(fig_softmax_50km_121, file = "OUTPUT/figures/figures_manuscript/Fig. 2/fig_softmax_50km_121.rds")

### Species 131 ----
# _____________
sp.id <- 131

load( file = sprintf("OUTPUT/softmax_out/Pred_matrix_%i_50km.Rdata", sp.id))
c_vals <- as.data.frame(t(cell_vals))

c_vals$cell <- str_remove(row.names(c_vals), "c_")

c_vals_long <-c_vals %>% 
  pivot_longer(!cell, names_to = 'cause', values_to = 'prob')

c_vals_main <- c_vals_long %>% 
  group_by(cell) %>% 
  filter(prob == max(prob))

g50 <- readRDS("DATA/surv_grid_50.RDS") # read in grid 50
g50_probs <- merge(g50, c_vals_main, by.x = 'cell_50', by.y = 'cell')

# align coordinate systems
g50t <- st_transform(g50_probs, crs  = st_crs(NoAm)) 

# Create the figure:
fig_softmax_50km_131 <-
  ggplot(data = NoAm) +
  geom_sf(fill = 'white') +
  geom_sf(data = g50t,
          aes(fill = cause), colour = 'grey80') +  # replace "ESTIMATE" by the name of the column based on which cells should be coloured (cause-related)
  coord_sf(xlim = c(-100, -73), ylim = c(24, 43), expand = FALSE) +
  
  # scale_fill_viridis_d(name = "Agent of\nmortality") +
  # scale_colour_viridis_d(name = "Agent of\nmortality") + 
  
  scale_colour_manual(name = "Agent of\nmortality",
                      values = mypalette) +
  scale_fill_manual(name = "Agent of\nmortality",
                    values = mypalette) +  
  labs(x=NULL, y = NULL)+
  annotate(geom = 'text', x= -97, y = 38, size = 14/.pt,
           label = expression(~italic("P. taeda")))+ 
  mytheme+
  theme(panel.background = element_rect(fill = 'lightblue')) 

fig_softmax_50km_131 + coord_sf(xlim = c(-100, -73), ylim = c(24, 40), expand = FALSE)

# Save the figure object and the figure:
ggsave(fig_softmax_50km_131, 
       file = "OUTPUT/figures/figures_manuscript/Fig. 2/fig131_softmax_50km.pdf",
       width = 8, height = 7)

saveRDS(fig_softmax_50km_131, file = "OUTPUT/figures/figures_manuscript/Fig. 2/fig_softmax_50km_131.rds")

## Assemble Fig. 2's multiple panels into final Fig. 2 ----
# ****************************************************

# Theme for figures:
mytheme <- theme_classic() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        plot.subtitle = element_text(size = 15),
        plot.caption = element_text(size = 11))

# Load maps of K_tmp:
fig111spa_50 <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 2/figspa111_50km.rds")
fig121spa_50 <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 2/figspa121_50km.rds")
fig131spa_50 <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 2/figspa131_50km.rds")

# Load maps of dominant mortality agents:
fig111_softmax <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 2/fig_softmax_50km_111.rds")
fig121_softmax <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 2/fig_softmax_50km_121.rds")
fig131_softmax <- readRDS("OUTPUT/figures/figures_manuscript/Fig. 2/fig_softmax_50km_131.rds")

# Prepare each of the subplots:

# Left-hand side plots (maps of K):
fig111spa_50_2 <- (fig111spa_50 +
                     annotate(geom = 'text', x= -97, y = 38, size = 14/.pt,
                              label = expression(~italic("P. elliottii"))) +
                     # labs(tag = "a") +
                     coord_sf(xlim = c(-100, -73),
                              ylim = c(24, 40), expand = FALSE)) +
  labs(tag = "a") + 
  theme(legend.position = c(.9, .22), # tweak position to get legend aligned with legend of corresponding softmax fig.
        plot.tag = element_text(size = 15, face = "bold"), 
        plot.margin = margin(0, 0, 0, 0), axis.title = element_blank())

fig121spa_50_2 <- fig121spa_50 +
  annotate(geom = 'text', x= -97, y = 38, size = 14/.pt,
           label = expression(~italic("P. palustris"))) +
  # labs(tag = "c") +
  coord_sf(xlim = c(-100, -73),
           ylim = c(24, 40), expand = FALSE) +
  guides(colour = "none", fill = "none") +
  labs(tag = "c") + 
  theme(plot.tag = element_text(size = 15, face = "bold"), 
        plot.margin = margin(0, 0, 0, 0), axis.title = element_blank())

fig131spa_50_2 <- fig131spa_50 +
  annotate(geom = 'text', x= -97, y = 38, size = 14/.pt,
           label = expression(~italic("P. taeda"))) +
  # labs(tag = "e") +
  coord_sf(xlim = c(-100, -73),
           ylim = c(24, 40), expand = FALSE) +
  guides(colour = "none", fill = "none") +
  labs(tag = "e") + 
  theme(plot.tag = element_text(size = 15, face = "bold"), 
        plot.margin = margin(0, 0, 0, 0), axis.title = element_blank())

# Right-hand side subplots (mortality agents):
fig111_softmax_2 <- fig111_softmax +
  coord_sf(xlim = c(-100, -73),
           ylim = c(24, 40), expand = FALSE) +
  labs(tag = "b") +
  theme(legend.position = c(.87, .28),
        plot.tag = element_text(size = 15, face = "bold"), 
        plot.margin = margin(0, 0, 0, 0), axis.title = element_blank())

fig121_softmax_2 <- fig121_softmax + 
  coord_sf(xlim = c(-100, -73),
           ylim = c(24, 40), expand = FALSE) +
  guides(fill = "none") +
  labs(tag = "d") +
  theme(plot.tag = element_text(size = 15, face = "bold"), 
        plot.margin = margin(0, 0, 0, 0), axis.title = element_blank())

fig131_softmax_2 <- fig131_softmax +
  coord_sf(xlim = c(-100, -73),
           ylim = c(24, 40), expand = FALSE) +
  guides(fill = "none") +
  labs(tag = "f") +
  theme(plot.tag = element_text(size = 15, face = "bold"), 
        plot.margin = margin(0, 0, 0, 0), axis.title = element_blank())

# Combine plots into a 2x3 grid:
right_column <- (fig111_softmax_2 / fig121_softmax_2 / fig131_softmax_2) & theme(plot.margin = margin(0, 0, 0, 0))
left_column <- (fig111spa_50_2 / fig121spa_50_2 / fig131spa_50_2) & theme(plot.margin = margin(0, 0, 0, 0))

# Create an empty plot with a vertical line only as tall as the combined plot
vertical_line <- ggplot() + 
  geom_vline(xintercept = 0, linetype = "solid", size = 1, color = "black") +
  theme_void() +
  coord_cartesian(ylim = c(0, 1)) +  # Adjust ylim to match the number of rows
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine left and right columns with the vertical line in between
combined_plots <- (left_column | vertical_line | right_column) +
  plot_layout(widths = c(0.49, 0.02, 0.49))

# Add subtitles manually, centered above each column
left_title <- ggplot() + 
  geom_text(aes(x = 0.5, y = 0.5, label = "Species' annual survival variation across their range"), 
            size = 5, fontface = "bold", hjust = 0.5) +
  theme_void() + 
  theme(plot.margin = margin(0, 0, 0, 0))

right_title <- ggplot() + 
  geom_text(aes(x = 0.5, y = 0.5, label = "Spatial variations in the dominance of mortality agents"),
            size = 5, fontface = "bold", hjust = 0.5) +
  theme_void() + 
  theme(plot.margin = margin(0, 0, 0, 0))

# Combine the titles with the plot, ensuring full width for the titles
final_plot <- wrap_plots(
  left_title, right_title, nrow = 1, widths = c(0.49, 0.49)
) / combined_plots +
  plot_layout(heights = c(0.05, 0.95))

# Print the final plot
print(final_plot)

# Save final figure 2:
ggsave(plot = final_plot,
       filename = "OUTPUT/figures/figures_manuscript/Fig. 2/Fig2_merge.pdf", 
       width = 14, height = 14, units = "in")

ggsave(plot = final_plot,
       filename = "OUTPUT/figures/figures_manuscript/Fig. 2/Fig2_merge.png", 
       width = 14, height = 14, units = "in", dpi = 400)

# Note that Figure S5 is generated with the same code as above, but using the grain of 100 x 100 km
# instead of 50 x 50 km. This code is not repeated here to avoid redundancy.


# Test of association between spatially structured K_L variation and cells' dominant mortality agent ----
# **************************************************************************************************

# For Pinus taeda only, since the other species do not present clear spatial patterns in their variation.

## Prepare the data for the model ----
# _______________________________

# Load most likely agent per cell:
sp.id <- 131

load(file = sprintf("OUTPUT/softmax_out/Pred_matrix_%i_50km.Rdata", sp.id))
c_vals <- as.data.frame(t(cell_vals))

c_vals$cell <- row.names(c_vals)

c_vals_long <- c_vals %>% 
  pivot_longer(!cell, names_to = 'agent', values_to = 'prob')

c_vals_main <- c_vals_long %>% 
  group_by(cell) %>% 
  filter(prob == max(prob))

# Merge c_vals_main and data_fig2 (fitted K_tmp against selected MEMs):
# Note that some cells may not have agent information. We drop these (no match):
d_merge <- left_join(rename(c_vals_main, plot = cell), 
                     data_fig2) # created further above, see "Test for spatial structure and select SWM..." in outline

## Create and run the model ----
# ________________________

# Test whether value (the spatially structured variation in K_tmp) is associated to the mortality agents.

# We used an index variable, to avoid using a level of reference as intercept, and to ensure a same 
# number of parameters per level (an intercept per level).

mod_formula <- bf(value ~ 0 + agent)

mod_priors <- c(prior(normal(0, 1), class = b)) # weakly informative priors

d_merge$value_std <- scale(d_merge$value) # standardise the response variable

mod <- brm(formula = mod_formula, 
           family = gaussian(), 
           data = d_merge, 
           prior = mod_priors, 
           iter = 1000, chains = 3, cores = 3, seed = 1)

summary(mod)

## Test hypotheses / contrasts of interest from the model fit ----
# ___________________________________________________________
my_hyp <- c("agentCompetition > agentWeather",
            "agentCompetition > agentInsects",
            "agentCompetition > agentUnknown")

(hyp_test <- hypothesis(mod, 
                        hypothesis = my_hyp, 
                        class = "b", 
                        alpha = 0.1,
                        seed = 2))

### Fig. S7: Create and save the figure ----
# _____________________________________

# plot(hyp_test)

# Make own figure:
data_fig3 <- hyp_test$samples %>% 
  rename("Competition > Weather" = H1,
         "Competition > Insects" = H2,
         "Competition > Unknown" = H3) %>% 
  pivot_longer(cols = 1:3, names_to = "hypoth", values_to = "posterior")

fig_hyp <- data_fig3 %>% 
  ggplot(aes(posterior, hypoth, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  labs(x = "Contrasts' posterior",
       title = expression(atop("Evidence that the spatial patterns of higher mortality of"~italic("P."~taeda)~"are more                       ",
                               "associated with mortality agents Weather, Insects, and Unknown, than with Competition"))) +
  mytheme +
  theme(plot.title = element_text(size = 11, hjust = 0),
        # plot.title.position = "plot",
        axis.title.y = element_blank())

fig_hyp

# Save:
ggsave(fig_hyp,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S7_test_spa_pattern_vs_agent/Fig_S7_hypoth_test_50km.pdf", 
       width = 8.5, height = 5, units = "in")

ggsave(fig_hyp,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S7_test_spa_pattern_vs_agent/Fig_S7_hypoth_test_50km.png", 
       width = 8.5, height = 6, units = "in", dpi = 400)


# Visualize results of the temporal softmax models (Fig. 4, Fig. S3, and Figs. S9-S10) ----
# ************************************************************************************
# ************************************************************************************

# Theme for figures:
mytheme <- theme_bw() +
  theme(panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.5, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"), # 12
        axis.title = element_text(size = 12, colour = "black"), # 14
        legend.text = element_text(size = 10, colour = "black"), # 12
        legend.title = element_text(size = 12, colour = "black"), # 14
        plot.title = element_text(size = 16, colour = "black"),
        plot.subtitle = element_text(size = 15, colour = "black"),
        plot.caption = element_text(size = 11, colour = "black"),
        strip.text = element_text(size = 11, colour = "black"))

# Palette for mortality causes:
mypalette <- c("#4aef94",
               "#c842ab",
               "#b9e055",
               "#004cb8",
               "#ffb84f",
               "#9174b8",
               "#1b4c00")

# Running the code below at the grain of 50 km generates Fig. 4 and Fig. S3.
# Running it for the model outputs based on the grain of 100 km generate Figs. S9 and S10.

## Species 111 (P. elliottii) ----
# ***************************

sp.id <- 111

condeff_111 <- read_csv(file = sprintf("OUTPUT/softmax_out/Cond_effects_year_%i_50km.csv", sp.id))

# # Figure of prob. agent ~ K_tmp or inv_logit(K_tmp):
# # __________________________________________________
# condeff_111 %>% 
#   mutate(K_tmp = boot::inv.logit(K_tmp)) %>% # Silence to keep logit of surv. prob.
#   ggplot(aes(K_tmp, estimate)) +
#   geom_line(aes(colour = group), linewidth = 1) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
#   scale_colour_manual(values = mypalette) +
#   scale_fill_manual(values = mypalette) +
#   scale_x_continuous(breaks = seq(0.88, 1, 0.04)) +
#   # labs(x = "Logit survival prob.",
#   labs(x = "Survival prob.",
#        y = "Prob. of mortality cause") +
#   facet_wrap(~year_0) +
#   mytheme

## Species 121 (P. palustris) ----
# ***************************

sp.id <- 121

condeff_121 <- read_csv(file = sprintf("OUTPUT/softmax_out/Cond_effects_year_%i_50km.csv", sp.id))

# # Figure of prob. agent ~ K_tmp or inv_logit(K_tmp):
# # __________________________________________________
# condeff_121 %>% 
#   mutate(K_tmp = boot::inv.logit(K_tmp)) %>% # Silence to keep logit of surv. prob.
#   ggplot(aes(K_tmp, estimate)) +
#   geom_line(aes(colour = group), linewidth = 1) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
#   scale_colour_manual(values = mypalette) +
#   scale_fill_manual(values = mypalette) +
#   scale_x_continuous(breaks = seq(0.88, 1, 0.04)) +
#   # labs(x = "Logit survival prob.",
#   labs(x = "Survival prob.",
#        y = "Prob. of mortality cause") +
#   facet_wrap(~year_0) +
#   mytheme

## Species 131 (P. taeda) ----
# **********************

sp.id <- 131

condeff_131 <- read_csv(file = sprintf("OUTPUT/softmax_out/Cond_effects_year_%i_50km.csv", sp.id))

# # Figure of prob. agent ~ K_tmp or inv_logit(K_tmp):
# # __________________________________________________
# condeff_131 %>% 
#   mutate(K_tmp = boot::inv.logit(K_tmp)) %>% # Silence to keep logit of surv. prob.
#   ggplot(aes(K_tmp, estimate)) +
#   geom_line(aes(colour = group), linewidth = 1) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
#   scale_colour_manual(values = mypalette) +
#   scale_fill_manual(values = mypalette) +
#   scale_x_continuous(breaks = seq(0.88, 1, 0.04)) +
#   # labs(x = "Logit survival prob.",
#   labs(x = "Survival prob.",
#        y = "Prob. of mortality cause") +
#   facet_wrap(~year_0) +
#   mytheme

## Three species together ----
# ***********************

condeff_allsp <- condeff_111 %>% 
  mutate(sp = "P. elliottii") %>% 
  bind_rows(condeff_121 %>% 
              mutate(sp = "P. palustris")) %>% 
  bind_rows(condeff_131 %>% 
              mutate(sp = "P. taeda")) %>% 
  select(sp, everything()) %>% 
  # Back-transform survival prob. to inverse logit scale:
  mutate(K_tmp = boot::inv.logit(K_tmp))

### Fig. 4: Predictions at two contrasted K_tmps ----
# **********************************************

unique(condeff_allsp$K_tmp)
K_sel <- unique(condeff_allsp$K_tmp)[c(16, 50)] # K  = 0.92, 0.97

fig4 <- condeff_allsp %>% 
  # Filter the two contrasted K_tmp values:
  filter(K_tmp %in% K_sel) %>% 
  mutate(K_tmp = as_factor(round(K_tmp, 2)),
         group = as_factor(group)) %>% 
  # To change the order of the levels of the factor "mortality cause":
  mutate(group = fct_relevel(group, c("Competition", "Weather", "Insects", "Fire", 
                                      "Disease", "Animal", "Unknown"))) %>% 
  # Remove some levels (here "Animals", without interest):
  filter(group %in% c("Competition", "Disease", "Fire", 
                      "Insects", "Weather", "Unknown")) %>% 
  # # In case we'd want "Unknown" hidden to better focus on other agents:
  # filter(group != "Unknown") %>%
  # Figure:
  ggplot(aes(year_0, estimate)) +
  
  # For an barplot (uncertainty not represented):
  # geom_col(aes(fill = K_tmp), position = "dodge") +
  
  # For point-interval plot:
  geom_pointinterval(aes(ymin = conf.low, ymax = conf.high, fill = K_tmp, colour = K_tmp),
                     shape = 21, position = "dodge", linewidth = 1, size = 2) +
  
  scale_fill_manual(name = expression(Under~an~average~survival~"("~italic("K")~")"~of), 
                    values = c("tomato3", "steelblue")) +
  scale_colour_manual(name = expression(Under~an~average~survival~"("~italic("K")~")"~of), 
                      values = c("tomato3", "steelblue")) +
  
  scale_x_continuous(breaks = seq(2003, 2018, 2)) +
  facet_grid(cols = vars(sp), rows = vars(group)) +
  labs(x = "Year", y = "Prob. of mortality agent") +
  mytheme +
  # Remove facet strips and tweak some aspects of the theme:
  theme(strip.background = element_blank(),
        strip.text.x = element_text(face = "italic"),
        # strip.text.y = element_text(size = 10),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank())

fig4

# Save:
ggsave(plot = fig4,
       filename = "OUTPUT/figures/figures_manuscript/Fig. 4/Fig4.pdf", width = 9, height = 7.5, units = "in")  

ggsave(plot = fig4,
       filename = "OUTPUT/figures/figures_manuscript/Fig. 4/Fig4.png", width = 9, height = 7.5, units = "in", dpi = 400)  


### Fig. S3: Facet_grid with sp as cols and years as rows ----
# ********************************************************

fig_causes_K_continuous <- condeff_allsp %>% 
  ggplot(aes(K_tmp, estimate)) +
  geom_line(aes(colour = group), linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3) +
  scale_colour_manual(name = "Agent", values = mypalette) +
  scale_fill_manual(name = "Agent", values = mypalette) +
  scale_x_continuous(breaks = seq(0.88, 1, 0.04)) +
  # Add vertical lines at values of K_L selected for Fig. 4:
  geom_vline(xintercept = 0.92, lty = 2) +
  geom_vline(xintercept = 0.97, lty = 2) +
  # labs(x = "Logit survival prob.",
  labs(x = "Survival prob. (K)",
       y = "Prob. of mortality agent") +
  facet_grid(cols = vars(sp), rows = vars(year_0)) + 
  mytheme +
  # Adapt facet strips:
  theme(strip.background = element_blank(),
        strip.text.x = element_text(face = "italic"),
        strip.text.y = element_text(size = 10),
        legend.position = "bottom")

# Save:
ggsave(plot = fig_causes_K_continuous,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S3_temporal_softmax_continuous_K/FigS3.pdf",
       width = 10, height = 14)

ggsave(plot = fig_causes_K_continuous,
       filename = "OUTPUT/figures/figures_supp_manuscript/Fig_S3_temporal_softmax_continuous_K/FigS3.png",
       width = 10, height = 14, dpi = 400)



