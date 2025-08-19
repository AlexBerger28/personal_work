library(gtrendsR)
library(zctaCrosswalk)
library(tidycensus)
library(readr)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(tidyr)
library(httr)
library(jsonlite)
library(stringdist)
library(tibble)
library(MatchIt)
library(cobalt)


fetch_all_studies <- function(page_size = 1000, delay = 0.5) {
  base_url <- "https://clinicaltrials.gov/api/v2/studies"
  next_token <- NULL
  all_studies <- list()
  
  repeat {
    query_params <- list(pageSize = page_size)
    if (!is.null(next_token)) {
      query_params$pageToken <- next_token
    }
    
    res <- GET(base_url, query = query_params)
    if (status_code(res) != 200) {
      warning("Failed request: ", status_code(res))
      break
    }
    
    content_raw <- content(res, as = "text", encoding = "UTF-8")
    parsed <- fromJSON(content_raw, flatten = TRUE)
    
    if (!"studies" %in% names(parsed)) break
    
    all_studies <- append(all_studies, list(parsed$studies))
    
    if (is.null(parsed$nextPageToken)) break
    next_token <- parsed$nextPageToken
    
    Sys.sleep(delay)  # Be respectful to the server
  }
  
  bind_rows(all_studies)
}

# Run the function to fetch all studies
df <- fetch_all_studies()


NCT00517179 <- df %>% filter(protocolSection.identificationModule.nctId == 'NCT00517179')
search_keyword <- function(df, keyword) {
  matches <- list()
  
  for (col in names(df)) {
    column_data <- df[[col]]
    
    if (is.list(column_data)) {
      # Flatten each element and search inside
      for (i in seq_along(column_data)) {
        entry <- column_data[[i]]
        if (!is.null(entry)) {
          entry_text <- tryCatch(paste(capture.output(str(entry)), collapse = " "), error = function(e) "")
          if (grepl(keyword, entry_text, ignore.case = TRUE)) {
            matches <- append(matches, list(data.frame(row = i, column = col, value = entry_text)))
          }
        }
      }
    } else {
      # Simple text search in regular columns
      match_rows <- which(grepl(keyword, as.character(column_data), ignore.case = TRUE))
      if (length(match_rows) > 0) {
        match_data <- data.frame(row = match_rows, column = col, value = as.character(column_data[match_rows]))
        matches <- append(matches, list(match_data))
      }
    }
  }
  
  if (length(matches) == 0) {
    return(data.frame(row = integer(0), column = character(0), value = character(0)))
  }
  
  do.call(rbind, matches)
}

# Safely extract NCT ID
nct_id <- df$protocolSection.identificationModule.nctId


# Try to extract areas
areas <- tryCatch(
  df$protocolSection.contactsLocationsModule.locations,
  error = function(e) NULL
)

# If arm_groups is NULL or empty, skip
if (is.null(areas) || length(areas) == 0) {
  return(NULL)
}

# Extract labels and return as a data frame
x1 <- tibble(
  nct_id = nct_id,
  label = sapply(areas, function(x) x$zip)
)


results_pre_processed <- df %>% 
  rename(nct_id = protocolSection.identificationModule.nctId,
         title = protocolSection.identificationModule.briefTitle,
         disease_area = protocolSection.conditionsModule.conditions) %>% 
  dplyr::select(nct_id,title,disease_area)


x2 <- x1 %>%
  filter(map_lgl(label, ~ is.character(.x) || is.list(.x))) %>%  # keep only valid entries
  #mutate(label = map(label, as.character)) %>%
  unnest_longer(label) %>%
  left_join(results_pre_processed, by = 'nct_id') %>%
  rename(zcta = label ) %>%
  group_by(zcta) %>% 
  summarise(n=n())  %>%
  mutate(zcta4 = str_sub(zcta, 1, 4)) ##fall back join for university zip codes non residential


county <- read.csv('DMA_FIPS_County_Mapping.csv') %>% 
  mutate(county_fips = sprintf("%05d", FIPS))
acs_vars <- c(
  # Socioeconomic
  med_inc     = "B19013_001",
  pop         = "B01003_001",
  ed_total    = "B15003_001",
  ed_bach     = "B15003_022",
  ed_mast     = "B15003_023",
  ed_prof     = "B15003_024",
  ed_doc      = "B15003_025",
  
  ##income
  hh_total     = "B19001_001",
  hh_100_124   = "B19001_013",
  hh_125_149   = "B19001_014",
  hh_150_199   = "B19001_015",
  hh_200_249   = "B19001_016",
  hh_250_up    = "B19001_017",

  
  uninsured = "B27001_017",
  total_ins = "B27001_001",
  
  # Age 0–17
  m_0_5     = "B27001_003",
  m_6_17    = "B27001_004",
  f_0_5     = "B27001_018",
  f_6_17    = "B27001_019",
  
  # Age 18–64
  m_18_24   = "B27001_005",
  m_25_34   = "B27001_006",
  m_35_44   = "B27001_007",
  m_45_54   = "B27001_008",
  m_55_64   = "B27001_009",
  
  f_18_24   = "B27001_020",
  f_25_34   = "B27001_021",
  f_35_44   = "B27001_022",
  f_45_54   = "B27001_023",
  f_55_64   = "B27001_024",
  
  # Age 65+
  m_65_74   = "B27001_010",
  m_75_up   = "B27001_011",
  f_65_74   = "B27001_025",
  f_75_up   = "B27001_026",
  
  
  # Uninsured children
  unins_m_0_5 = "B27001_005",
  unins_m_6_17 = "B27001_008",
  unins_f_0_5 = "B27001_020",
  unins_f_6_17 = "B27001_023",
  
  # Uninsured adults
  unins_m_18_24 = "B27001_011",
  unins_m_25_34 = "B27001_014",
  unins_m_35_44 = "B27001_017",
  unins_m_45_54 = "B27001_020",
  unins_m_55_64 = "B27001_023",
  
  unins_f_18_24 = "B27001_026",
  unins_f_25_34 = "B27001_029",
  unins_f_35_44 = "B27001_032",
  unins_f_45_54 = "B27001_035",
  unins_f_55_64 = "B27001_038",
  
  # Uninsured 65+
  unins_m_65_74 = "B27001_026",
  unins_m_75_up = "B27001_029",
  unins_f_65_74 = "B27001_041",
  unins_f_75_up = "B27001_044",
  
  total_civ = "B27020_001",
  public_only = "B27020_003",
  public_and_private = "B27020_004",
  
  # Foreign-born
  foreign_born = "B05002_013",
  total_born   = "B05002_001",
  
  # Race
  white_alone  = "B02001_002",
  total_race   = "B02001_001",
  
  # Median Age
  med_age      = "B01002_001",
  
  no_english    = "B16004_005",   # speaks English "not well" or "not at all"
  total_english = "B16004_001",
  
  veterans      = "B21001_002",
  total_vet     = "B21001_001",
  
  disabled      = "B18101_001",  # total civilian noninstitutional pop
  disabled_pop  = "B18101_002",  # total with a disability
  
  broadband_total = "B28002_001",
  broadband_yes   = "B28002_004",
  
  commute_time    = "B08303_001" , # mean travel time to work
  
  # Housing units
  total_units     = "B25024_001",
  sf_detached     = "B25024_003",
  sf_attached     = "B25024_004",
  units_2         = "B25024_005",
  units_3_4       = "B25024_006",
  units_5_9       = "B25024_007",
  units_10_19     = "B25024_008",
  units_20_49     = "B25024_009",
  units_50plus    = "B25024_010",
  
  # Tenure
  occ_units       = "B25003_001",
  owner_occ       = "B25003_002",
  renter_occ      = "B25003_003"
)


county_data <- get_acs(
  geography = "county",
  variables = acs_vars,
  year = 2023,
  survey = "acs5",
  output = "wide",
  cache = TRUE
) %>%
  transmute(
    county_fips    = GEOID,
    med_inc        = med_incE,
    pop            = popE,
    
    # Education
    ed_total       = ed_totalE,
    ed_bachplus    = ed_bachE + ed_mastE + ed_profE + ed_docE,
    #pct_bachplus   = 100 * ed_bachplus / ed_total,
    mean_commute   = commute_timeE,
    
    total_0_17 = m_0_5E + m_6_17E + f_0_5E + f_6_17E,
    total_18_64 = m_18_24E + m_25_34E + m_35_44E + m_45_54E + m_55_64E +
      f_18_24E + f_25_34E + f_35_44E + f_45_54E + f_55_64E,
    total_65_plus = m_65_74E + m_75_upE + f_65_74E + f_75_upE,

    uninsured = uninsuredE,
    total_ins = total_insE,
   # Percent households earning over $100k
   hh_total      = hh_totalE,
   hh_100k_plus  = hh_100_124E + hh_125_149E + hh_150_199E + hh_200_249E + hh_250_upE,

  
    
    # Median Age
    med_age        = med_ageE,
    
    # Middle housing as % of all units
    middle_units = units_2E + units_3_4E + units_5_9E,
    total_units = total_unitsE,
 #   pct_middle_housing = 100 * middle_units / total_unitsE,
    
    # Renter vs. owner
    renter_occ = renter_occE,
    owner_occ = owner_occE,
    occ_units  = occ_unitsE ,
 
 no_english = no_englishE,
 total_english = total_englishE,
 veterans = veteransE,
 total_vet = total_vetE,
 disabled_pop = disabled_popE,
 disabled = disabledE,
 broadband_yes = broadband_yesE,
 broadband_total = broadband_totalE
    
  ) 



dma_level_stats <- county_data %>% 
  left_join(county,by='county_fips') %>% 
  group_by(GOOGLE_DMA) %>% 
  summarize(    pct_owner = sum(owner_occ) / sum(occ_units),
                pct_over100k = sum(hh_100k_plus) / sum(hh_total),
                pct_uninsured = sum(uninsured)/sum(total_ins),
                pct_over65 = sum(total_65_plus)/sum(total_0_17 + total_18_64 + total_65_plus),
               # pct_middle_housing = sum(middle_units) / sum(total_units),
                pct_bachplus = sum(ed_bachplus)/sum(ed_total),
               pct_no_english = sum(no_english)/ sum(total_english),
               pct_veterans   = sum(veterans) / sum(total_vet),
               pct_disabled   = sum(disabled_pop) / sum(disabled),
               pct_broadband  = sum(broadband_yes) / sum(broadband_total),
               mean_commute_hours   = mean(mean_commute/60),
                pop = sum(pop)
                )


crosswalk4 <- zcta_crosswalk %>%
  mutate(zcta4 = str_sub(zcta, 1, 4)) %>%  # fallback for non-res ZIPs
  select(zcta4, county_fips) %>%
  group_by(zcta4) %>%
  slice(1) %>%  # keep only the first row per zcta4 group
  ungroup()

x3 <- x2 %>%  ##zip code level
  left_join(crosswalk4,by='zcta4') %>% ##county_level
  left_join(county,by='county_fips')  ##now add in DMA

x4 <- x3 %>% 
  group_by(GOOGLE_DMA) %>% 
  summarise(dma_count_trials= n())

google_search <- read.csv("healing_hotspots.csv", stringsAsFactors = FALSE) %>% 
  rename(GOOGLE_DMA=geoName) %>% 
  dplyr::select(GOOGLE_DMA,average) %>% 
  rename(healing_search=average)

x5 <- x4 %>% 
  inner_join(google_search,by="GOOGLE_DMA") %>% 
  filter(!is.na(GOOGLE_DMA)) %>% 
  filter(GOOGLE_DMA != c('Jackson TN'))
plot_df <- x5 %>%
  mutate(
    dma_count_trials_z = scale(dma_count_trials),
    healing_search_z = scale(healing_search)
  ) %>% 
  
  pivot_longer(cols = c(dma_count_trials_z, healing_search_z),
               names_to = "metric",
               values_to = "value")


order_by_healing <- x5 %>%
  mutate(healing_search_z = scale(healing_search)) %>%
  arrange(desc(healing_search_z)) %>%
  pull(GOOGLE_DMA) 

ggplot(plot_df, aes(x = factor(GOOGLE_DMA, levels = order_by_healing), y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) ifelse(seq_along(x) %% 2 == 1, x, "")
  ) +
  labs(
    title = "DMAs Sorted by Healing Search Interest (Z-Score)",
    x = "DMA",
    y = "Z-Score"
  ) +
  theme_minimal()


filtered_plot_df <- plot_df %>%
  mutate(row = row_number()) %>%
  filter(row == 1 | row %% 2 == 0)

filtered_plot_df <- plot_df %>%
  mutate(row = row_number()) %>%
  filter(row == 1 | row %% 2 == 1)

ggplot(filtered_plot_df, aes(x = factor(GOOGLE_DMA, levels = order_by_healing), y = value, fill = metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_x_discrete(
    labels = function(x) ifelse(x %in% filtered_plot_df$GOOGLE_DMA, x, "")
  ) +
  labs(
    title = "DMAs Sorted by Healing Search Interest (Z-Score)",
    x = "DMA",
    y = "Z-Score"
  ) +
  theme_minimal()



x6 <- x5 %>% 
  left_join(dma_level_stats,by='GOOGLE_DMA') %>% 
  mutate(trials = 1000*dma_count_trials/pop)
cor_df <- x6 %>%
  select(healing_search, trials, pct_owner, pct_over100k,
         pct_bachplus, pop) %>%
  cor(use = "complete.obs")

print(cor_df)
model <- lm(healing_search ~ 
              pct_over100k + 
              pct_over65 +
              trials +
              pct_owner  + pct_uninsured + pct_bachplus + pct_no_english + 
              pct_veterans + pct_disabled + pct_broadband + mean_commute_hours,
            data = x6)  

summary(model)




x7 <- x6 %>%
  filter(trials <  quantile(trials, 0.95, na.rm = TRUE) ) %>% 
  mutate(treated = ifelse(trials >= quantile(trials, 0.5, na.rm = TRUE), 1, 0)) %>%
  drop_na(treated, pct_over100k, pct_over65, pct_owner,
           pct_uninsured, pct_bachplus, healing_search)

ps_model <- matchit(
  treated ~ pct_over100k + 
    pct_over65 +
    pct_owner  + pct_uninsured + pct_bachplus + pct_no_english + 
    pct_veterans + pct_disabled + pct_broadband + mean_commute_hours + pop,
  data = x7,
  method = "nearest",
  distance = "logit",
  caliper = .1,             # remove poor matches
  discard = "both"           # remove units with non-overlapping propensity scores
)

summary(ps_model)


love.plot(ps_model, binary = "std", 
          var.order = "unadjusted", 
          stat = "mean.diffs", 
          threshold = 0.2)

matched_x7 <- match.data(ps_model)

lm_out <- lm(healing_search ~ trials, data = matched_x7)
summary(lm_out) ##look at number of trials 
lm_out <- lm(healing_search ~ treated, data = matched_x7) 
summary(lm_out) ## just split high and low regions of trials

##distribution of trials
hist(matched_x7$trials, breaks = 20)



mod_interact <- lm(healing_search ~ trials * pct_over100k +
                     trials * pct_bachplus +
                     trials * pct_veterans +
                     trials * pct_over65 + 
                     trials * pct_veterans + 
                     trials * pct_broadband,
                   data = matched_x7)

summary(mod_interact)
library(interactions)

interact_plot(mod_interact,
              pred = trials,
              modx = pct_bachplus,
              interval = TRUE,
              plot.points = TRUE,
             x.label = c('Trials per 1000 People'),
             y.label = c('Healing Search Score'),
             legend.main = c('% with At Least Bachelors'))



hist(matched_x7$pct_over65, breaks = 20, main = "Distribution of % over 65", xlab = "% 65")
abline(v = median(matched_x7$pct_over65, na.rm = TRUE), col = "red", lty = 2)
quantile(matched_x7$pct_over65, probs = seq(0, 1, 0.1), na.rm = TRUE)




# Only look at matched subclasses with both treated and control
matched_subclasses <- matched_x7 %>%
  group_by(subclass) %>%
  filter(any(treated == 1) & any(treated == 0))

# Summarize by metro and treatment
summary_by_dma <- matched_subclasses %>%
  group_by(subclass, GOOGLE_DMA, treated) %>%
  summarize(
    pct_over65 = mean(pct_over65),
    pct_over100k = mean(pct_over100k),
    pct_bachplus = mean(pct_bachplus),
    healing_search = mean(healing_search),
    .groups = 'drop'
  )

# Join treated and control rows within the same subclass
joined <- summary_by_dma%>%
  pivot_wider(names_from = treated, values_from = c(pct_over65, pct_over100k, pct_bachplus, healing_search, GOOGLE_DMA), names_prefix = "treated_") %>%
  filter(
    abs(pct_over65_treated_1 - pct_over65_treated_0) < 1,
    abs(pct_over100k_treated_1 - pct_over100k_treated_0) < 1,
    abs(pct_bachplus_treated_1 - pct_bachplus_treated_0) < 1,
    healing_search_treated_1 < healing_search_treated_0  # low in treated, high in control
  )


##that was too subclass specific


dma_ranked <- matched_x7 %>%
  mutate(
    rank_over65 = rank(pct_over65),
    rank_over100k = rank(pct_over100k),
    rank_bachplus = rank(pct_bachplus)
  ) %>%
  rowwise() %>%
  mutate(rank_score = mean(c(rank_over65, rank_over100k, rank_bachplus))) %>%
  ungroup()
treated <- dma_ranked %>% filter(treated == 1)
control <- dma_ranked %>% filter(treated == 0)

# Do a fuzzy join to find nearest control for each treated
library(fuzzyjoin)

matched_pairs <- dma_ranked %>%
  filter(!is.na(rank_score)) %>%
  select(GOOGLE_DMA, treated, healing_search, pop, pct_over65, pct_over100k, pct_bachplus, rank_score, trials) %>%
  fuzzy_inner_join(
    dma_ranked %>% filter(treated == 0),
    by = c("rank_score" = "rank_score"),
    match_fun = function(x, y) abs(x - y) <= 2
  ) %>%
  filter(treated.x == 1 & treated.y == 0) %>%
  select(rank_score.x,rank_score.y,
         trials.x,trials.y,
    GOOGLE_DMA.x,  GOOGLE_DMA.y,
    pop.x,pop.y,
    pct_over65.x, pct_over65.y,
    pct_over100k.x, pct_over100k.y,
    pct_bachplus.x, pct_bachplus.y,
    healing_search.x, healing_search.y
  ) %>%
  arrange(abs(healing_search.x - healing_search.y)) %>% 
  mutate(diff = healing_search.y - healing_search.x) %>% 
  filter(pop.x/pop.y < 1.5 & pop.x/pop.y > .5)



# Create a label to identify each pair
matched_pairs <- matched_pairs %>%
  mutate(
    dma_pair = paste0(GOOGLE_DMA.x, " vs. ", GOOGLE_DMA.y)
  ) %>%
  arrange(desc(diff)) %>%
  mutate(dma_pair = factor(dma_pair, levels = unique(dma_pair))) %>% 
  filter(abs(diff) > 3)

ggplot(matched_pairs, aes(y = dma_pair)) +

  geom_segment(aes(x = healing_search.x, xend = healing_search.y, y = dma_pair, yend = dma_pair), 
               color = "gray70", size = 1) +
  
  geom_point(aes(x = healing_search.x), color = "blue", size = 3) +
  geom_point(aes(x = healing_search.y), color = "red", size = 3) +
  
  geom_tile(aes(x = max(healing_search.x, healing_search.y) + 8, fill = pct_bachplus.x), 
            width = 4, height = 0.6) +
  
  scale_fill_gradient(low = "lightyellow", high = "darkgreen", name = "% Bach+ (Treated DMA)") +
  
  labs(
    title = "Healing Search Gaps Between Matched Metro Pairs",
    subtitle = "Sorted by size of healing search gap (Control – Treated).\nRight color bar = % with bachelor’s degree (Treated DMA)",
    x = "Healing Search Score",
    y = "",
    caption = "Matched on age, income, education, and population size.\nBlue = Treated (more trials), Red = Control"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  ) +
  xlim(NA, max(matched_pairs$healing_search.x, matched_pairs$healing_search.y, na.rm = TRUE) + 15)

matched_pairs <- matched_pairs %>%
  mutate(
    dma_pair = paste0(GOOGLE_DMA.x, " vs. ", GOOGLE_DMA.y),
    dma_pair = factor(dma_pair, levels = dma_pair[order(-diff)]),
    trials_label_x =  round(trials.x, 2),
    trials_label_y = round(trials.y, 2)
  )

# Plot
ggplot(matched_pairs, aes(y = dma_pair)) +
  # Line connecting healing search scores
  geom_segment(aes(x = healing_search.x, xend = healing_search.y, yend = dma_pair), 
               color = "gray70", size = 1) +
  
  # Healing search score dots
  geom_point(aes(x = healing_search.x, color = "Treated"), size = 3) +
  geom_point(aes(x = healing_search.y, color = "Control"), size = 3) +
  
  # Trial labels above each dot
  geom_text(aes(x = healing_search.x, label = trials_label_x), 
            vjust = -1, hjust = 0, size = 3.2, color = "blue") +
  geom_text(aes(x = healing_search.y, label = trials_label_y), 
            vjust = -1, hjust = 1, size = 3.2, color = "red") +
  
  # Education gradient strip
  geom_tile(aes(x = max(healing_search.x, healing_search.y, na.rm = TRUE) + 8, fill = pct_bachplus.x), 
            width = 4, height = 0.6) +
  
  # Color legends
  scale_color_manual(
    name = "DMA Type",
    values = c("Treated" = "blue", "Control" = "red")
  ) +
  scale_fill_gradient(
    low = "lightyellow", high = "darkgreen",
    name = "% Bachelor’s (Treated DMA)"
  ) +
  
  labs(
    title = "Healing Search Score Gaps Between Matched Metro Pairs",
    x = "Healing Search Score",
    y = "",
    caption = "Matched on age, income, education, and population size"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right"
  ) +
  xlim(NA, max(matched_pairs$healing_search.x, matched_pairs$healing_search.y, na.rm = TRUE) + 20)

x8 <- matched_x7 %>%
  mutate(income_group = ifelse(pct_over100k >= median(pct_over100k,na.rm=TRUE), "high", "low")) %>% 
  mutate(uninsured_group = ifelse(pct_uninsured >= median(pct_uninsured, na.rm = TRUE), "high", "low"))
mod_simple <- lm(healing_search ~ trials*pct_uninsured , data = x8)
summary(mod_simple)



abbr_n  <- 30     # first N letters for point labels
padding <- 20    # extra space to the right of 70 for labels/strip

xr <- range(c(matched_pairs$healing_search.x, matched_pairs$healing_search.y), na.rm = TRUE)
label_pad <- diff(xr) * 0.015 + 0.25

mp <- matched_pairs %>%
  mutate(
    diff = healing_search.x - healing_search.y,
    pair_label = sprintf("%s vs %s", GOOGLE_DMA.x, GOOGLE_DMA.y),
    # point labels (abbreviated)
    GOOGLE_DMA_abbrev.x = substr(GOOGLE_DMA.x, 1, abbr_n),
    GOOGLE_DMA_abbrev.y = substr(GOOGLE_DMA.y, 1, abbr_n),
    # which side is High (blue)
    high_on_left = healing_search.x <= healing_search.y,
    # precompute label positions & justification so they don't overlap the dots
    x_high_label = healing_search.x + ifelse(high_on_left, -label_pad,  label_pad),
    x_low_label  = healing_search.y + ifelse(high_on_left,  label_pad, -label_pad),
    high_hjust   = ifelse(high_on_left, 1, 0),
    low_hjust    = ifelse(high_on_left, 0, 1)
  ) %>%
  arrange(diff) %>%
  mutate(pair_label = factor(pair_label, levels = pair_label))

ggplot(mp, aes(y = pair_label)) +
  # dumbbell
  geom_segment(aes(x = healing_search.x, xend = healing_search.y, yend = pair_label),
               color = "gray70", size = 1) +
  geom_point(aes(x = healing_search.x, color = "High"), size = 3) +
  geom_point(aes(x = healing_search.y, color = "Low"),  size = 3) +
  
  # DATA LABELS over points (abbreviated, auto-flipped)
  geom_text(aes(x = x_high_label, label = GOOGLE_DMA_abbrev.x, color = "High", hjust = high_hjust),
            vjust = -0.25, size = 3.2, show.legend = FALSE) +
  geom_text(aes(x = x_low_label,  label = GOOGLE_DMA_abbrev.y,  color = "Low",  hjust = low_hjust),
            vjust = -0.25, size = 3.2, show.legend = FALSE) +
  
  # single vertical BA+ strip
  geom_tile(aes(x = 70 + padding/2, fill = pct_bachplus.x),
            width = 4, height = 0.6) +
  
  scale_color_manual(
    name = "GOOGLE_DMA Trial Availability",
    values = c("High" = "blue", "Low" = "red"),
    limits = c("High", "Low")
  ) +
  scale_fill_gradient(low = "lightyellow", high = "darkgreen",
                      name = "% Bachelor’s (High DMA)") +
  scale_x_continuous(limits = c(0, 70 + padding), breaks = seq(0, 70, 10),
                     expand = expansion(mult = c(0, 0))) +
  
  labs(
    title = "Healing Search Score Gaps Between Matched Metro Pairs",
    subtitle = "Sorted by signed High – Low difference; DMA labels on points",
    x = "Healing Search Score", y = NULL,
    caption = "Matched on age, income, education, and population size"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    axis.text.y = element_blank(),   # hide y-axis labels
    axis.ticks.y = element_blank(),  # hide y ticks
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  )
