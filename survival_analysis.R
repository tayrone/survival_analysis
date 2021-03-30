library(RTN)
library(RTNsurvival)
library(tidyverse)
library(survival)
library(ranger)
library(ggfortify)
library(gdata)
library(survminer)

theme_set(theme_light())

load("../rdatas/g3_survival.RData")

interest_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")


#---- Gets required data for analysis ----

enrich_scores <- tnsGet(rtns, "regulonActivity")

enrich_scores <- as.data.frame(enrich_scores$differential)

survival_data <- tnsGet(rtns, "survivalData")


#---- Just to make sample names more readable ----

rownames(enrich_scores) <- str_remove(rownames(enrich_scores), ".CEL")
rownames(survival_data) <- str_remove(rownames(survival_data), ".CEL")


#---- Some data wrangling ----

score_stats <- enrich_scores %>% 
  select(all_of(interest_regs)) %>% 
  summarise_all(list("mean" = mean, "median" = median, "iqr" = IQR))  %>% 
  gather("method", "value") %>% 
  separate(method, c("regulon", "method"), sep = "_") %>% 
  spread(regulon, value)

score_stats <- map(transpose(score_stats, .names = score_stats$method), 
                   as.vector)


calculate_status <- function(regulon, stat){
  enrich_scores %>% 
    select(all_of(regulon)) %>% 
    rownames_to_column(var = "sample") %>% 
    mutate(stat_value = score_stats[[stat]][[regulon]],
           status = case_when(get(regulon) >= stat_value ~ 1,
                              get(regulon) < stat_value ~ -1))
}

interest_regs <- set_names(interest_regs, interest_regs)

median_status <- map(interest_regs, calculate_status, stat = "median")

mean_status <- map(interest_regs, calculate_status, stat = "mean")

iqr_status <- map(interest_regs, calculate_status, stat = "iqr")


#---- Calculate regulon status considering different statistical measures ----


add_survival <- function(status){
  
status_and_survival <- 
  survival_data %>%
    rownames_to_column(var = "sample") %>%
    select(sample, time, event) %>%
    inner_join(status)

}


median_status <- map(median_status, add_survival)
mean_status <- map(mean_status, add_survival)
iqr_status <- map(iqr_status, add_survival)


#---- Generates the KM estimation and subsequent plots ----

status <- list("mean" = mean_status, 
               "median" = median_status, 
               "iqr" = iqr_status)

for(stat in names(status)){

  for(regulon in interest_regs){
  
    regulon_data <- status[[stat]][[regulon]]
    
    
    km <- with(regulon_data, Surv(time = time, event = event))
    
    km_fit <- survfit(Surv(time = time, event = event) ~ status, 
                      data = regulon_data)
    
    surv_pvalue(km_fit)$pval.txt
    surv_median(km_fit)
    print(km_fit)
    
    summary(km_fit, times = c(1,30,60,90*(1:10)))
    
    
    plot <- autoplot(km_fit, surv.linetype = "solid", conf.int = FALSE,
                     censor.size = 3, surv.size = 0.5, main = "Kaplan Meier Estimate", 
                     xlab = "Time (Years)", ylab = "Survival Rate") 
    
    plot + labs(color = "Regulon Activity Status") +
      scale_color_discrete(labels = c("Up", "Down")) +
      labs(tag = paste0("Log-rank p-value: ", round(surv_pvalue(km_fit)$pval, 3))) +
      theme(plot.tag = element_text(size = 10, face = "plain"), 
            plot.tag.position = c(.803, .985))
    
    
    ggsave(paste0("./km_plots/", stat, "/", 
                 regulon, ".png"))
  
  }

}
