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

status <- as.data.frame(enrich_scores$status)
enrich_scores <- as.data.frame(enrich_scores$differential)

survival_data <- tnsGet(rtns, "survivalData")


#---- Just to make sample names more readable ----

rownames(enrich_scores) <- str_remove(rownames(enrich_scores), ".CEL")
rownames(survival_data) <- str_remove(rownames(survival_data), ".CEL")
rownames(status) <- str_remove(rownames(status), ".CEL")


#---- Some data wrangling ----

score_stats <- enrich_scores %>% 
  select(all_of(interest_regs)) %>% 
  summarise_all(list("mean"= mean, "median"= median, "iqr" = IQR))  %>% 
  gather("method", "value") %>% 
  separate(method, c("regulon", "method"), sep = "_") %>% 
  spread(regulon, value)

score_stats <- map(transpose(score_stats, .names = score_stats$method), 
                   as.vector)

enrich_scores %>% 
  select(CAMTA1) %>% 
  mutate(iqr_value = score_stats[["median"]][["CAMTA1"]],
         status = case_when(CAMTA1 > iqr_value ~ 1,
                            near(CAMTA1, iqr_value) ~ 0,
                            CAMTA1 < iqr_value ~ -1))

calculate_status <- function(regulon_scores, stat){
  print(names(regulon_scores))
  
}

map(enrich_scores, calculate_status, stat = stats_median)



# interest_status <- status %>% 
#   select(all_of(interest_regs)) %>% 
#   rownames_to_column(var = "sample")
#   
# status_and_survival <- survival_data %>% 
#   rownames_to_column(var = "sample") %>% 
#   select(sample, time, event) %>% 
#   inner_join(interest_status)


#---- Generates the KM estimation and subsequent plots ----

for(regulon in interest_regs){

  km <- with(status_and_survival, Surv(time = time, event = event))
  
  km_fit <- survfit(Surv(time = time, event = event) ~ get(regulon), 
                    data = status_and_survival)
  
  surv_pvalue(km_fit)$pval.txt
  surv_median(km_fit)
  print(km_fit)
  
  summary(km_fit, times = c(1,30,60,90*(1:10)))
  
    
  plot <- autoplot(km_fit, surv.linetype = "solid", conf.int = FALSE,
          censor.size = 3, surv.size = 0.5, main = "Kaplan Meier Estimate", 
          xlab = "Time (Years)", ylab = "Survival Rate") 
  
  plot + labs(color = "Regulon Activity Status") +
    scale_color_discrete(labels = c("Up", "Neutral", "Down")) +
    labs(tag = paste0("Log-rank p-value: ", round(surv_pvalue(km_fit)$pval, 3))) +
    theme(plot.tag = element_text(size = 10, face = "plain"), 
          plot.tag.position = c(.803, .985))
  
  ggsave(paste0("./km_plots/", regulon, ".png"))

}

  
