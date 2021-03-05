library(RTNsurvival)
library(tidyverse)
library(survival)
library(ranger)
library(ggfortify)
library(gdata)

theme_set(theme_light())

load("../rdatas/g3_survival.RData")

interest_regs <- c("BHLHE41", "CAMTA1", "ZNF365", "KCNIP3", "RFX4", "SOX2", 
                    "NACC2", "ZNF385B", "NR1D1", "LHX4")

#----

enrich_scores <- tnsGet(rtns, "regulonActivity")

status <- as.data.frame(enrich_scores$status)
enrich_scores <- as.data.frame(enrich_scores$differential)

survival_data <- tnsGet(rtns, "survivalData")


#---- Just to make sample names more readable ----

rownames(enrich_scores) <- str_remove(rownames(enrich_scores), ".CEL")
rownames(survival_data) <- str_remove(rownames(survival_data), ".CEL")
rownames(status) <- str_remove(rownames(status), ".CEL")


#----

interest_status <- status %>% 
  select(interest_regs) %>% 
  rownames_to_column(var = "sample")
  
status_and_survival <- survival_data %>% 
  rownames_to_column(var = "sample") %>% 
  select(sample, time, event) %>% 
  inner_join(interest_status)


#----

complete_data %>% 
  select(sample, time, event, SOX2) %>% 
  mutate()



#----

km <- with(status_and_survival, Surv(time, event))

km_fit <- survfit(Surv(time, event) ~ CAMTA1, data = status_and_survival)

summary(km_fit, times = c(1,30,60,90*(1:10)))

# km_fit$n.censor[km_fit$n.censor == 0] <- NA_integer_
# 
# ggplot(km_fit) +
#   geom_step(aes(time, surv, color = strata)) +
#   geom_point(aes(time*n.censor, surv*n.censor), size = 1, shape = 3, color = "black") +
#   scale_y_continuous(limits = c(0, 1), labels = scales::percent)
  
plot <- 
autoplot(km_fit, surv.linetype = "solid", conf.int = FALSE,
         censor.size = 3, surv.size = 0.5,
         main = "Kaplan Meier Estimate",
         xlab = "Time (Years)",
         ylab = "Survival Rate") 

plot + labs(color = "Title")


