library(tidyverse)

alrtd_p = 0.3
pn_given_alrtd_p = 0.1
pn_given_no_alrtd_p = 0.001
sample_size = 10000

sim = tibble(
  id = 1:sample_size,
  alrtd_status = rbinom(n = sample_size, size = 1, prob = alrtd_p),
  pn_status = rbinom(n = sample_size, size = 1, prob = alrtd_status*pn_given_alrtd_p+(1-alrtd_status)*pn_given_no_alrtd_p)
)

symptom_screening_sens = 0.95
symptom_screening_spec = 0.3

case_review_sens = 0.95
case_review_spec = 0.99

urine_collected_sens = 0.5
urine_collected_spec = 0.5

uad_sens = 0.8
uad_spec = 0.99

sim = sim %>% mutate(
  
  symptom_screening_result = rbinom(n = sample_size, size = 1, alrtd_status * symptom_screening_sens + (1-alrtd_status)*(1-symptom_screening_spec)),
  
  case_review_result = ifelse(
    symptom_screening_result == 0 , NA, 
    rbinom(n = sample_size, size = 1, alrtd_status * case_review_sens + (1-alrtd_status)*(1-case_review_spec)))
  ,
      
  urine_collected_result = ifelse(
    case_review_result == 0, NA, 
    rbinom(n = sample_size, size = 1, alrtd_status * urine_collected_sens + (1-alrtd_status)*(1-urine_collected_spec))
  ),  
  
  uad_pn_result = ifelse(
    urine_collected_result == 0, NA, 
    rbinom(n = sample_size, size = 1, pn_status * uad_sens + (1-pn_status)*(1-uad_spec))
  )
  
)

sim %>% summarise(
  across(c(-id), ~ sum(.x,na.rm = TRUE)/ sum(!is.na(.x)))
)

sim %>% summarise(
  across(c(-id), ~ sum(alrtd_status * !is.na(.x),na.rm = TRUE)/ sum(!is.na(.x)))
)

# prior probability
sim %>% summarise(
  across(c(-id), ~ sum(pn_status * !is.na(.x),na.rm = TRUE)/ sum(!is.na(.x)))
)

# posterior probability
sim %>% summarise(
  across(c(-id), ~ sum(pn_status * (!is.na(.x) & .x==1))/ sum(!is.na(.x) & .x==1))
)

sim %>% with(table(pn_status, uad_pn_result,useNA = "always"))

sim %>% summarise(
  p_uad_positive = sum(uad_pn_result,na.rm=TRUE) / sum(!is.na(uad_pn_result)),
  p_pn_given_no_uad = sum(pn_status & is.na(uad_pn_result)) / sum(is.na(uad_pn_result))
)


# alluvial

sim_all = sim %>% 
  mutate(across(c(-id), as.logical)) %>%
  group_by(alrtd_status,pn_status,symptom_screening_result,case_review_result,urine_collected_result,uad_pn_result) %>% summarise(freq = n())

ggplot(sim_all, aes(axis1 = pn_status, axis2 = alrtd_status, axis3 = symptom_screening_result, axis4 = case_review_result, axis5 = urine_collected_result, axis6 = uad_pn_result, y = freq)) +
  ggalluvial::geom_alluvium(aes(fill = pn_status))+
  ggalluvial::geom_stratum()

sim_san = sim %>% ggsankey::make_long(pn_status,alrtd_status,symptom_screening_result,case_review_result,urine_collected_result,uad_pn_result)
sim_san = sim %>% pivot_longer(cols = )

sim_san = sim_san %>% filter(!is.na(node) & !is.na(next_node))

ggplot(sim_san, aes(x = x, 
                      next_x = next_x, 
                      node = node, 
                      next_node = next_node,
                      fill = factor(node))) +
  ggsankey::geom_sankey()
  
