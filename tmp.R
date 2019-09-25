

#2 models are weirdly different

mdl_withSibs <- 
  inova2 %>%
  filterMuts(uniqGenome) %>%
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr,
           batch) %>% 
  count() %>% 
  lme4::lmer(n ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears + (1|familyNr) + (1|batch),
             data = .)


mdl_noSibs <- 
  inovaDNMs102::allDNMs %>%
  dplyr::rename(
    SampleID = sampleName,
    Variant = alt, 
    Chr = chrom,
    Position = start) %>%
  full_join(inovaDNMs102::trios, 
            by=c(SampleID = "PROBANDID")) %>%
  dplyr::rename(familyNr = FAMILYID) %>% 
  left_join(inovaDNMs102::qcStats %>%
              dplyr::select(pipeline, STUDY_ID) %>%
              mutate(STUDY_ID=toupper(STUDY_ID), batch=pipeline),
            by=c(SampleID="STUDY_ID")) %>%
  mutate(end = NULL)%>%
  filterMuts(uniqGenome) %>%
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr,
           batch) %>% 
  count() %>% 
  lme4::lmer(n ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears + (1|familyNr) + (1|batch),
             data = .)



mdl_withSibs %>% broom::tidy()
mdl_withSibs %>% broom::glance()

mdl_noSibs %>% broom::tidy()
mdl_noSibs %>% broom::glance()





bli <- 
  mdl_withSibs %>% 
  lme4::ranef() %>% 
  pluck(1) %>% 
  rownames_to_column() %>% 
  dplyr::rename(withAll = `(Intercept)`)

bla <- 
  mdl_noSibs %>% 
  lme4::ranef() %>% 
  pluck(1) %>% 
  rownames_to_column() %>% 
  dplyr::rename(NoSibGroups = `(Intercept)`)

full_join(bla, bli,
          by = "rowname") %>% 
  ggplot(aes(x=NoSibGroups,y=withAll)) + 
  geom_point()  
## something weird happens to the fit here



mdl_noSibs %>% broom::augment() %>% 
  full_join(mdl_withSibs %>% broom::augment(),
            by=c("n", "fathersAgeAtConceptionInYears", "mothersAgeAtConceptionInYears", "batch", "familyNr"))






