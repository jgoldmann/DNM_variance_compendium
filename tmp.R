
bli <- 
  inova2 %>%
  filterMuts(uniqGenome) %>%
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr) %>% 
  count() %>% 
  lme4::lmer(n ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears + (1|familyNr),
                  data = .) %>% 
  lme4::ranef() %>% 
  pluck(1) %>% 
  rownames_to_column() %>% 
  dplyr::rename(withAll = `(Intercept)`)

bla <- 
  inovaDNMs102::allDNMs %>%
  dplyr::rename(
    SampleID = sampleName,
    Variant = alt, 
    Chr = chrom,
    Position = start) %>%
  full_join(inovaDNMs102::trios, 
            by=c(SampleID = "PROBANDID")) %>%
  dplyr::rename(familyNr = FAMILYID) %>%
  mutate(end = NULL)%>%
  filterMuts(uniqGenome) %>%
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr) %>% 
  count() %>% 
  lme4::lmer(n ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears + (1|familyNr),
             data = .) %>% 
  lme4::ranef() %>% 
  pluck(1) %>% 
  rownames_to_column() %>% 
  dplyr::rename(NoSibGroups = `(Intercept)`)

full_join(bla, bli,
          by = "rowname") %>% 
  ggplot(aes(x=NoSibGroups,y=withAll)) + 
  geom_point()  ## something weird happens to the fit here
