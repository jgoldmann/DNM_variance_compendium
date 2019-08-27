

library(tidyverse)
library(drake)
library(here)
library(rtracklayer)
library(rptR)
library(inovaDNMs)
library(inovaDNMs102)
library(dnmVarianceComponents)
library(ggbeeswarm)
library(foreach)
library(poissonMutSim)
library(car)


inova1 <- 
  inovaDNMs::obs %>%
  dplyr::rename(
    SampleID = internalID,
    Variant = Variant, 
    Chr = Chromosome,
    Position = Start.position)  %>%
  full_join(inovaDNMs::trios %>% 
              filter(in.data), 
            by=c(SampleID = "Trio.ID")) %>%
  dplyr::rename(familyNr = Family.Study.ID)

inova2 <- 
  inovaDNMs102::allDNMs %>%
  dplyr::rename(
    SampleID = sampleName,
    Variant = alt, 
    Chr = chrom,
    Position = start) %>%
  full_join(inovaDNMs102::trios, 
            by=c(SampleID = "PROBANDID")) %>%
  dplyr::rename(familyNr = FAMILYID) %>%
  mutate(end = NULL)

#sasani DNMs
sasani <- 
  read_delim("https://github.com/quinlan-lab/ceph-dnm-manuscript/raw/master/data/second_gen.dnms.txt",
             delim = "\t",
             col_types = cols(chrom = col_character(),
                              new_family_id = col_character())) %>%
  bind_rows(read_delim("https://github.com/quinlan-lab/ceph-dnm-manuscript/raw/master/data/third_gen.dnms.txt",
                       delim = "\t",
                       col_types = cols(chrom = col_character(),
                                        new_family_id = col_character()))) %>%
  bind_rows(read_delim("https://github.com/quinlan-lab/ceph-dnm-manuscript/raw/master/data/post-pgcs.dnms.txt",
                       delim = "\t",
                       col_types = cols(chrom = col_character(),
                                        new_family_id = col_character()))) %>%
  bind_rows(read_delim("https://github.com/quinlan-lab/ceph-dnm-manuscript/raw/master/data/gonosomal.dnms.txt",
                       delim = "\t",
                       col_types = cols(chrom = col_character(),
                                        new_family_id = col_character()))) %>%
  mutate(Chr = paste0("chr", chrom)) %>%
  dplyr::rename(
    SampleID = new_sample_id,
    Variant = ref, 
    Position = start,
    fathersAgeAtConceptionInYears = paternal_age_at_birth,
    mothersAgeAtConceptionInYears = maternal_age_at_birth,
    familyNr = new_family_id
  )

autism1 <- dnmVarianceComponents::dnmDataAut
autism2 <- dnmVarianceComponents::dnmDataAut2


# --------------


#comparison of observed family residual sum of squares to bootstrapped values

simpMod <- lm(unique_snv ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears, autism1)

obs <- 
  autism1 %>%
  mutate(resi = residuals(simpMod)) %>%
  lm(resi~familyNr, data=.) %>% 
  deviance()

sim <- 
  foreach(1:10000,
          .combine = c) %do% {
  autism1 %>%
    mutate(resi = sample(residuals(simpMod))) %>%
    lm(resi~familyNr, data=.) %>% 
              deviance()
          }

boxplot(sim)
points(y=obs,x=1, pch = 18, col = "red", cex=2)

ggplot(data.frame(simulations=sim),
       aes(x=sim)) + 
  geom_density() + 
  geom_vline(xintercept = obs,
             col = "red")

(table(obs <= sim)/10000)["FALSE"]  #81% of bootstrapped results have a lower deviance



compareFamilyRSS <- function(autism1) {
  simpMod <- lm(n ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears, autism1)
  obs <- 
    autism1 %>%
    mutate(resi = residuals(simpMod)) %>%
    lm(resi~familyNr, data=.) %>% 
    deviance()
  sim <- 
    foreach(1:10000,
            .combine = c) %do% {
              autism1 %>%
                mutate(resi = sample(residuals(simpMod))) %>%
                lm(resi~familyNr, data=.) %>% 
                deviance()
            }
  (table(obs <= sim)/10000)["FALSE"]
}

filterMuts(inova1, uniqGenome) %>% 
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr) %>% 
  count() %>%
  ungroup() %>% 
  compareFamilyRSS()  #0.0871

filterMuts(inova2, uniqGenome) %>% 
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr) %>% 
  count() %>%
  ungroup() %>% 
  compareFamilyRSS()

filterMuts(sasani %>% dplyr::select(-end), uniqGenome) %>% 
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr) %>% 
  count() %>%
  ungroup() %>% 
  compareFamilyRSS() #0.1286

compareFamilyRSS(autism1 %>% mutate(n=unique_snv))
compareFamilyRSS(autism2 %>% mutate(n=snvs_unique))
