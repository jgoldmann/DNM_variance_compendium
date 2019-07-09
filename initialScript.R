library(tidyverse)
library(drake)
library(here)
library(rtracklayer)
library(rptR)


#gather datasets
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
sasani <- 
  DNMsets::DNMs_Sas %>%
  dplyr::rename(
    SampleID = Trio.ID,
    Variant = Variant, 
    Chr = Chromosome,
    Position = Start.position
  )


autism1 <- dnmVarianceComponents::dnmDataAut
autism2 <- dnmVarianceComponents::dnmDataAut2

#unique genome fraction
stopifnot(file.exists(here("data/master_repeat_file.bed.gz")))
stopifnot(file.exists(here("data/recent_repeat_file.bed.gz")))

uniqGenome <- 
  import(here("data/master_repeat_file.bed.gz")) %>% 
  gaps()
seqlevelsStyle(uniqGenome) <- "UCSC"
recRepGenome <- 
  import(here("data/recent_repeat_file.bed.gz"))
seqlevelsStyle(uniqGenome) <- "UCSC"


#functions for cleaning data
filterMuts <- 
  function(muts, ## needs columns "SampleID", Variant, Chr, Position
           filter){
    muts %>% 
      makeGRangesFromDataFrame(seqnames.field = "Chr", 
                               start.field = "Position", 
                               end.field = "Position", 
                               ignore.strand = TRUE, 
                               keep.extra.columns = TRUE) %>% 
      subsetByOverlaps(filter) %>%
      as.data.frame() %>%
      as.tibble() %>%
      mutate(Chr=seqnames, Position=start)
  }

calcRelativeVarCorsWithCi <- function(dnms, 
                                      withBatch=TRUE,
                                      nsim = 500) {
  if(withBatch) {
    estim <-
      rpt(n ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears + (1|familyNr) + (1|batch),
          data = dnms,
          grname=c("familyNr", "batch", "Fixed", "Residual"),
          nboot=nsim,
          npermut = 0,
          adjusted = FALSE,
          datatype = "Gaussian",
          parallel = TRUE)
  } else {
    estim <-
      rpt(n ~ fathersAgeAtConceptionInYears + mothersAgeAtConceptionInYears + (1|familyNr),
          data = dnms,
          grname=c("familyNr", "Fixed", "Residual"),
          nboot=nsim,
          npermut = 0,
          adjusted = FALSE,
          datatype = "Gaussian",
          parallel = TRUE)
  }
  res <- 
    estim %>%
    with(bind_cols(estimate=t(.$R), .$CI_emp %>% rownames_to_column())) %>% 
    dplyr::rename(factor=rowname, 
                  relVar.total=estimate, 
                  lower.total=`2.5%`, 
                  upper.total=`97.5%`)
  return(res)
}




filterAndCalculate <- function(sasani, uniqGenome) {
  filterMuts(sasani, uniqGenome) %>% 
    group_by(SampleID, 
             fathersAgeAtConceptionInYears, 
             mothersAgeAtConceptionInYears, 
             familyNr) %>% 
    count() %>% 
    calcRelativeVarCorsWithCi(withBatch = FALSE) %>%
    filter(factor == "familyNr")
}


# collect the results

familyEstimates <-
  bind_rows(
    filterAndCalculate(inova1, uniqGenome),
    filterAndCalculate(inova2, uniqGenome),
    filterAndCalculate(sasani, uniqGenome),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut %>% 
        mutate(n = unique_snv), 
      withBatch = FALSE) %>%
      filter(factor == "familyNr"),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut2 %>% 
        mutate(n = snvs_unique), 
      withBatch = FALSE)  %>%
      filter(factor == "familyNr")
  ) %>%
  bind_rows(
    filterAndCalculate(inova1, recRepGenome),
    filterAndCalculate(inova2, recRepGenome),
    filterAndCalculate(sasani, recRepGenome),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut %>% 
        mutate(n = recent_snv), 
      withBatch = FALSE) %>%
      filter(factor == "familyNr"),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut2 %>% 
        mutate(n = snvs_recent), 
      withBatch = FALSE)  %>%
      filter(factor == "familyNr")
  ) %>%
  mutate(set = rep(c("inova1", "inova2", "sasani", "autism1", "autism2"),2)) %>%
  mutate(region = c(rep("unique", 5), rep("recentReps", 5)))


ggplot(familyEstimates,
       aes(x=set, 
           y=relVar.total,
           ymin=lower.total,
           ymax=upper.total)) + 
  geom_pointrange()


