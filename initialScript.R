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
seqlevelsStyle(recRepGenome) <- "UCSC"
ancientRepGenome <- setdiff(gaps(uniqGenome), recRepGenome)
allGenome <- SeqinfoForUCSCGenome("hg19") %>% as("GRanges")


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
  bind_rows(
    filterAndCalculate(inova1, ancientRepGenome),
    filterAndCalculate(inova2, ancientRepGenome),
    filterAndCalculate(sasani, ancientRepGenome),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut %>% 
        mutate(n = ancient_snv), 
      withBatch = FALSE) %>%
      filter(factor == "familyNr"),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut2 %>% 
        mutate(n = snvs_ancient), 
      withBatch = FALSE)  %>%
      filter(factor == "familyNr")
  ) %>%
  bind_rows(
    filterAndCalculate(inova1, allGenome),
    filterAndCalculate(inova2, allGenome),
    filterAndCalculate(sasani, allGenome),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut %>% 
        mutate(n = SNVs), 
      withBatch = FALSE) %>%
      filter(factor == "familyNr"),
    calcRelativeVarCorsWithCi(
      dnmVarianceComponents::dnmDataAut2 %>% 
        mutate(n = snvs_all), 
      withBatch = FALSE)  %>%
      filter(factor == "familyNr")
  ) %>%
  mutate(set = rep(c("inova1", "inova2", "sasani", "autism1", "autism2"), 4)) %>%
  mutate(region = c(rep("unique", 5), rep("recentReps", 5), rep("ancientReps", 5), rep("allGenome", 5)))


ggplot(familyEstimates,
       aes(x=set, 
           y=relVar.total,
           ymin=lower.total,
           ymax=upper.total)) + 
  geom_pointrange() + 
  facet_wrap(~region, ncol=1) +
  theme_bw()

ggsave(here("estimates.png"))



#some augmentative statistics as proxy for reliability
familyStats <- function(inova1, allGenome) {
  familys <- 
    filterMuts(inova1, allGenome) %>% 
    group_by(familyNr) %>% 
    summarise(no.offspring = length(unique(SampleID))) %>% 
    filter(no.offspring>=2) 
  no.muts <- 
    filterMuts(inova1, allGenome) %>% 
    filter(familyNr %in% familys$familyNr) %>% 
    nrow()
  no.familys <- 
    familys %>% 
    nrow()
  no.offspring <- 
    familys %>%
    with(sum(no.offspring))
  tibble(no.familys=no.familys, no.muts.pp=no.muts/no.offspring)
}
familyStatsFrTbl <- function(autism1, tblHdr = "unique_snv") {
  no.offspring <- sum(autism1[, tblHdr]>0)
  no.familys <- autism1$familyNr %>% unique() %>% length()
  no.muts <- sum(autism1[, tblHdr])
  tibble(no.familys=no.familys, no.muts.pp=no.muts/no.offspring)
}


familyDf <- 
  bind_rows(
    familyStats(inova1, uniqGenome),
    familyStats(inova2, uniqGenome),
    familyStats(sasani, uniqGenome),
    familyStatsFrTbl(autism1, "unique_snv"),
    familyStatsFrTbl(autism2, "snvs_unique"),
    familyStats(inova1, recRepGenome),
    familyStats(inova2, recRepGenome),
    familyStats(sasani, recRepGenome),
    familyStatsFrTbl(autism1, "recent_snv"),
    familyStatsFrTbl(autism2, "snvs_recent"),
    familyStats(inova1, ancientRepGenome),
    familyStats(inova2, ancientRepGenome),
    familyStats(sasani, ancientRepGenome),
    familyStatsFrTbl(autism1, "ancient_snv"),
    familyStatsFrTbl(autism2, "snvs_ancient"),
    familyStats(inova1, allGenome),
    familyStats(inova2, allGenome),
    familyStats(sasani, allGenome),
    familyStatsFrTbl(autism1, "SNVs"),
    familyStatsFrTbl(autism2, "snvs_all")
  )

familyEstimates <- 
  bind_cols(
    familyEstimates,
    familyDf
  )

write_delim(familyEstimates,
            here("estimates.tsv"),
            delim = "\t")
