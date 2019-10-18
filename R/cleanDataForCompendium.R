 

library(tidyverse)
library(here)

library(rtracklayer)

library(inovaDNMs)
library(inovaDNMs102)
library(dnmVarianceComponents)


# load datasets as in result overview -------

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
  dplyr::rename(familyNr = Family.Study.ID) %>%
  mutate(familyNr = if_else(familyNr %in% c("101-112", "101-689"), "101-112/101-689", familyNr)) %>% # properly label sibling pair (see Wendy's mail 2nd oct 2016)
  mutate(PROBANDID=paste0(SampleID, "-03")) %>%
  left_join(inovaDNMs::seqReactions[,c("software", "name")] %>%
              mutate(batch = substr(software, 1,6)),
            by=c(Proband.Study.ID="name"))

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
  mutate(end = NULL) %>% 
  left_join(inovaDNMs102::siblings %>% 
              filter(!SiblingPairId == "Fam12" & !SiblingPairId == "Sib14") %>% #clean incomplete families
              dplyr::select(-PROBANDID) %>% #drop double entries
              distinct(), 
            by=c("familyNr"="FAMILYID")) %>%
  mutate(familyNr = if_else(mode=="siblings" & !is.na(mode), SiblingPairId, familyNr)) %>% 
  left_join(inovaDNMs102::qcStats %>%
              dplyr::select(pipeline, STUDY_ID) %>%
              mutate(STUDY_ID=toupper(STUDY_ID), batch=pipeline),
            by=c(SampleID="STUDY_ID"))

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
    Variant = alt, 
    Position = start,
    fathersAgeAtConceptionInYears = paternal_age_at_birth,
    mothersAgeAtConceptionInYears = maternal_age_at_birth,
    familyNr = new_family_id
  ) %>%
  filter(mut != "indel")

autism1 <- dnmVarianceComponents::dnmDataAut
autism2 <- dnmVarianceComponents::dnmDataAut2


# get genomic data for filtering mutations ------

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



# filtering function ---------
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



# un-clutter and write out sample datasets ----------

samples_c1 <- 
  inova1 %>%
  filterMuts(uniqGenome) %>% 
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr,
           batch) %>% 
  count() %>%
  ungroup() %>%
  dplyr::rename(no.DNMs = n)
write_delim(samples_c1,
            here("compendium/data/samples_c1.tsv"),
            delim = "\t")

samples_c2 <- 
inova2 %>%
  filterMuts(uniqGenome) %>% 
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr,
           batch) %>% 
  count() %>%
  ungroup() %>%
  dplyr::rename(no.DNMs = n)
write_delim(samples_c2,
            here("compendium/data/samples_c2.tsv"),
            delim = "\t")


samples_c3 <- 
  sasani %>% 
  dplyr::select(-end) %>%
  filterMuts(uniqGenome) %>% 
  group_by(SampleID, 
           fathersAgeAtConceptionInYears, 
           mothersAgeAtConceptionInYears, 
           familyNr) %>% 
  count() %>%
  ungroup() %>%
  mutate(batch = NA_character_) %>%
  dplyr::rename(no.DNMs = n)
write_delim(samples_c3,
            here("compendium/data/samples_c3.tsv"),
            delim = "\t")

samples_c4 <- 
  autism1 %>%
  transmute(
    SampleID = id,
    fathersAgeAtConceptionInYears,
    mothersAgeAtConceptionInYears,
    familyNr,
    n = unique_snv,
    batch = as.character(batch)) %>%
  dplyr::rename(no.DNMs = n)
write_delim(samples_c4,
            here("compendium/data/samples_c4.tsv"),
            delim = "\t")


samples_c5 <- 
  autism2 %>%
  transmute(
    SampleID = sample,
    fathersAgeAtConceptionInYears,
    mothersAgeAtConceptionInYears,
    familyNr,
    n = snvs_unique,
    batch = as.character(batch)) %>%
  dplyr::rename(no.DNMs = n)
write_delim(samples_c5,
            here("compendium/data/samples_c5.tsv"),
            delim = "\t")



# un-clutter and write out mutation datasets ----------

mutations_c1 <- 
  inova1 %>%
  transmute(
    Chr,
    Position,
    SampleID,
    Reference,
    Variant, 
    Phase = parent,
    Substitution = if_else(is.CpG, 
                           paste(substitution, "CpG"), 
                           as.character(substitution))
  ) %>% 
  filterMuts(uniqGenome)
write_delim(mutations_c1,
            here("compendium/data/mutations_c1.tsv"),
            delim = "\t")

mutations_c2 <- 
  inova2 %>%
  transmute(
    Chr,
    Position,
    SampleID,
    Reference = ref,
    Variant, 
    Phase = parent,
    Substitution = if_else(substitution == "C - T" & substr(surrounding,3,3) == "G", 
                           paste(substitution, "CpG"), 
                           as.character(substitution))
  ) %>% 
  filterMuts(uniqGenome)
write_delim(mutations_c2,
            here("compendium/data/mutations_c2.tsv"),
            delim = "\t")

mutations_c3 <- 
  sasani %>%
  transmute(
    Chr = paste0("chr", chrom),
    Position,
    SampleID,
    Reference = ref,
    Variant, 
    Phase = phase,
    Substitution = mut
  ) %>% 
  filterMuts(uniqGenome)
write_delim(mutations_c3,
            here("compendium/data/mutations_c3.tsv"),
            delim = "\t")
