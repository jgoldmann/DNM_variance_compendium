
library(tidyverse)
library(here)
library(ggpubr)
library(cowplot)

allEstimates <- read_delim(here("dnmvar-compendium/analyses/varianceComponentEstimations/allEstimates.tsv"),
                           delim = "\t")

familyEstimates <- read_delim(here("dnmvar-compendium/analyses/varianceComponentEstimations/familyEstimates.tsv"),
                              delim = "\t")


ggplot(familyEstimates %>%
         mutate(set = fct_inorder(set)),
       aes(x=set, 
           y=relVar.total,
           ymin=lower.total,
           ymax=upper.total)) +  
  geom_errorbar(aes(size=no.childrenPerFamily), col = "grey") +
  geom_point(aes(size=no.childrenPerFamily), col="#008000", shape=18) + 
  geom_hline(yintercept = 0) + 
  labs(x=NULL, y="Family variance component") + 
  ggpubr::theme_pubr() + 
  geom_text(y=0.9, 
            cex = 2.8,
            aes(label=paste0("no. families: ",no.familys, "\n",
                             "no. children per family ", round(no.childrenPerFamily, 2))))  +
  ylim(-0.06,1)



ggdotchart(familyEstimates %>% mutate(set = fct_inorder(set)),
           x="set",
           y="relVar.total",
           rotate = TRUE,
           ggtheme = theme_pubr()) + 
  theme_cleveland() + 
  geom_errorbar(aes(ymin=lower.total,
                    ymax=upper.total))


ggplot(familyEstimates,
       aes(x=set,
           y=relVar.total,
           ymin=lower.total,
           ymax=upper.total)) +
  geom_linerange() +
  geom_point(aes(size=no.childrenPerFamily), shape=18) + 
  coord_flip() + 
  theme_pubclean() + 
  geom_hline(yintercept = 0) + 
  geom_text(aes(label=round(relVar.total,3)), nudge_x = 0.2) +
  scale_size(range = c(2,8), guide = FALSE) +
  labs(x=NULL, y="Family variance component")


ggplot(familyEstimates,
       aes(x=set,
           y=relVar.total,
           ymin=lower.total,
           ymax=upper.total)) +
  geom_linerange() +
  geom_point(aes(size=no.childrenPerFamily), shape=18) + 
  coord_flip() + 
  theme_pubclean() + 
  geom_hline(yintercept = 0) + 
  geom_hline(yintercept = 0.055, lty=2, col="green3") +
  geom_label(aes(label=round(relVar.total,3)), nudge_x = 0.2) +
  scale_size(range = c(2,8), guide = FALSE) +
  labs(x=NULL, y="Family variance component") + 
  geom_text(y=0.35, 
            cex = 2.8,
            aes(label=paste0("no. multi-offspring families: ",no.familys, "\n",
                             "no. children per family ", round(no.childrenPerFamily, 2))))


a <-
  ggplot(familyEstimates,
         aes(x=set,
             y=relVar.total,
             ymin=lower.total,
             ymax=upper.total)) +
  geom_linerange() +
  geom_point(aes(size=no.childrenPerFamily), shape=18) + 
  coord_flip() + 
  theme_pubclean() + 
  geom_hline(yintercept = 0) + 
  geom_text(aes(label=round(relVar.total,3)), nudge_x = 0.2) +
  scale_size(range = c(2,8), guide = FALSE) +
  labs(x=NULL, y="Family variance component")


b <- 
  familyEstimates %>% 
  transmute(includesBatch = as.character(includesBatch), 
            no.familys, 
            no.childrenPerFamily = round(no.childrenPerFamily, 2), 
            no.muts.pp=round(no.muts.pp,2)) %>%
  ggtexttable(rows = NULL, 
              theme = ttheme("blank"))

plot_grid(a,b)


ggplot(familyEstimates,
       aes(x=set,
           y=relVar.total,
           ymin=lower.total,
           ymax=upper.total)) +
  geom_col(aes(alpha=no.childrenPerFamily)) +
  geom_errorbar() + 
  theme_pubclean() + 
  geom_text(aes(label=round(relVar.total,3)), nudge_x = 0.2, nudge_y = 0.01) +
  scale_size(range = c(2,8), guide = FALSE) +
  labs(x=NULL, y="Family variance component") +
  geom_hline(yintercept = 0.055, lty=2, col="green3") +
  geom_text(y=-0.03, 
            cex = 2.8,
            aes(label=paste0("multi-offspring families: ",no.familys, "\n",
                             "children per family ", round(no.childrenPerFamily, 2)))) + 
  ylim(-0.05,0.4)

