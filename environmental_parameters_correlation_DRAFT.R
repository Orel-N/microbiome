## Correlation between environmental parameters
library(rstatix);packageVersion("rstatix")
library(dplyr);packageVersion("dplyr")
library(tidyr);packageVersion("tidyr")
library(tidyverse);packageVersion("tidyverse")


## Env. parameter coorelation test 

```{r envpar-selfcorr}
#define scale function
scale_par <- function(x) scale(x, center = FALSE, scale = TRUE)[,1]

#import nv.par. data
metadata.raw <- read.csv("../Data/Consumed_nutrients_PS99_2.csv", sep = ",", dec = ".", header = TRUE)

env.par <- c("Temperature","Salinity","Chla_conc", "d.NO3.", "d.PO4.", "d.SiO3.", "NH4")

#scale all the parameters and extract surface samples
metadata.scaled.SRF <- metadata.raw %>% 
  mutate_at(env.par, scale_par)%>%
  as.data.frame()

#correlation between the env parameters
envpar_corr <- metadata.scaled.SRF %>% 
  select(env.par)%>% 
  cor_mat(method = "pearson")

#check the p-values
envpar_corr.pvalues<- envpar_corr %>% cor_get_pval()

#plot
envpar_corr %>%
  cor_reorder() %>%
  pull_lower_triangle() %>%
  cor_plot(label = TRUE)




#####################################
#Figure 6: RDA ordination of the fitted model of bacterial composition data constrained by environmental variables
####################################
#subset by fraction and remove NAs in metadata
BAC_FL.no.na <- BAC_pruned.vst %>% 
  subset_samples(
    !is.na(dNO3) & 
      !is.na(dSiO3) &
      !is.na(dPO4) &
      !is.na(ChlA)&
      Fraction =="0.22")
##remove unobserved OTU
BAC_FL.no.na <- prune_taxa(taxa_sums(BAC_FL.no.na)>0,BAC_FL.no.na)

#extract and scale the env. parameters
BAC_FL.env <- data.frame(sample_data(BAC_FL.no.na))[c("Temperature", "Salinity", "dNO3", "dSiO3", "dPO4", "ChlA")]  
BAC_FL.env <- as.data.frame(scale(BAC_FL.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_FL.otu <- t(otu_table(BAC_FL.no.na))

#Particles-associated fraction
BAC_PA.no.na <- BAC_pruned.vst %>% 
  subset_samples(
    !is.na(dNO3) & 
      !is.na(dSiO3) &
      !is.na(dPO4) &
      !is.na(ChlA)&
      Fraction =="3")
##remove unobserved OTU
BAC_PA.no.na <- prune_taxa(taxa_sums(BAC_PA.no.na)>0,BAC_PA.no.na)

#extract and scale the env. parameters
BAC_PA.env <- data.frame(sample_data(BAC_PA.no.na))[c("Temperature", "Salinity", "dNO3", "dSiO3", "dPO4", "ChlA")]  
BAC_PA.env <- as.data.frame(scale(BAC_PA.env,center = FALSE, scale = TRUE))

#extract OTU tables from Phyloseq object
BAC_PA.otu <- t(otu_table(BAC_PA.no.na))


#RDA analysis
BAC_FL.rda.all <- rda (BAC_FL.otu ~ ., data = BAC_FL.env) # model including all variables 
BAC_PA.rda.all <- rda (BAC_PA.otu ~ ., data = BAC_PA.env) # model including all variables 

#generate an RDA plot 
#FL
BAC_FL.rda.scores <- vegan::scores(BAC_FL.rda.all,display=c("sp","wa","lc","bp","cn"))
BAC_FL.rda.sites <- data.frame(BAC_FL.rda.scores$sites)
BAC_FL.rda.sites$Sample.ID <- as.character(rownames(BAC_FL.rda.sites))
sample_data(BAC_FL.no.na)$Sample.ID <- as.character(rownames(sample_data(BAC_FL.no.na)))
sample_data(BAC_FL.no.na)$Type <- as.character(sample_data(BAC_FL.no.na)$Type )
BAC_FL.rda.sites <- BAC_FL.rda.sites %>%
  left_join(sample_data(BAC_FL.no.na))

#Draw biplots
BAC_FL.rda.arrows<- BAC_FL.rda.scores$biplot*5
colnames(BAC_FL.rda.arrows)<-c("x","y")
BAC_FL.rda.arrows <- as.data.frame(BAC_FL.rda.arrows)
BAC_FL.rda.evals <- 100 * (BAC_FL.rda.all$CCA$eig / sum(BAC_FL.rda.all$CCA$eig))

#Plot 
BAC_FL.rda.plot <- ggplot() +
  geom_point(data = BAC_FL.rda.sites, aes(x = RDA1, y = RDA2, fill = Region), 
             shape =21, size = 4) +
  geom_text(data = BAC_FL.rda.sites,aes(x = RDA1, y = RDA2,label = StationName), 
            nudge_y= -0.3,size=3)+
  labs(x = sprintf("RDA1 [%s%%]", round(BAC_FL.rda.evals[1], 2)), 
       y = sprintf("RDA2 [%s%%]", round(BAC_FL.rda.evals[2], 2))) +
  scale_fill_manual(values = reg_colours) +
  #scale_x_reverse()+ 
  theme(legend.position = "none")+
  geom_segment(data=BAC_FL.rda.arrows, aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.5)+
  geom_text(data=as.data.frame(BAC_FL.rda.arrows*1.2),
            aes(x, y, label = rownames(BAC_FL.rda.arrows)),color="black",alpha=0.5)