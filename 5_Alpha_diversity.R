##################################################
#Alpha diversity plots and statistic: Dataset 2020
##################################################

#Load packages
library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)
library(iNEXT)

#Set theme
theme_set(theme_bw())

#Import collor palette
source("./scripts/Color_palettes.R")


############
#Import data
############

ps20 <- readRDS("./data/phyloseqRaw20.RDS")
psTaxF20 <- readRDS("./data/phyloseqTaxFiltered20.RDS")
psPrevF20 <- readRDS("./data/phyloseqPrevFiltered20.RDS")


#################
#Explore dataset
#################

#Number of reads, after QC
ps20
sum(sample_sums(ps20))
mean(sample_sums(ps20))
min(sample_sums(ps20))
max(sample_sums(ps20))
sample_sums(ps20)


################################
#Supplementary: Rarefaction plot
################################

ps20.iNEXT <- iNEXT(as.data.frame(otu_table(ps20)), q=0, datatype="abundance", knots = 40)

#Plot all together
ggiNEXT(ps20.iNEXT, type = 1)
ggiNEXT(ps20.iNEXT, type = 2)

#Plot separately Depth~Location
ps20.iNEXT.rare <-fortify(ps20.iNEXT, type=1)
ps20.iNEXT.meta <- as(sample_data(ps20), "data.frame")
ps20.iNEXT.meta$site <- rownames(ps20.iNEXT.meta)

ps20.iNEXT.rare$Location <- ps20.iNEXT.meta$Location[match(ps20.iNEXT.rare$site, ps20.iNEXT.meta$site)]
ps20.iNEXT.rare$Depth <- ps20.iNEXT.meta$Depth[match(ps20.iNEXT.rare$site, ps20.iNEXT.meta$site)]
ps20.iNEXT.rare$Season <- ps20.iNEXT.meta$Season[match(ps20.iNEXT.rare$site, ps20.iNEXT.meta$site)] 

ps20.iNEXT.rare.point <- ps20.iNEXT.rare[which(ps20.iNEXT.rare$method == "observed"),]
ps20.iNEXT.rare.line <- ps20.iNEXT.rare[which(ps20.iNEXT.rare$method != "observed"),]
ps20.iNEXT.rare.line$method <- factor (ps20.iNEXT.rare.line$method,
                                      c("interpolated", "extrapolated"),
                                      c("interpolation", "extrapolation"))
iNEXT.rare.line <- data.frame()
iNEXT.rare.point<- data.frame()

iNEXT.rare.line <- rbind(iNEXT.rare.line, ps20.iNEXT.rare.line)
iNEXT.rare.point <- rbind(iNEXT.rare.point, ps20.iNEXT.rare.point)

iNEXT.rare.line$Season <- factor(iNEXT.rare.line$Season, levels = c("autumn", "winter", "spring", "summer"))
iNEXT.rare.point$Season <- factor(iNEXT.rare.point$Season, levels = c("autumn", "winter", "spring", "summer"))
iNEXT.rare.line$Depth <- factor(iNEXT.rare.line$Depth, levels = c("surface", "bottom"))
iNEXT.rare.point$Depth <- factor(iNEXT.rare.point$Depth, levels = c("surface", "bottom"))
iNEXT.rare.line$Location <- factor(iNEXT.rare.line$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine", "OS-Marine", "SM-Outfall"))
iNEXT.rare.point$Location <- factor(iNEXT.rare.point$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine", "OS-Marine", "SM-Outfall"))

rare.p <- ggplot(iNEXT.rare.line, aes(x=x, y=y))+
  geom_line(aes(linetype = method, colour = Season), lwd = 1, data= iNEXT.rare.line)+
  geom_point(size =3, colour = "black", data= iNEXT.rare.point)+
  labs(x = "Sample size", y = "Species richness")+
  xlim(0,3e5)+
  theme_classic(base_size = 12)+theme(legend.position="top", axis.text.x = element_text(angle = 50, hjust = 1))+
  facet_grid(Depth~Location, scales = "fixed",as.table = TRUE)+
  geom_hline(aes(yintercept=-Inf)) + 
  geom_vline(aes(xintercept=-Inf)) +
  scale_color_manual(values = season.col)+
  coord_cartesian(clip="off")
rare.p

ggsave("./output_graphs/16S_20_rarefactions.pdf", rare.p)


############################
#Alpha diversity calculation
############################

#Calculate alpha diversity
ps20_alpha <- estimate_richness(ps20, measures = c("Observed", "Chao1","Shannon", "InvSimpson"))
ps20_alpha$Evenness <- ps20_alpha$Shannon/log(ps20_alpha$Observed)
ps20_alpha$Row.names<-rownames(ps20_alpha)
alpha_df_20 <- data.frame(sample_data(ps20), ps20_alpha)

#Edit table
alpha_df_20$Season <- factor(alpha_df_20$Season, levels = c("winter", "spring", "summer", "autumn"))
alpha_df_20$Location <- factor(alpha_df_20$Location, levels = c("R-Estuary-1", "R-Estuary-2", "NS-Marine","OS-Marine", "SM-Outfall"))
alpha_df_20$Depth <- factor(alpha_df_20$Depth, levels = c("bottom", "surface"))


#####################################
#Supplementary: Alpha diversity plots
#####################################

#Plot Chao1
Chao.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=Chao1, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Chao1, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  scale_colour_manual(values = location.col) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Chao.20.p

#Plot Shannon d.i.
Shannon.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=Shannon, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Shannon, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  scale_colour_manual(values = location.col) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Shannon.20.p

#Plot Evenness
Evenness.20.p <- ggplot(data = alpha_df_20) +
  geom_boxplot(mapping=aes(x=Location, y=Evenness, color=Location)) +
  geom_point(mapping=aes(x=Location, y=Evenness, color=Location, shape = Season), size = 3, stat = 'identity') +
  facet_wrap(~Depth) +
  scale_colour_manual(values = location.col) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1))
Evenness.20.p

alpha.20 <- ggarrange(Chao.20.p, Shannon.20.p, Evenness.20.p, nrow=1, common.legend = TRUE)
alpha.20

ggsave("./output_graphs/AlphaDiversity.pdf", alpha.20)


##########################
#Alpha diversity statistic
##########################

#Shapiro-Wilk's test
alpha_df_20 %>% shapiro_test(Chao1, Shannon, Evenness)

#ANOVA
alpha_df_20 %>% anova_test(Chao1 ~ Season + Location + Depth)
alpha_df_20 %>% tukey_hsd(Chao1 ~ Season)

alpha_df_20 %>% anova_test(Shannon ~ Season + Location + Depth)
alpha_df_20 %>% tukey_hsd(Shannon ~ Season)

alpha_df_20 %>% anova_test(Evenness ~ Season + Location + Depth)
alpha_df_20 %>% tukey_hsd(Evenness ~ Season + Location)

#ANOVA group by depth
alpha_df_20 %>% group_by(Depth) %>% anova_test(Chao1 ~ Season+Location)
alpha_df_20 %>% group_by(Depth) %>% tukey_hsd(Chao1 ~ Season+Location)

alpha_df_20 %>% group_by(Depth) %>% anova_test(Shannon ~ Season + Location)
alpha_df_20 %>% group_by(Depth) %>% tukey_hsd(Shannon ~ Season)

alpha_df_20 %>% group_by(Depth) %>% anova_test(Evenness ~ Season + Location)
alpha_df_20 %>% group_by(Depth) %>% tukey_hsd(Evenness ~ Season + Location)



# CLEAN UP #################################################

# Clear plots
dev.off()

# Clear environment
rm(list = ls()) 

# Clear console
cat("\014")
