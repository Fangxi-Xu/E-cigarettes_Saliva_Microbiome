library(readxl)
library(phyloseq); packageVersion("phyloseq")
library(ggplot2);packageVersion("ggplot2") 
library(ggpubr); packageVersion("ggpubr")
#import data function

import_data <- function()


# Set Themes for Figures
alpha.theme <- theme(axis.title.y=element_text(size=16, color = "black"),
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y=element_text(size=14, color="black"),
                     strip.text = element_text(size = 12,  color="black"),
                     plot.title = element_text(hjust = 0.5, size=12, face = "bold"),
                     legend.position="top",
                     legend.text=element_text(size=16, color="black"),
                     legend.title = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA),
                     )


beta.theme <- theme_classic() + 
              theme(plot.title = element_text(size = 12, color = "black"),
                    axis.text.x=element_text(size=12, color = "black"),
                    axis.title.x=element_text(size=12, color = "black"),
                    axis.title.y=element_text(size=12, color = "black"),
                    axis.text.y=element_text(size=12,  color = "black"),
                    legend.title=element_text(size=12, color = "black"), 
                    legend.text=element_text(size=12, color = "black"),
                    legend.position="bottom",
                    plot.background = element_rect(fill = "transparent",colour = NA),
              )

stacked_bar.theme <- theme_classic()+
  theme(axis.text.x=element_text(size=12, angle = 10, color = "black", face = "bold"))+
  theme(axis.text.y=element_text(size=12, color = "black", face = "bold"))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y=element_text(size=16, color = "black"))+
  theme(legend.text=element_text(size=14, color = "black"))+
  theme(legend.title=element_text(size=16, color = "black"))+
  theme(legend.position="right")


tern.theme <- theme_ggtern(base_size = 16, base_family = "")+
  theme(axis.title.y=element_text(size=16, color = "black"),
                     axis.text.x=element_blank(),
                     axis.title.x=element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.text.y=element_text(size=16, color="black"),
                     strip.text = element_text(size = 16,  color="black"),
                     plot.title = element_text(hjust = 0.5, size=16, face = "bold"),
                     legend.position="top",
                     legend.text=element_text(size=14, color="black"),
                     legend.title = element_blank(),
                     plot.background = element_rect(fill = "transparent",colour = NA),
)

#lefse.theme <- theme(axis.text.x=element_text(size=12, color = "black"))+
#  theme(axis.title.x=element_text(size=12, color = "black"))+
#  theme(axis.title.y=element_text(size=12, color = "black"))+
#  theme(axis.text.y=element_text(size=12, color = "black"))+
#  theme(legend.title=element_text(size=12, color = "black"), 
#        legend.text=element_text(size=12,  color = "black"),
#        axis.text.y=element_text(size=12, color="black"),
#        strip.text = element_blank(),
#        legend.position="bottom")