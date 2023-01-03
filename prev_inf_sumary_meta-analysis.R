##Prev inf meta-analysis
#was created a spreadsheet (past_inf_meta) using all results from forest plots/meta-analysis

meta_pinf <- read.csv("FILE_PATH/past_inf_meta.csv")

library(ggplot2)
library(scales)

gg1 <-  ggplot(meta_pinf, aes(x = factor(severity, level = c("infection", "symptomatic", "severe")), y = variant)) +
  geom_tile(aes(fill = efficacy, height=1, width = 0.7)) + 
  geom_text(aes(label = ve_low_upp), size = 5) +
  scale_fill_gradient(low = "peachpuff2", high = "royalblue", na.value = "white", breaks=c(0, 25, 50, 75, 100), limits=c(0,100)) +
  labs(x = " ", y= " ", fill = "Previous infection protection") + 
  theme_grey(base_size=12) +
  theme(legend.text=element_text(face="bold"), 
        legend.key.size = unit(4, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(3, 'cm'), #change legend key width
        #legend.title = element_text(size=20),
        plot.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.ticks = element_blank())+
  scale_x_discrete(labels=c("Panel A: Protection against reinfection ", "Panel B: Protection against symptomatic disease", "Panel C: Protection against severe disease")) +
  theme(legend.position="bottom") +
  ggtitle(" ", subtitle = " ") +
  theme(plot.title = element_text(hjust = 0.5))    
print(gg1)

