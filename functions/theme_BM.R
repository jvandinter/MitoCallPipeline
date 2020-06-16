# Theme_BM
theme_BM <- function (base_size = 20, base_family = "") 
{
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.ticks = element_line(colour = "black"),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.backgroun = element_blank(), 
          axis.line	= element_line(colour = "black", size = 0.5),
          axis.ticks.x = element_line(colour = "black", size = 0.5),
          axis.text = element_text(size = rel(0.8), colour = "black"), 
          strip.background = element_blank(),
          legend.key = element_rect(fill = "white", colour = "white"),
    )
}
