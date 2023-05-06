library(ggplot2)
library(tidyverse)
library(stringr)
library(ggrepel)

# test data
# data <- res
# label.gene <- NULL
# label.gene <- c("Nfkbia","Bcr","Cdkn1a","Cdk1","Cdkn3","Foxs1",
#                 "Foxc2","Mmp14","Mmp2","Mmp23","Mmp12","Melk",
#                 "Mmp8","Mmp9","Timp1","Pdgfc","Pdgfrl","Ets1",
#                 "Angptl7","Tgfbr3","Tgfbr1","Tgfb3","Tgfbr2",
#                 "Tgfb1","Ctgf","Pcna")


volcanoPlotter <- function(data, label.gene) {
  
  data <- data[order(data$P.Value),]
  # toPlot
  toPlot <- data.frame(
    gene = data$symbol,
    p.value = -log10(data$P.Value),
    lfc = data$logFC,
    expression = data$baseMean)
  
  toPlot <- toPlot %>%
    mutate(response = ifelse(
      toPlot$lfc > 0 & toPlot$p.value > 1.30103,
      yes = "induced",
      no = ifelse(
        toPlot$lfc < 0 & toPlot$p.value > 1.30103,
        yes = "suppressed",
        no = "none"
      )
    ))
  toPlot <- toPlot[!is.na(toPlot$response),]
  
  # labeled gene
  if (length(label.gene) == 0) {
    tolabelled <- rbind(top_n(toPlot[toPlot$p.value > 1.30103 &
                                       toPlot$response == "induced",], n = 15, wt = p.value), 
                        top_n(toPlot[toPlot$p.value > 1.30103 & 
                                       toPlot$response == "suppressed",], n = 15, wt = p.value))
  }else{
    tolabelled <- toPlot[toPlot$gene %in% label.gene,]
    }
  
  # xlimit
  xlimit <- ifelse(max(toPlot$lfc) > min(toPlot$lfc) * -1,
                   max(toPlot$lfc) + 0.25,
                   (min(toPlot$lfc) * -1) + 0.25)
  
  toPlot$response <-
    factor(toPlot$response, levels = c("suppressed", "none", "induced"))

  volc <- ggplot(toPlot, aes(x = lfc, y = p.value)) +
    geom_point(
      aes(fill = response, colour = response),
      size = 2.5,
      alpha = 0.72,
      na.rm = T,
      pch = 21
    ) +
    scale_fill_brewer(palette = "Set3", direction = -1) +
    scale_colour_brewer(palette = "Set2", direction = -1) +
    theme_bw(base_size = 16) + # change theme
    theme(legend.position = c(0.85,0.85), # legend position
          legend.background = element_rect(fill = NA, colour = NA))+ # legend background
    xlab(expression(log[2]("Treatment" / "Untreated"))) + # x-axis label
    ylab(expression(-log[10]("p-value"))) + # y-axis label
    geom_vline(xintercept = c(-2, 2), colour = "darkgrey") + # Add cutoffs
    geom_hline(yintercept = 1.3, colour = "darkgrey") + # Add cutoffs
    geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
    scale_x_continuous(limits = c(-xlimit, xlimit)) +
    geom_point(
      data = tolabelled,
      color = "black",
      size = 2.5,
      shape = 21
    )+
    geom_text_repel(
      data = tolabelled,
      mapping = aes(label = gene),
      size = 3.5,
      fontface = 'italic',
      color = '#5A5A5A',
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.5, "lines")
    )
  volc
}
