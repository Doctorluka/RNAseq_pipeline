my_volcanoplot <- function(diff,
                           logFC_cutoff = 0.58,
                           pvalue_cutoff = 0.05,
                           label = F,
                           label.n = 7,
                           label.color = "#A58533",
                           color = c("blue", "grey","#D31638"),
                           save = F,
                           filename = NULL,
                           save.width = 6, 
                           save.height = 6){
  #test
  # diff = HYP_res
  # logFC_cutoff = 0.58
  # pvalue_cutoff = 0.05
  # label.n = 7
  # label.color = "#A58533"
  
  library(ggplot2)
  library(tidyverse)
  library(stringr)
  
  diff <- diff %>% na.omit() #去除NA值
  
  p <- ggplot(diff, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(alpha=0.5, size=3.5, 
               aes(color = regulated))+
    ylab("-log10(Pvalue)") +
    scale_color_manual(values = color) +
    geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(pvalue_cutoff),lty=4,col="black",lwd=0.8) +
    theme_bw()+
    theme(panel.grid.major=element_blank(),  #消除背景线条
          panel.grid.minor=element_blank())+
    theme(legend.justification = c(1,0),  #调整legend
          legend.position = c(1,0),
          legend.background = element_rect(fill = 'white', colour = 'black'))
  
  if (label) {
    a <- diff %>% 
      dplyr::filter(regulated != "unchanged") %>% 
      dplyr::arrange(desc(abs(logFC)))
    a1 <- c(head(a$gene_id,label.n))
    
    b <- diff %>% 
      dplyr::filter(regulated != "unchanged") %>% 
      dplyr::arrange(P.Value)
    b1 <- c(head(b$gene_id,label.n))
    
    label.text <- c(a1, b1)
    label.text <- label.text[!duplicated(label.text)]
    
    for_label <- diff %>% dplyr::filter(gene_id %in% label.text)
    
    p_label <- p +
      geom_point(size = 3.7, shape = 1, data = for_label, color = 'black')+
      ggrepel::geom_label_repel(
        aes(label = gene_id), 
        data = for_label, 
        color = label.color)
    p_label
    
    plot(p_label)
  }else{
    plot(p)
  }
  
  if (save) {
    ggsave(p_label, file = filename, 
           width = save.width, height = save.height)
  }
  
}
