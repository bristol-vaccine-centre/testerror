
# set up theme
ggplot2::theme_set(
  ggplot2::theme_bw(base_size = 8,base_line_size = 0.25,base_rect_size = 0.25)+
       ggplot2::theme(
         axis.text.x.bottom = ggplot2::element_text(angle = 30,vjust=1,hjust=1),
         legend.position = "bottom"
       )
)

default_table = function(df) {
  return(df %>% huxtable::huxtable() %>% 
    huxtable::theme_article())
}
