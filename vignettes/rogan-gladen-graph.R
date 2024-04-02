library(tidyverse)
library(patchwork)

## Supp 2 figure 1 ----

# not vectorised
# .pos_test_p(1, 0, 0.8, 0.9975, 100, 30, 800)
# .pos_test_p(1, 0, 0.8, 0.9975, 100)
.pos_test_p = function(test, prev, sens, spec, size, sens_n=NA, spec_n=NA, ...) {
  # This is calculating the probability N positives are observed
  # which is the probability of 0 TP * N FP + 1 TP*(N-1) FP + ... N TP * 0 FP
  sum(
    # TRUE POSITIVES:
    (if (is.na(sens_n)) {
      dbinom(test:0, prob = sens*prev, size=size) 
    } else {
      if (prev==0) {
        c(rep(0,test),1)
      } else {
        extraDistr::dbbinom(test:0, alpha = sens*prev*sens_n, beta=(1-sens*prev)*sens_n, size=size)
      }
    })
    *
    (if (is.na(spec_n)) {
      # FALSE POSITIVES:
      dbinom(0:test, prob = (1-spec)*(1-prev), size=size)
    } else {
      if (prev==1) {
        c(1,rep(0,test))
      } else {
        extraDistr::dbbinom(0:test, alpha = (1-spec)*(1-prev)*spec_n, beta=(1-(1-spec)*(1-prev))*spec_n, size=size)
      }
    })
  )
}

.gg_hide_X_axis= function () {
  ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(), 
                 axis.text.x.bottom = ggplot2::element_blank(), axis.text.x.top = ggplot2::element_blank())
}

.gg_hide_Y_axis=function () {
  ggplot2::theme(axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(), 
                 axis.text.y.left = ggplot2::element_blank(), axis.text.y.right = ggplot2::element_blank())
}


.apparent_prevalence = function(p, sens, spec) {
    p * sens + (1 - p) * (1 - spec)
}

# rogan_gladen_plot()
# This is with uncertainty in sens and spec. This is expressed as the sample size
# of the sample that informs sensitivity or specificity 
# rogan_gladen_plot(sens_n = 30, spec_n=100, examples=c(0.05,0.3))
# a UAD example
# rogan_gladen_plot(sens = 0.8, spec = 0.9975, size = 1000, examples = c(0.002,0.05), spec_n=800, xlim=0.1)
# rogan_gladen_plot(sens = 0.8, spec = 0.9975, size = 1000, examples = c(0.002,0.05), sens_n = 120, spec_n=800,xlim=0.1)
# rogan_gladen_plot(sens = 0.8, spec = 0.9975, size = 1000, examples = c(0.002,0.05), spec_n=800, xlim=0.1)
# rogan_gladen_plot(sens = 0.5, spec = 0.99, size = 1000, examples = c(0,0.01,0.05), xlim=0.1)
rogan_gladen_plot = function(
  sens = 0.7,
  spec = 0.96,
  size = 100,
  examples = c(0.01, 0.1, 0.2, 0.4),
  sens_n=NA, 
  spec_n=NA,
  xlim = 1
) {

  tmp = crossing(
    prev=seq(0,xlim,length.out=101),
    test=seq(0,floor(xlim*size),1),
  ) %>% mutate(
    sens = sens,
    spec = spec,
    size = size,
    sens_n = sens_n,
    spec_n = spec_n,
    ap = .apparent_prevalence(prev,sens,spec),
    test_n = test/size,
  ) %>% mutate(
    p_test = purrr::pmap_dbl(., .pos_test_p, .progress = TRUE)
  )
  
  marks = tibble(
    p_mark = round(examples/(xlim/100))*(xlim/100)
  ) %>% mutate(
    ap_mark = .apparent_prevalence(p_mark, sens, spec),
    test_mark = round(ap_mark*size),
    binom::binom.confint(test_mark,size,method="wilson"),
    i = row_number() %% 2 + 1
  )
  
  ylim = xlim*size # .apparent_prevalence(xlim, sens, spec)
  
  thres = (1 - spec)/(1 - sens + 1 - spec)
  
  # Colours for example curves
  p_cols = scales::brewer_pal(palette = "Dark2")(length(examples))
  # p_cols = c("blue","red")
  names(p_cols) = as.character(marks$p_mark)
  test_cols = p_cols
  names(test_cols) = as.character(marks$test_mark)
  
  marginal_x = tmp %>% inner_join(marks, by = c("test"="test_mark"))
  marginal_y = tmp %>% inner_join(marks, by = c("prev"="p_mark"))
  
  app2 = ggplot(tmp, aes(x=prev,y=ap*size))+
    geom_tile(aes(x=prev, y=test, fill = p_test),width = xlim/101*1.025, height=1.025, data = tmp, inherit.aes = FALSE, linewidth=0)+
    geom_abline(colour="grey40",slope = size)+
    geom_line()+
    scale_fill_gradient(trans="sqrt",high = "#202020" , low = "#FFFFFF00", guide="none",name="probability", oob=scales::squish, limits=c(0,0.15))+
    scale_y_continuous(sec.axis = sec_axis( trans=~./size, name="apparent prevalence (AP)", breaks=round(unique(marks$ap_mark)*1000)/1000, labels = sprintf("%1.1f%%",round(unique(marks$ap_mark)*1000)/10)), expand = c(0, 0))+
    scale_x_continuous(expand = c(0, 0), sec.axis = sec_axis( trans=~.*size, name=sprintf("disease count (N=%d)",size), breaks = round(unique(marks$p_mark*size))),labels = ~ paste0(.x*100,"%"))+
    theme(axis.title.y.left = element_blank(),axis.title.x.bottom = element_blank())+
    # .gg_hide_X_axis()+
    geom_hline(data = marks, mapping=aes(yintercept=ap_mark*size, color = as.character(p_mark)), linetype="dotted")+
    geom_vline(data = marks, mapping=aes(xintercept=p_mark, color = as.character(p_mark)), linetype="dashed")+
    geom_point(data = marks, mapping=aes(x=p_mark, y=ap_mark*size, color = as.character(p_mark)), size=1)+
    geom_point(x=thres,y=thres*size,colour="white", size=1)+
    # ggrepel::geom_label_repel(data=tibble(x=thres,y=thres*size,label=sprintf("(%1.2g,%1.2g)",thres,thres)), mapping=aes(x=x,y=y,label=label), size = 6/ggplot2:::.pt,colour="black",label.padding = unit(0.1, "lines"),label.r = unit(0, "lines"), direction = "x", label.size=0.1, inherit.aes = FALSE)+
    # geom_errorbarh(data= marks, mapping=aes(y=ap_mark*size, x=mean, xmin=lower, xmax=upper, colour = as.character(p_mark)))+
    scale_color_manual(values = p_cols, guide="none")+
    coord_cartesian(clip="on", xlim=c(0,xlim), ylim=c(0,ylim))
  
  mx = ggplot(marginal_x, aes(x=prev, y=p_test/max(p_test)*0.2, colour=as.character(test),fill=as.character(test)))+geom_area(position=position_identity(),alpha=0.1)+
    scale_color_manual(values = test_cols, name = "count", aesthetics = c("colour","fill"), guide="none")+
    geom_vline(data=marks, mapping=aes(xintercept=p_mark, colour=as.character(test_mark)), linetype = "dashed")+
    geom_rect(data=marks, mapping=aes(ymin=-0.005-(0.01*i), ymax=+0.0025-(0.01*i), xmin=lower, xmax=mean, fill = as.character(test_mark), colour = as.character(test_mark)), alpha=0.2, inherit.aes = FALSE)+
    geom_rect(data=marks, mapping=aes(ymin=-0.005-(0.01*i), ymax=+0.0025-(0.01*i), xmin=mean, xmax=upper, fill = as.character(test_mark), colour = as.character(test_mark)), alpha=0.2, inherit.aes = FALSE)+
    # ggplot2::geom_text(data=marks,mapping=ggplot2::aes(x=p_mark,y=Inf,label=p_mark,colour=as.character(test_mark)), size = 8/ggplot2:::.pt, vjust=1, hjust=1.1, angle=30, inherit.aes = FALSE)+
    scale_x_continuous(expand = c(0, 0), name="true prevalence (P)",breaks = marks$p_mark,labels = ~ paste0(.x*100,"%"))+
    scale_y_reverse(expand = c(0.005, 0.005),breaks = NULL)+.gg_hide_Y_axis()+
    theme(
      axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0))
    )+
    coord_cartesian(clip="on", xlim=c(0,xlim))
  
  
  my = ggplot(marginal_y, aes(x=test, y=p_test/max(p_test)*0.2, colour=as.character(prev),fill=as.character(prev)))+
    scale_color_manual(values = p_cols, name="prevalence",aesthetics = c("colour","fill"), guide="none")+
    geom_bar(stat="identity",position=position_identity(),colour=NA,alpha=0.1)+
    geom_step(direction="mid")+
    # ggplot2::geom_text(data=marks,mapping=ggplot2::aes(x=ap_mark*size,y=Inf,label=round(ap_mark*size),colour=as.character(p_mark)), size = 10/ggplot2:::.pt, vjust=0, hjust=1.1, angle=60, inherit.aes = FALSE)+
    geom_vline(data=marks, mapping=aes(xintercept=ap_mark*size,colour=as.character(p_mark)), linetype = "dotted")+
    ggplot2::scale_x_continuous(expand = c(0, 0), name=sprintf("test positive count (N=%d)",size), breaks=marks$ap_mark*size, labels = sprintf("%d",round(marks$ap_mark*size)))+
    coord_flip(clip="on", xlim=c(0,ylim))+
    scale_y_reverse(expand = c(0.01,0),breaks = NULL)+.gg_hide_X_axis()+
    theme(
      axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))
    )
  
  sens_label = if (is.na(sens_n)) {
    sprintf("sens: %1.4g",sens)
  } else {
    sprintf("sens: %1.4g\n[%1.4g\u2013%1.4g]",sens, qbeta(0.025, sens_n*sens, sens_n*(1-sens) ), qbeta(0.975, sens_n*sens, sens_n*(1-sens) ))
  }
  
  spec_label = if (is.na(spec_n)) {
    sprintf("spec: %1.4g",spec)
  } else {
    sprintf("spec: %1.4g\n[%1.4g\u2013%1.4g]",spec, qbeta(0.025, spec_n*spec, spec_n*(1-spec)), qbeta(0.975, spec_n*sens, spec_n*(1-spec)))
  }
  
  thres_label = sprintf("\nAP=P: %1.2g\n(n=%d)",thres,round(thres*size))
  
  info = ggplot()+annotate("text",x = 0,y=0,label=paste0(c(sens_label,spec_label, thres_label),collapse="\n"), size = 6/ggplot2:::.pt)+
    .gg_hide_X_axis()+.gg_hide_Y_axis()+theme_void()
  
  sf1 = my + app2 + info + mx + plot_layout(design =
  "ABBBB
  ABBBB
  ABBBB
  ABBBB
  CDDDD
  ",guides="collect")
  
  sf1
}
