import_sparra_expr <- function(f) {
  sink_marker <- "***************"
  f2 <- readChar(f, file.info(f)$size)
  parse(text = strsplit(gsub("\r", "", f2), sink_marker, fixed = TRUE)[[1]][2])
}




##' build_diff
##' Prepares a data frame for a ggplot object to compare differences using linear interpolation.
##' 
##' @param df data frame
##' @param xvar name of variable to consider as 'x': interpolate over evenly spaced values of this variable.
##' @return data frame using (common) interpolated x values rather than arbitrary x values
##' @examples 
##' # Only used internally
build_diff = function(df, xvar) {
  models = levels(df$Model)
  
  # Range of column 'xvar' of df
  xrange = range(df[xvar])
  
  # Equal-spaced sequence over xvar
  dfx = data.frame(seq(xrange[1], xrange[2], length = 100)[2:99])
  
  # Interpolate xvar and yvar at equal-spaced sequence for each model
  df2=data.frame()
  xnames=setdiff(colnames(df),c(xvar,"Model"))
  for(m in models) {
    d0=data.frame(Model=m,dfx[,1])
    for (yvar in xnames) d0=cbind(d0,suppressWarnings(
      approx(df[which(df$Model==m),xvar],df[which(df$Model==m),yvar],dfx[,1],rule = 2)$y
    ))
    df2=rbind(df2,d0)
  }
  colnames(df2)=c("Model",xvar,xnames)
  df2$Model=factor(df2$Model,levels=models)
  
  return(df2)
}



### Dumped from fairness package



##' roc_2panel
##' Draws a ROC curve (with legend) with a second panel underneath showing sensitivity difference.
##' 
##' @param rocs list of sparraROC objects
##' @param labels labels to use in legend; default to names of rocs.
##' @param col line colours
##' @param xy_col line colour for x-y line, defaults to red
##' @param highlight if non-null, add a point at this cutoff
##' @param yrange_lower y range for lower plot. If NULL, generates automatically
##' @param legend_title title for legend, defaults to nothing
##' @return Invisibly returns plot as ggplot object
##' @export
##' @examples 
##' # See vignette
roc_2panel_gg=function(rocs,labels=names(rocs),col=1:length(rocs),
                    xy_col="black",
                    highlight=NULL,yrange_lower=NULL,legend_title="") {
  
  # Fix colours
  if (length(col)>length(rocs)) col=col[1:length(rocs)]
  
  # Compose data frame for ggplot
  df=data.frame()
  for (i in 1:length(rocs)) 
    df=rbind(df,
             data.frame(Model=labels[i],
                        sens=rocs[[i]]$sens[1,],
                        spec=rocs[[i]]$spec[1,]))
  df$Model=factor(df$Model,levels=unique(df$Model))
  
  # Assemble top figure
  df1=df
  df1$spec=1-df1$spec
  p1 = ggplot(df1) +
    geom_path(aes(x = .data$spec, y = .data$sens, col = .data$Model), linewidth = 0.4) +
    xlim(0, 1) + ylim(0, 1) +
    xlab("") + ylab("Sensitivity") +
    theme_minimal(base_size = 8) + theme(legend.justification = c(1,0),
                                         legend.position = c(1,0),
                                         legend.spacing = unit(0, "npc"),
                                         legend.margin = margin(),
                                         legend.background = 
                                           element_rect(fill = "white", linewidth=0, colour = "white"))
  
  
  
  
  # Set up data frame for lower panel
  df2=build_diff(df,"spec")
  df2$spec=1-df2$spec
  bx=df2$sens[which(df2$Model==labels[1])]
  for (i in 1:length(labels)) {
    w=which(df2$Model==labels[i])
    df2$sens[w]= df2$sens[w]-bx
  }
  colnames(df2)[3]="dsens"
  
  p2 = ggplot(df2) + 
    geom_path(aes(x = .data$spec, y = .data$dsens, col = .data$Model), linewidth = 0.4) +
    xlim(0, 1) + 
    xlab("1 - Specificity") + ylab(expression(paste(Delta," Sensitivity"))) +
    theme_minimal(base_size = 8) + theme(legend.position = "none")
  
  
  ## Modify styles
  # Line colours
  p1=p1 + scale_color_manual(values=col,name=legend_title)
  p2=p2 + scale_color_manual(values=col)
  # A-B line
  p1=p1 + geom_abline(slope=1,color=xy_col,linetype=2)
  # Scale for lower plot
  if (!is.null(yrange_lower)) p2=p2 + ylim(yrange_lower[1],yrange_lower[2])
  
  # Highlights
  if (!is.null(highlight)) {
    xy=data.frame(t(data.frame(lapply(rocs,function(x) {
      w0=which.min(abs(x$cutoffs-highlight))
      return(c(1-x$spec[w0],x$sens[w0]))
    }))))
    colnames(xy)=c("x","y")
    p1 = p1 + geom_point(aes(.data$x,.data$y),data=xy,color=col) 
    
  }
  
  # Combine figures
  p = p1 / p2 + plot_layout(heights = c(3,1)) 
  
  suppressWarnings(print(p))
  invisible(p)
}








##' prc_2panel
##' Draws a PRC curve (with legend) with a second panel underneath showing precision difference.
##' 
##' @param prcs list of sparraPRC objects.
##' @param labels labels to use in legend
##' @param col line colours
##' @param highlight if non-null, draw a point at a particular cutoff
##' @param yrange_lower y range for lower plot. If NULL, generates automatically
##' @param legend_title title for legend, defaults to nothing
##' @return Silently return ggplot object
##' @export
##' @examples 
##' # See vignette
prc_2panel_gg=function(prcs,labels=names(prcs),col=1:length(prcs),
                    highlight=NULL,yrange_lower=NULL,legend_title="") {
  
  # Fix colours
  if (length(col)>length(prcs)) col=col[1:length(prcs)]
  
  # Compose data frame for ggplot
  df=data.frame()
  for (i in 1:length(prcs)) 
    df=rbind(df,
             data.frame(Model=labels[i],
                        sens=prcs[[i]]$sens[1,],
                        ppv=prcs[[i]]$ppv[1,]))
  df$Model=factor(df$Model,levels=unique(df$Model))
  
  # Assemble top figure
  p1 = ggplot(df) +
    geom_path(aes(x = .data$sens, y = .data$ppv, col = .data$Model), linewidth = 0.4) +
    xlim(0, 1) + ylim(0, 1) +
    xlab("") + ylab("Precision") +
    theme_minimal(base_size = 8) + theme(legend.justification = c(1,0),
                                         legend.position = c(1,0),
                                         legend.spacing = unit(0, "npc"),
                                         legend.margin = margin(),
                                         legend.background = 
                                           element_rect(fill = "white", linewidth=0, colour = "white"))
  
  
  
  
  
  # Set up data frame for lower panel
  df2=build_diff(df,"sens")
  bx=df2$ppv[which(df2$Model==labels[1])]
  for (i in 1:length(labels)) {
    w=which(df2$Model==labels[i])
    df2$ppv[w]= df2$ppv[w]-bx
  }
  colnames(df2)[3]="dppv"
  
  p2 = ggplot(df2) + 
    geom_path(aes(x = .data$sens, y = .data$dppv, col = .data$Model), linewidth = 0.4) +
    xlim(0, 1) + 
    xlab("Recall") + ylab(expression(paste(Delta," Precision"))) +
    theme_minimal(base_size = 8) + theme(legend.position = "none")
  
  
  ## Modify styles
  # Line colours
  p1=p1 + scale_color_manual(values=col,name=legend_title)
  p2=p2 + scale_color_manual(values=col)
  # Scale for lower plot
  if (!is.null(yrange_lower)) p2=p2 + ylim(yrange_lower[1],yrange_lower[2])
  
  # Highlights
  if (!is.null(highlight)) {
    xy=data.frame(t(data.frame(lapply(prcs,function(x) {
      w0=which.min(abs(x$cutoffs-highlight))
      return(c(x$sens[w0],x$ppv[w0]))
    }))))
    colnames(xy)=c("x","y")
    p1 = p1 + geom_point(aes(.data$x,.data$y),data=xy,color=col) 
    
  }
  
  # Combine figures
  p = p1 / p2 + plot_layout(heights = c(3,1)) 
  
  suppressWarnings(print(p))
  invisible(p)
  
}



##' cal_2panel
##' Draws calibration curves (with legend) with a second panel underneath showing predicted differences.
##' 
##' @param cals list of calibration objects, output from getcal(). 
##' @param labels labels to use in legend
##' @param col line colours
##' @param xy_col line colour for x-y line, defaults to phs-magenta
##' @param ci_col colours to draw confidence intervals on lower panel; NA to not draw. 
##' @param highlight if non-null, highlight a particular value
##' @param yrange_lower y range for lower plot. If NULL, generates automatically
##' @param legend_title title for legend, defaults to nothing
##' @return Silently return ggplot object
##' @export
##' @examples 
##' # See vignette
cal_2panel_gg=function(cals,labels,col=1:length(cals),
                    xy_col="black",
                    ci_col=col,highlight=NULL,yrange_lower=NULL,
                    legend_title="") {
  
  # Fix colours
  if (length(col)>length(cals)) col=col[1:length(cals)]
  
  
  # Compose data frame for ggplot
  df=data.frame()
  for (i in 1:length(cals)) 
    df=rbind(df,
             data.frame(Model=labels[i],
                        obs=cals[[i]]$x,
                        exp=cals[[i]]$y,
                        du=cals[[i]]$upper,
                        dl=cals[[i]]$lower))
  df$Model=factor(df$Model,levels=unique(df$Model))
  
  # Assemble top figure
  if (!is.null(ci_col)) ccol=ci_col[as.factor(df$Model)]
  p1 = ggplot(df)
  if (!is.null(ci_col)) p1=p1 + geom_ribbon(aes(x = .data$obs, ymin = .data$dl, 
                                                ymax = .data$du,fill=.data$Model),alpha = 0.25)
  p1=p1 +
    geom_path(aes(x = .data$obs, y = .data$exp, col = .data$Model), linewidth = 0.4) +
    xlim(0, 1) + ylim(0, 1) +
    xlab("") + ylab("Observed") +
    theme_minimal(base_size = 8) + theme(legend.justification = c(1,0),
                                         legend.position = c(1,0),
                                         legend.spacing = unit(0, "npc"),
                                         legend.margin = margin(),
                                         legend.background = 
                                           element_rect(fill = "white", linewidth=0, colour = "white"))
  
  
  
  
  
  # Set up data frame for lower panel
  df2=build_diff(df,"obs")
  for (i in 1:length(labels)) {
    w=which(df2$Model==labels[i])
    df2$exp[w]= df2$exp[w]-df2$obs[w]
    df2$du[w]=df2$exp[w]-df2$obs[w] + (df2$du[w]-df2$exp[w])
    df2$dl[w]=df2$exp[w]-df2$obs[w] + (df2$dl[w]-df2$exp[w])
  }
  colnames(df2)[3]="dexp"
  
  p2 = ggplot(df2)
  if (!is.null(ci_col)) p2=p2 + geom_ribbon(aes(x = .data$obs, ymin = .data$dl, 
                                                ymax = .data$du,fill=.data$Model),alpha = 0.25)
  p2 = p2 +  geom_path(aes(x = .data$obs, y = .data$dexp,col=.data$Model), linewidth = 0.4) +
    xlim(0, 1) + 
    xlab("Expected") + ylab(expression(paste(Delta," from calibrated"))) +
    theme_minimal(base_size = 8) + theme(legend.position = "none")
  
  
  ## Modify styles
  # Line colours
  cval=c("r1"="red","r2"="black")
  p1=p1 + scale_color_manual(values=col,name=legend_title) 
  p2=p2 + scale_color_manual(values=col)
  
  if (!is.null(ci_col)) {
    p1=p1 + scale_fill_manual(values=ci_col,name=legend_title) 
    p2=p2 + scale_fill_manual(values=ci_col)
  }
  
  # A-B line
  p1=p1 + geom_abline(slope=1,color=xy_col,linetype=2)
  # Scale for lower plot
  if (!is.null(yrange_lower)) p2=p2 + ylim(yrange_lower[1],yrange_lower[2])
  
  # Highlights
  if (!is.null(highlight)) {
    xy=data.frame(t(data.frame(lapply(cals,function(x) {
      w0=which.min(abs(x$x-highlight))
      return(c(x$x[w0],x$y[w0]))
    }))))
    colnames(xy)=c("x","y")
    p1 = p1 + geom_point(aes(.data$x,.data$y),data=xy,color=col) 
    
  }
  
  # Combine figures
  p = p1 / p2 + plot_layout(heights = c(3,1)) 
  
  suppressWarnings(print(p))
  invisible(p)
}  

