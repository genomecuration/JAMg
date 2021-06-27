
cq_dataframe <- read.csv(cq_input, header = T, sep = "\t")

not_na <- which(cq_dataframe$Median != "N/A")
cq_dataframe <- cq_dataframe[not_na,]

cq_dataframe$lnMedian <- log(cq_dataframe$Median)

description_categories <- as.factor(sort(unique(cq_dataframe$CONTIG_DESCRIPTION)))


size_categories <- as.integer(quantile(cq_dataframe$CONTIG_SIZE, probs = c(0.25, 0.5, 0.75, 0.99))) 

# find which contig sizes in each quartile
for (row in 1:nrow(cq_dataframe)){
  size <- cq_dataframe[row, "CONTIG_SIZE"]
  
  if (size < size_categories[1]){
    cq_dataframe[row, "size"]  <- paste(paste(min(cq_dataframe$CONTIG_SIZE), size_categories[1], sep ="-"),"bp")
  }else if (size < size_categories[2]){
    cq_dataframe[row, "size"] <- paste(paste(size_categories[1], size_categories[2], sep ="-"),"bp")
  }else if (size < size_categories[3]){
    cq_dataframe[row, "size"] <- paste(paste(size_categories[2], size_categories[3], sep ="-"),"bp")
  }else if (size < size_categories[4]){
    cq_dataframe[row, "size"] <- paste(paste(size_categories[3], size_categories[4], sep ="-"),"bp")
  }else {
    cq_dataframe[row, "size"] <- paste(paste(size_categories[4], max(cq_dataframe$CONTIG_SIZE), sep ="-"),"bp")
  }
}


#human friendly labels for the different sizes
xaxis_vector <- c(paste(paste(min(cq_dataframe$CONTIG_SIZE), size_categories[1], sep ="-"),"bp"), 
                  paste(paste(size_categories[1], size_categories[2], sep ="-"),"bp"), 
                  paste(paste(size_categories[2], size_categories[3], sep ="-"),"bp"), 
                  paste(paste(size_categories[3], size_categories[4], sep ="-"),"bp"), 
                  paste(paste(size_categories[4], max(cq_dataframe$CONTIG_SIZE), sep ="-"),"bp"))

# first do the figures with different size categories

##violin plot with Median values in y
cq_violin_median <- cq_dataframe %>%
  plot_ly(
    type = "violin"
  )

cq_violin_median <- cq_violin_median %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[1]],
    y = ~ Median[cq_dataframe$size == xaxis_vector[1]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[1]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    # color = I("yellow"),
    name = xaxis_vector[1]
  )

cq_violin_median <- cq_violin_median %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[2]],
    y = ~ Median[cq_dataframe$size == xaxis_vector[2]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[2]],
    points = "all",
    pointpos = 0,
    jitter = 1,
   # color = I("yellow"),
    name = xaxis_vector[2]
  )
cq_violin_median <- cq_violin_median %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[3]],
    y = ~ Median[cq_dataframe$size == xaxis_vector[3]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[3]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    #color = I("green"),
    name = xaxis_vector[3]
  )
cq_violin_median <- cq_violin_median %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[4]],
    y = ~ Median[cq_dataframe$size == xaxis_vector[4]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[4]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    #color = I("red"),
    name = xaxis_vector[4]
  )
cq_violin_median <- cq_violin_median %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[5]],
    y = ~ Median[cq_dataframe$size == xaxis_vector[5]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[5]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    #color = I("purple"),
    name = xaxis_vector[5]
  )

cq_violin_median <- cq_violin_median %>%
  layout(
    xaxis = list(
      title = "Contig size",
      type = "category"
    ),
    yaxis = list(
      title = "Median ratio of windows",
      zeroline = T
    ),
    title = list(
      text = paste(
        c("Identification of male VS female (CQ)", species_name, dataset_type)
      )
    )
  )



cq_violin_lnmedian <- cq_dataframe %>%
  plot_ly(
    type = "violin"
  )

cq_violin_lnmedian <- cq_violin_lnmedian %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[1]],
    y = ~ lnMedian[cq_dataframe$size == xaxis_vector[1]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[1]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    # color = I("blue"),
    name = xaxis_vector[1]
  )
cq_violin_lnmedian <- cq_violin_lnmedian %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[2]],
    y = ~ lnMedian[cq_dataframe$size == xaxis_vector[2]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[2]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    # color = I("yellow"),
    name = xaxis_vector[2]
  )
cq_violin_lnmedian <- cq_violin_lnmedian %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[3]],
    y = ~ lnMedian[cq_dataframe$size == xaxis_vector[3]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[3]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    #color = I("green"),
    name = xaxis_vector[3]
  )
cq_violin_lnmedian <- cq_violin_lnmedian %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[4]],
    y = ~ lnMedian[cq_dataframe$size == xaxis_vector[4]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[4]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    #color = I("red"),
    name = xaxis_vector[4]
  )
cq_violin_lnmedian <- cq_violin_lnmedian %>%
  add_trace(
    x = ~ size[cq_dataframe$size == xaxis_vector[5]],
    y = ~ lnMedian[cq_dataframe$size == xaxis_vector[5]],
    type = 'violin',
    box = list(visible = T),
    meanline = list(visible = T),
    marker = list(opacity = 0.2),
    opacity = 0.5,
    text = ~ CONTIG_SIZE[cq_dataframe$size == xaxis_vector[5]],
    points = "all",
    pointpos = 0,
    jitter = 1,
    #color = I("purple"),
    name = xaxis_vector[5]
  )

cq_violin_lnmedian <- cq_violin_lnmedian %>%
  layout(
    xaxis = list(
      title = "Contig size",
      type = "category"
    ),
    yaxis = list(
      title = "lnMedian ratio of windows",
      zeroline = T
    ),
    title = list(
      text = paste(
      c("Identification of male VS female (CQ)", species_name, dataset_type)
      )
    )
  )


dir.create(file.path(getwd(), "figures"), showWarnings = FALSE)

htmlwidgets::saveWidget(as_widget(cq_violin_median), paste0(file.path(getwd(), "figures"), "/cq_", output_file_base, "_Median_w", window_size, ".html"), selfcontained = T, libdir = "library")

htmlwidgets::saveWidget(as_widget(cq_violin_lnmedian), paste0(file.path(getwd(), "figures"), "/cq_", output_file_base, "_lnMedian_w", window_size,  ".html"), selfcontained = T, libdir = "library")


# if we have contig descriptions, then produce these graphs

if (length(description_categories) > 1){
  
  cq_violin_median <- cq_dataframe %>%
    plot_ly(
      type = "violin",
      x = ~ CONTIG_DESCRIPTION,
      y = ~ Median,
      split = ~ CONTIG_DESCRIPTION,
      box = list(visible = T),
      meanline = list(visible = T),
      marker = list(opacity = 0.2),
      opacity = 0.5,
      text = ~ CONTIG_DESCRIPTION,
      points = "all",
      pointpos = 0,
      jitter = 1
    )
  

  cq_violin_median <- cq_violin_median %>%
    layout(
      xaxis = list(
        title = "Contig category",
        type = "category"
      ),
      yaxis = list(
        title = "Median ratio of windows",
        zeroline = T
      ),
      title = list(
        text = paste(
          c("Identification of male VS female (CQ)", species_name, dataset_type)
        )
      )
    )
  
  cq_violin_lnmedian <- cq_dataframe %>%
    plot_ly(
      type = "violin",
      x = ~ CONTIG_DESCRIPTION,
      y = ~ lnMedian,
      split = ~ CONTIG_DESCRIPTION,
      box = list(visible = T),
      meanline = list(visible = T),
      marker = list(opacity = 0.2),
      opacity = 0.5,
      text = ~ CONTIG_DESCRIPTION,
      points = "all",
      pointpos = 0,
      jitter = 1
    )
  
   cq_violin_lnmedian <- cq_violin_lnmedian %>%
    layout(
      xaxis = list(
        title = "Contig size",
        type = "category"
      ),
      yaxis = list(
        title = "lnMedian ratio of windows",
        zeroline = T
      ),
      title = list(
        text = paste(
          c("Identification of male VS female (CQ)", species_name, dataset_type)
        )
      )
    )
  
  
  dir.create(file.path(getwd(), "figures"), showWarnings = FALSE)
  
  htmlwidgets::saveWidget(as_widget(cq_violin_median), paste0(file.path(getwd(), "figures"), "/cq_", output_file_base, "_Median_w", window_size, "_description.html"), selfcontained = T, libdir = "library")
  
  htmlwidgets::saveWidget(as_widget(cq_violin_lnmedian), paste0(file.path(getwd(), "figures"), "/cq_", output_file_base, "_lnMedian_w", window_size,  "_description.html"), selfcontained = T, libdir = "library")
  
}
  


