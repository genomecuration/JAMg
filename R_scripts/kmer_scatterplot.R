female_data <-read.csv(kmer_female_input,header=F,sep="\t")
male_data <-read.csv(kmer_male_input,header=F,sep="\t")
#id_data <- read.csv(id_file, header = F, sep = "\t")

kmer_dataframe <- merge(female_data, male_data, by=c("V1","V2"),all=T)
#kmer_dataframe <- merge(kmer_dataframe, id_data, by="V1", all = T)
kmer_dataframe[is.na(kmer_dataframe)]<-0
#kmer_dataframe <- filter(kmer_dataframe, V2 >= 10000, .preserve = TRUE)

colnames(kmer_dataframe) <- c("ID", "Size", "Female", "Male")
remaining_ind <- which(abs(kmer_dataframe$Male - kmer_dataframe$Female) > 0.05)
kmer_dataframe_sub <- kmer_dataframe[remaining_ind,]

kmer_dataframe_sub$Ratio <-kmer_dataframe_sub$Female / kmer_dataframe_sub$Male

fig_all <- kmer_dataframe %>%
  plot_ly(
    type = 'scatter',
    mode = 'markers',
    x = ~Female,
    y = ~Male,
    size = ~log10(Size),
    color = "red",
    alpha = 0.2,
    text = ~ID,
    hovertemplate = ~paste(
      "<b>%{text}</b><br><br>",
      "%{yaxis.title.text}: %{y}<br>",
      "%{xaxis.title.text}: %{x}<br>",
      "Size: ",
      Size,
      #"<br>Description: ",
      #Description,
      "<extra></extra>"
    )
  ) 

fig_all <- fig_all %>%
  layout(
    title = paste(
      c("Identification of male VS female (kmer)", species_name, dataset_type)
    )
  )

fig_sub <- kmer_dataframe_sub %>%
  plot_ly(
    type = 'scatter',
    mode = 'markers',
    x = ~Female,
    y = ~Male,
    size = ~log10(Size),
    color = "red",
    alpha = 0.2,
    text = ~ID,
    hovertemplate = ~paste(
      "<b>%{text}</b><br><br>",
      "%{yaxis.title.text}: %{y}<br>",
      "%{xaxis.title.text}: %{x}<br>",
      "Size: ",
      Size,
      # "<br>Description: ",
      #Description, 
      "<extra></extra>"
    )
  ) 

fig_sub <- fig_sub %>%
  layout(
    title = paste(
      c("Identification of male VS female (kmer) subset", species_name, dataset_type)
    )
  )

dir.create(file.path(getwd(), "figures"), showWarnings = FALSE)

htmlwidgets::saveWidget(as_widget(fig_all), paste0(file.path(getwd(), "figures"), "/kmer_", output_file_base, ".html"), selfcontained = T, libdir = "library")

htmlwidgets::saveWidget(as_widget(fig_sub), paste0(file.path(getwd(), "figures"), "/kmer_", output_file_base, "_subset.html"), selfcontained = T, libdir = "library")