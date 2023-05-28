library(tidyverse, quietly = T)
library(magrittr)
library(ggpubr)
library(ComplexUpset)
library(enrichR)
library(eulerr)


setEnrichrSite("Enrichr")


### Color palette

cRed <- "#f94144"
cOrange1 <- "#f3722c"
cOrange2 <- "#f8961e"
cOrange3 <- "#f9844a"
cYellow <- "#f9c74f"
cGreen1 <- "#90be6d"
cGreen2 <- "#43aa8b"
cGreen3 <- "#4d908e"
cBlue <- "#577590"
cTeal <- "#277da1"



### Processing MMCC data

readQuantReportsMMCC <- function(fileRoot1, ...) {


  LOQ_report <-  read.csv(paste0(fileRoot1, "_QuantOptimization.csv"))

  quant_unrefined <- read.csv(paste0(fileRoot1, "_QuantResults.csv"))  %>%
    mutate(Total.Area.Fragment = gsub("#N/A", 0, Total.Area.Fragment)) %>%
    mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
    mutate(FileRoot = fileRoot1) %>%
    mutate(MSMethod = paste0(Analyzer, ", ", IsolationWindow, " Th, ",
                             InjectionTime, " ms, ", type))%>%
    mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
    mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
    mutate(Transitions = "All transitions")%>%
    mutate(SampleLabel = paste0(Analyte.Concentration, "%")) %>%
    full_join(select(LOQ_report, Peptide, Original.Limit.of.Quantification), by = "Peptide")

  quant_refined <- read.csv(paste0(fileRoot1, "_QuantResultsRefined.csv"))  %>%
    mutate(Total.Area.Fragment = gsub("#N/A", 0, Total.Area.Fragment)) %>%
    mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
    mutate(FileRoot = paste0(fileRoot1, "_Refined")) %>%
    mutate(MSMethod = paste0(Analyzer, ", ", IsolationWindow, " Th, ",
                             InjectionTime, " ms, ", type))%>%
    mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
    mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
    mutate(Transitions = "LOQ optimized transitions") %>%
    mutate(SampleLabel = paste0(Analyte.Concentration, "%")) %>%
    full_join(select(LOQ_report, Peptide, Optimized.Limit.of.Quantification), by = "Peptide")

  colnames(quant_unrefined)[colnames(quant_unrefined) == 'Original.Limit.of.Quantification'] <- 'LOQ'

  colnames(quant_refined)[colnames(quant_refined) == 'Optimized.Limit.of.Quantification'] <- 'LOQ'

  quant_full <- bind_rows(quant_unrefined, quant_refined) %>%
    arrange(Analyte.Concentration) %>%
    mutate(SampleLabel = factor(SampleLabel, unique(SampleLabel)))

  files <- list(...)

  for (f in files) {

    LOQ_report <-  read.csv(paste0(f, "_QuantOptimization.csv"))

    quant_unrefined <- read.csv(paste0(f, "_QuantResults.csv"))  %>%
      mutate(Total.Area.Fragment = gsub("#N/A", 0, Total.Area.Fragment)) %>%
      mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
      mutate(FileRoot = f) %>%
      mutate(MSMethod = paste0(Analyzer, ", ", IsolationWindow, " Th, ",
                               InjectionTime, " ms, ", type))%>%
      mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
      mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
      mutate(Transitions = "All transitions")%>%
      mutate(SampleLabel = paste0(Analyte.Concentration, "%")) %>%
      full_join(select(LOQ_report, Peptide, Original.Limit.of.Quantification), by = "Peptide")

    quant_refined <- read.csv(paste0(f, "_QuantResultsRefined.csv"))  %>%
      mutate(Total.Area.Fragment = gsub("#N/A", 0, Total.Area.Fragment)) %>%
      mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
      mutate(FileRoot = paste0(f, "_Refined")) %>%
      mutate(MSMethod = paste0(Analyzer, ", ", IsolationWindow, " Th, ",
                               InjectionTime, " ms, ", type))%>%
      mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
      mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
      mutate(Transitions = "LOQ optimized transitions") %>%
      mutate(SampleLabel = paste0(Analyte.Concentration, "%")) %>%
      full_join(select(LOQ_report, Peptide, Optimized.Limit.of.Quantification), by = "Peptide")

    colnames(quant_unrefined)[colnames(quant_unrefined) == 'Original.Limit.of.Quantification'] <- 'LOQ'

    colnames(quant_refined)[colnames(quant_refined) == 'Optimized.Limit.of.Quantification'] <- 'LOQ'

    quant_full <- bind_rows(quant_full, quant_unrefined, quant_refined)
  }

  quant_full <- quant_full %>%
    drop_na(LOQ) %>%
    drop_na(Transitions)
  
  return(quant_full)
}


plotCVsMMCCAll <- function(QuantReports, TransSet, overlap = FALSE,
                           file1 = "FILE1", file2 = "FILE2") {

  if (overlap == TRUE) {

    sharedPeps <- semi_join(filter(QuantReports,FileRoot == file1),
                            filter(QuantReports,FileRoot == file2),
                            by = "Peptide")

    QuantReports <- QuantReports %>%
      filter(FileRoot == file1 | FileRoot == file2) %>%
      filter(Peptide %in% sharedPeps$Peptide)


  }

  peptideSummary <- QuantReports %>%
    filter(Transitions == TransSet) %>%
    group_by(Description, Peptide,
             SampleLabel, Analyzer, IsolationWindow, InjectionTime,
             MSMethod, Precursor.Mz, type, Transitions, Analyte.Concentration) %>%
    summarise( NormCV = sd(as.numeric(Total.Area.Fragment))/mean(as.numeric(Total.Area.Fragment))*100)%>%
    ungroup() %>%
    arrange(Analyzer, desc(type),InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    arrange(desc(Analyte.Concentration)) %>%
    mutate(SampleLabel = factor(SampleLabel, unique(SampleLabel)))


  peptideSummary_summary <- peptideSummary %>%
    ungroup() %>%
    group_by(Description, Analyzer, SampleLabel, MSMethod, Transitions, Analyte.Concentration) %>%
    summarise(medCV = median(NormCV, na.rm = TRUE),
              n = n())



  ggplot(peptideSummary,
         aes(x = NormCV, fill = MSMethod)) +
    geom_histogram(color = "black", bins = 50) +
    scale_fill_manual(values = c(cRed, cOrange1, cYellow, cTeal)) +
    facet_grid(SampleLabel~MSMethod) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid.minor = element_blank())+
    xlab("% CV")+
    ylab("Peptides")+
    theme(strip.text = element_text(
      size = 8)) +
    geom_vline(data = peptideSummary_summary,
               aes(xintercept = medCV), linetype = "dashed", color = "black")+
    geom_text(data = peptideSummary_summary,
              aes(x = 60, label = paste0("Median CV = \n", round(medCV, 1), "%")),
              y = Inf, vjust = 2, size = 3.5) +
    coord_cartesian(xlim = c(0,100))

}


plotCVsMMCCProteins <- function(QuantReports, TransSet, overlap = FALSE,
                                file1 = "FILE1", file2 = "FILE2") {

  if (overlap == TRUE) {

    sharedPeps <- semi_join(filter(QuantReports,FileRoot == file1),
                            filter(QuantReports,FileRoot == file2),
                            by = "Protein")

    QuantReports <- QuantReports %>%
      filter(FileRoot == file1 | FileRoot == file2) %>%
      filter(Protein %in% sharedPeps$Protein)


  }


  peptideSummary <- QuantReports %>%
    filter(Transitions == TransSet)%>%
    group_by(Replicate, Description, Protein,
             SampleLabel, Analyzer, IsolationWindow, InjectionTime,
             MSMethod, type, Transitions, Analyte.Concentration) %>%
    summarize(Total.Area.Fragment = sum(as.numeric(Total.Area.Fragment))) %>%
    ungroup() %>%
    group_by(Description, Protein,
             SampleLabel, Analyzer, IsolationWindow, InjectionTime,
             MSMethod, type, Transitions, Analyte.Concentration) %>%
    summarise( NormCV = sd(as.numeric(Total.Area.Fragment))/mean(as.numeric(Total.Area.Fragment))*100)%>%
    ungroup() %>%
    arrange(Analyzer, desc(type), InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    arrange(desc(Analyte.Concentration)) %>%
    mutate(SampleLabel = factor(SampleLabel, unique(SampleLabel)))


  peptideSummary_summary <- peptideSummary %>%
    ungroup() %>%
    group_by(Description, Analyzer, SampleLabel, MSMethod, Transitions, Analyte.Concentration) %>%
    summarise(medCV = median(NormCV, na.rm = TRUE),
              n = n())



  ggplot(peptideSummary,
         aes(x = NormCV, fill = MSMethod)) +
    geom_histogram(color = "black", bins = 50) +
    scale_fill_manual(values = c(cRed, cOrange1, cYellow, cTeal)) +
    facet_grid(SampleLabel~MSMethod) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid.minor = element_blank())+
    xlab("% CV")+
    ylab("Proteins")+
    theme(strip.text = element_text(
      size = 8)) +
    geom_vline(data = peptideSummary_summary,
               aes(xintercept = medCV), linetype = "dashed", color = "black")+
    geom_text(data = peptideSummary_summary,
              aes(x = 60, label = paste0("Median CV = \n", round(medCV, 1), "%")),
              y = Inf, vjust = 2, size = 3.5) +
    coord_cartesian(xlim = c(0,100))


}


plotLOQMMCC <- function(QuantReports) {

  peptideSummary <- QuantReports %>%
    group_by(MSMethod, Peptide, Transitions,
             type, InjectionTime, IsolationWindow, Analyzer) %>%
    summarise( LOQ = mean(as.numeric(LOQ)))%>%
    ungroup() %>%
    filter(LOQ < 100) %>%
    arrange(Analyzer, desc(type), InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod)))


  peptideSummary_summary <- peptideSummary %>%
    ungroup() %>%
    filter(LOQ < 100) %>%
    group_by(MSMethod, Transitions) %>%
    summarise(medLOQ = median(LOQ, na.rm = TRUE),
              n = n())



  ggplot(peptideSummary,
         aes(x = LOQ, fill = MSMethod)) +
    geom_histogram(color = "black", bins = 30) +
    facet_grid(MSMethod~Transitions) +
    scale_fill_manual(values = c(cRed, cOrange1, cYellow, cTeal)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          legend.text=element_text(size=2), legend.title = element_text(size=8))+
    xlab("LLOQ (%)")+
    ylab("Peptides") +
    geom_vline(data = peptideSummary_summary,
               aes(xintercept = medLOQ), linetype = "dashed", color = "black")+
    geom_text(data = peptideSummary_summary,
              aes(x = 60, label = paste0(n,
                                         " quantifiable peptides \n Median LLOQ = ", round(medLOQ, 0), "%")),
              y = Inf, vjust = 2, size = 3)+
    theme(strip.text = element_text(
      size = 8))

}




plotLOQMMCCSummary <- function(QuantReports,
                               OTmultiplier = 60,
                               ASmultiplier = 60,
                               AxisLabel = "Peptides") {
  
  peptideSummary <- QuantReports %>%
    group_by(MSMethod, Peptide, Transitions,
             type, InjectionTime, IsolationWindow, Analyzer) %>%
    summarise( LOQ = mean(as.numeric(LOQ)))%>%
    ungroup() %>%
    arrange(Analyzer, desc(type), InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    mutate(Quant = if_else(LOQ < 100, "Quantitative", "Detectable"))%>%
    mutate(Quant = if_else(LOQ < 10, "10x dynamic range", Quant))%>%
    mutate(Quant = if_else(LOQ < 2, "50x dynamic range", Quant)) %>%
    mutate(Quant = factor(Quant, levels = c("Detectable", "Quantitative",
                                            "10x dynamic range",
                                            "50x dynamic range")))%>%
    group_by(Quant, MSMethod, Transitions, Analyzer) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    mutate(Count = if_else(Analyzer == "Orbitrap", Count/(OTmultiplier/60),
                           Count/(ASmultiplier/60)))
  
  
  ggplot(peptideSummary, 
         aes(alpha = Quant, y = Count, x = MSMethod)) +
    geom_col(color = "white", fill = "white", alpha = 1) +
    geom_col(fill = cBlue, color = "black") +
    scale_alpha_manual(values = c(0.25, 0.5, 0.75, 1),
                       name = "") +
    facet_wrap(~Transitions, nrow = 2) +
    theme_minimal() +
    ylab(AxisLabel)+
    theme(strip.text = element_text(
      size = 9), legend.text=element_text(size=8)) +
    theme(legend.position = "right", axis.title.x = element_blank(),
          axis.text.x = element_text(vjust = 1, hjust=1, size = 8, angle = 45))
  
  
}


writeLOQMMCCSummary <- function(QuantReports) {
  
  peptideSummary <- QuantReports %>%
    group_by(MSMethod, Peptide, Transitions,
             type, InjectionTime, IsolationWindow, Analyzer) %>%
    summarise( LOQ = mean(as.numeric(LOQ)))%>%
    ungroup() %>%
    arrange(Analyzer, desc(type), InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    mutate(Quant = if_else(LOQ < 100, "Quantitative", "Detectable"))%>%
    mutate(Quant = if_else(LOQ < 10, "10x dynamic range", Quant))%>%
    mutate(Quant = if_else(LOQ < 2, "50x dynamic range", Quant)) %>%
    mutate(Quant = factor(Quant, levels = c("Detectable", "Quantitative",
                                            "10x dynamic range",
                                            "50x dynamic range"))) %>%
    drop_na(LOQ) %>%
    drop_na(Transitions) %>%
    group_by(MSMethod, Transitions, Quant) %>%
    summarize(Peptides = n())  %>%
    mutate(PeptidesTotal = if_else(Quant == "Detectable", 
                               Peptides + lead(Peptides,3) +
                                 lead(Peptides,2) +
                                 lead(Peptides,1), Peptides))%>%
    mutate(PeptidesTotal = if_else(Quant == "Quantitative", 
                               Peptides + 
                                 lead(Peptides,2) +
                                 lead(Peptides,1), PeptidesTotal))%>%
    mutate(PeptidesTotal = if_else(Quant == "10x dynamic range", 
                               Peptides + 
                                 lead(Peptides,1), PeptidesTotal)) %>%
    ungroup() %>%
    select(-Peptides) %>%
    pivot_wider(names_from = Quant, values_from = PeptidesTotal)
  
  
  return(peptideSummary)
  
  
}



plotLOQMMCCHistogram <- function(QuantReports,
                                 Fileroot1, Label1,
                                 Fileroot2, Label2) {

  file1 <- QuantReports %>%
    filter(FileRoot == Fileroot1) %>%
    group_by(MSMethod, Peptide, Transitions,
             type, InjectionTime, IsolationWindow, Analyzer) %>%
    summarise( LOQ = mean(as.numeric(LOQ)))%>%
    ungroup() %>%
    filter(LOQ < 100) %>%
    select(Peptide, LOQ)%>%
    set_colnames(c("Peptide",  "LOQ.x")) %>%
    rbind(data.frame(Peptide = "test", LOQ.x = 1))

  file2 <- QuantReports %>%
    filter(FileRoot == Fileroot2) %>%
    group_by(MSMethod, Peptide, Transitions,
             type, InjectionTime, IsolationWindow, Analyzer) %>%
    summarise( LOQ = mean(as.numeric(LOQ)))%>%
    ungroup() %>%
    filter(LOQ < 100) %>%
    select(Peptide, LOQ) %>%
    set_colnames(c("Peptide",  "LOQ.y")) %>%
    rbind(data.frame(Peptide = "test", LOQ.y = 1))

  pairwise <- inner_join(file1, file2, by = "Peptide") %>%
    mutate(ratio = LOQ.x/LOQ.y) %>%
    mutate(Better = if_else(ratio > 1, Label2, Label1)) %>%
    mutate(Better = if_else(ratio == 1, "Tie", Better)) %>%
    mutate(Better = factor(Better, levels = c(Label1, "Tie", Label2))) %>%
    mutate(RAT = ratio) %>%
    mutate(ratio = log10(ratio)) %>%
    filter(abs(ratio) < 100)


  ggplot() +
    geom_histogram(data = pairwise, color = "black", bins = 50, aes(x = ratio, fill = Better)) +
    theme_minimal() +
    xlab(paste0("log10(LOQ in ", Label1, " / LOQ in ", Label2, ")"))+
    ylab("Peptides") +
    scale_fill_manual(values = c(cRed, "grey", cTeal),
                      name = "Better LLOQ in: ") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 0.8) +
    geom_vline(xintercept = median(pairwise$ratio), linetype = "dashed", color = "black", size = 1.2) +
    theme(legend.position = "bottom",
          legend.text=element_text(size=8), legend.title = element_text(size=9),
          axis.title.x = element_text(size = 8))+
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
    geom_text(y = Inf, vjust = 1, x = 1.3, 
              aes(label = paste0("Shared peptides = ", nrow(pairwise),
                                 "\n Median Ratio = ", round(median(pairwise$RAT), 2))))

}


upsetPlot <- function(QuantReports, Grouplevel) {


  if (Grouplevel == "Protein") {
    listInput <- QuantReports %>%
      select(MSMethod, Protein) %>%
      unique() %>%
      ungroup() %>%
      reshape2::dcast(Protein ~ MSMethod, fun.aggregate = function(x) 1L, fill = 0L)

    Methods <- c(unique(QuantReports$MSMethod))


    upset(listInput, Methods, name = "Method",
          width_ratio=0.1,
          base_annotations=list(
            'Intersection size'=upset_annotate('..count..',
                                               geom_bar(fill='#264653', color = "#264653"))+ ylab("Shared proteins")),
          matrix=(
            intersection_matrix(geom=geom_point(shape='circle filled', size=3)) +  scale_color_manual(
              values=c("Astral, 4 Th, 10 ms, DIA" = cOrange1,
                       "Astral, 2 Th, 3.5 ms, DIA" = cRed,
                       "Astral, 2 Th, 10 ms, dDIA" = cYellow,
                       "Orbitrap, 8 Th, 22 ms, DIA" = cTeal,
                       "Orbitrap, 2 Th, 23 ms, DIA" = cGreen2), guide = "none"
            )
          ), queries = list(upset_query(set = "Astral, 4 Th, 10 ms, DIA", fill = cOrange1),
                            upset_query(set = "Astral, 2 Th, 3.5 ms, DIA", fill = cRed),
                            upset_query(set = "Astral, 2 Th, 10 ms, dDIA", fill = cYellow),
                            upset_query(set = "Orbitrap, 8 Th, 22 ms, DIA", fill = cTeal),
                            upset_query(set = "Orbitrap, 2 Th, 23 ms, DIA", fill = cGreen2)),
          set_sizes=(
            upset_set_size()
            + scale_y_continuous(breaks = c(7000, 0))))
  }

  else {

    listInput <- QuantReports %>%
      select(MSMethod, Peptide) %>%
      unique() %>%
      ungroup() %>%
      reshape2::dcast(Peptide ~ MSMethod, fun.aggregate = function(x) 1L, fill = 0L)

    Methods <- c(unique(QuantReports$MSMethod))


    upset(listInput, Methods, name = "Method",
          width_ratio=0.1,
          base_annotations=list(
            'Intersection size'=upset_annotate('..count..',
                                               geom_bar(fill='#264653', color = "#264653"))+ ylab("Shared peptides")),
          matrix=(
            intersection_matrix(geom=geom_point(shape='circle filled', size=3)) +  scale_color_manual(
              values=c("Astral, 4 Th, 10 ms, DIA" = cOrange1,
                       "Astral, 2 Th, 3.5 ms, DIA" = cRed,
                       "Astral, 2 Th, 10 ms, dDIA" = cYellow,
                       "Orbitrap, 8 Th, 22 ms, DIA" = cTeal,
                       "Orbitrap, 2 Th, 23 ms, DIA" = cGreen2), guide = "none"
            )
          ), queries = list(upset_query(set = "Astral, 4 Th, 10 ms, DIA", fill = cOrange1),
                            upset_query(set = "Astral, 2 Th, 3.5 ms, DIA", fill = cRed),
                            upset_query(set = "Astral, 2 Th, 10 ms, dDIA", fill = cYellow),
                            upset_query(set = "Orbitrap, 8 Th, 22 ms, DIA", fill = cTeal),
                            upset_query(set = "Orbitrap, 2 Th, 23 ms, DIA", fill = cGreen2)),
          set_sizes=(
            upset_set_size()
            + scale_y_continuous(breaks = c(70000, 0)))
    )
  }



}


quantRatioDensity <- function(QuantReports, Concentration1, Concentration2, Alpha = 0.025,
                              multiplier = 130, short = FALSE) {

  if (short == TRUE) {
    cOrange1 = "#277da1"

  }

  results <- QuantReports %>%
    filter(Transitions == "LOQ optimized transitions") %>%
    arrange(Analyzer, desc(type), InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    filter(Analyte.Concentration == Concentration1 |
             Analyte.Concentration == Concentration2 |
             Analyte.Concentration == 0) %>%
    select(Peptide, Analyte.Concentration,MSMethod, Total.Area.Fragment) %>%
    mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment)) %>%
    group_by(Peptide, MSMethod, Analyte.Concentration) %>%
    summarize(Total.Area.Fragment = mean(Total.Area.Fragment)) %>%
    ungroup() %>%
    mutate(Conc = Analyte.Concentration) %>%
    mutate(Analyte.Concentration = if_else(Analyte.Concentration == Concentration1, "Sample1", "Sample2"))%>%
    mutate(Analyte.Concentration = if_else(Conc == 0, "Blank", Analyte.Concentration)) %>%
    select(-Conc) %>%
    pivot_wider(names_from = Analyte.Concentration, values_from = Total.Area.Fragment) %>%
    mutate(Ratio = (Sample2-Blank)/(Sample1-Blank)) %>%
    mutate(Group = "A") %>%
    group_by(MSMethod) %>%
    mutate(Sample1 = Sample1/max(Sample1)*100)


  labs <- c("")

  names(labs) <- c("A")

  scatter_plot <- ggplot(results, aes(x = log10(Sample1), y = log10(Ratio))) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    geom_bin2d(bins = 200) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    facet_grid(MSMethod~Group, labeller = labeller(Group = labs)) +
    theme(legend.position = c(0.9, 0.75), panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9, color = "white"),
          legend.title.align=0.5, legend.background = element_rect(fill = "white", color = "black",
                                                                   linetype = "solid")) +
    coord_cartesian(ylim = c(log10(Concentration2/Concentration1/multiplier),
                             log10(Concentration2/Concentration1* multiplier)),
                    xlim = c(-4.2, 2.2)) +
    ylab(paste0("log10(Signal in ", Concentration2, "% / Signal in ", Concentration1, "%)")) +
    xlab("log10(Relative Abundance)") +
    geom_hline(yintercept = log10(Concentration2/Concentration1), linetype = "dashed") +
    scale_fill_gradientn(colors = c("white", cBlue, cRed, cYellow), name = "Peptides")


  density_plot <- ggplot(results, aes(y = log10(Ratio), color = MSMethod)) +
    geom_density(size = 1) +
    theme_minimal() +
    geom_hline(yintercept = log10(Concentration2/Concentration1), linetype = "dashed") +
    facet_grid(MSMethod~Group, labeller = labeller(Group = labs)) +
    theme(legend.position = "none", axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9))  +
    scale_color_manual(values = c(cRed, cOrange1, cYellow, cTeal)) +
    coord_cartesian(ylim = c(log10(Concentration2/Concentration1/multiplier),
                             log10(Concentration2/Concentration1* multiplier)),
                    xlim = c(0,2.1)) +
    xlab("Density") +
    scale_x_continuous(breaks = seq(0, 10, 2))

  ggarrange(scatter_plot, density_plot,
            widths = c(1, 0.35))


}


quantRatioDensityProtein <- function(QuantReports, Concentration1, Concentration2, Alpha = 0.025,
                                     multiplier = 125, short = FALSE) {

  if (short == TRUE) {
    cOrange1 = "#277da1"

  }

  results <- QuantReports %>%
    filter(Transitions == "LOQ optimized transitions") %>%
    arrange(Analyzer, desc(type), InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    filter(Analyte.Concentration == Concentration1 |
             Analyte.Concentration == Concentration2|
             Analyte.Concentration == 0) %>%
    select(Protein, Analyte.Concentration,MSMethod, Total.Area.Fragment, Replicate)%>%
    group_by(Protein, Analyte.Concentration,MSMethod, Total.Area.Fragment, Replicate) %>%
    summarize(Total.Area.Fragment = sum(as.numeric(Total.Area.Fragment))) %>%
    ungroup() %>%
    group_by(Protein, MSMethod, Analyte.Concentration) %>%
    summarize(Total.Area.Fragment = mean(Total.Area.Fragment)) %>%
    ungroup()  %>%
    mutate(Conc = Analyte.Concentration) %>%
    mutate(Analyte.Concentration = if_else(Analyte.Concentration == Concentration1, "Sample1", "Sample2"))%>%
    mutate(Analyte.Concentration = if_else(Conc == 0, "Blank", Analyte.Concentration)) %>%
    select(-Conc) %>%
    pivot_wider(names_from = Analyte.Concentration, values_from = Total.Area.Fragment) %>%
    mutate(Ratio = (Sample2-Blank)/(Sample1-Blank)) %>%
    mutate(Group = "A") %>%
    mutate(Sample1 = Sample1/max(Sample1)*100)


  labs <- c("")

  names(labs) <- c("A")

  scatter_plot <- ggplot(results, aes(x = log10(Sample1), y = log10(Ratio))) +
    geom_bin2d(bins = 200) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    facet_grid(MSMethod~Group, labeller = labeller(Group = labs)) +
    theme(legend.position = "left", panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9, color = "white"),
          legend.title.align=0.5, legend.background = element_rect(fill = "white", color = "black", linetype = "solid")) +
    coord_cartesian(ylim = c(log10(Concentration2/Concentration1/multiplier),
                             log10(Concentration2/Concentration1* multiplier))) +
    ylab(paste0("log10(Signal in ", Concentration2, "% / Signal in ", Concentration1, "%)")) +
    xlab("log10(Relative Abundance)") +
    geom_hline(yintercept = log10(Concentration2/Concentration1), linetype = "dashed") +
    scale_fill_gradientn(colors = c("white", cBlue, cRed, cYellow), name = "Proteins")


  density_plot <- ggplot(results, aes(y = log10(Ratio), color = MSMethod)) +
    geom_density(size = 1) +
    theme_minimal() +
    geom_hline(yintercept = log10(Concentration2/Concentration1), linetype = "dashed") +
    facet_grid(MSMethod~Group, labeller = labeller(Group = labs)) +
    theme(legend.position = "none", axis.text.y=element_blank(),
          axis.ticks.y=element_blank(), axis.title.y = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 9))  +
    scale_color_manual(values = c(cRed, cOrange1, cYellow, cTeal)) +
    coord_cartesian(ylim = c(log10(Concentration2/Concentration1/multiplier),
                             log10(Concentration2/Concentration1* multiplier)),
                    xlim = c(0,6)) +
    xlab("Density")+
    scale_x_continuous(breaks = seq(0, 10, 2))


  ggarrange(scatter_plot, density_plot,
            widths = c(1, 0.35))


}


plotLOQMMCCSummaryGraphicalAbstract <- function(QuantReports,
                               OTmultiplier = 60,
                               ASmultiplier = 60,
                               AxisLabel = "Peptides") {
  
  peptideSummary <- QuantReports %>%
    group_by(MSMethod, Peptide, Transitions,
             type, InjectionTime, IsolationWindow, Analyzer) %>%
    filter(Transitions == "LOQ optimized transitions") %>%
    summarise( LOQ = mean(as.numeric(LOQ)))%>%
    ungroup() %>%
    arrange(Analyzer, desc(type), InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    mutate(Quant = if_else(LOQ < 100, "Quantitative", "Detectable"))%>%
    mutate(Quant = if_else(LOQ < 10, "10x dynamic range", Quant))%>%
    mutate(Quant = if_else(LOQ < 2, "50x dynamic range", Quant)) %>%
    mutate(Quant = factor(Quant, levels = c("Detectable", "Quantitative",
                                            "10x dynamic range",
                                            "50x dynamic range")))%>%
    filter(Quant != "Detectable") %>%
    group_by(MSMethod, Transitions, Analyzer) %>%
    summarize(Count = n()) %>%
    ungroup() %>%
    mutate(Count = if_else(Analyzer == "Orbitrap", Count/(OTmultiplier/60),
                           Count/(ASmultiplier/60))) %>%
    group_by(Analyzer) %>%
    summarize(Count = max(Count))
  
  
  ggplot(peptideSummary, 
         aes(x = Count, y = Analyzer, fill = Analyzer)) +
    geom_col(color = "black") +
    theme_minimal() +
    xlab(AxisLabel)+
    theme(axis.title.y = element_blank()) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c(cRed, cTeal)) +
    scale_x_continuous(breaks = c(50000, 100000))
  
  
}

### Processing plasma data

readQuantReports <- function(fileRoot1, ...) {


  quant_unrefined <- read.csv(paste0(fileRoot1, "_Quant.csv"))  %>%
    mutate(Normalized.Area = gsub("#N/A", 0, Normalized.Area)) %>%
    mutate(Library = gsub(".*\\_","", fileRoot1)) %>%
    mutate(Search = paste0(Library, " unrefined")) %>%
    mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
    mutate(FileRoot = fileRoot1) %>%
    mutate(MSMethod = paste0(IsolationWindow, " Th, ",
                             InjectionTime, " ms, ", type))%>%
    mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
    mutate(Total.Background.Fragment = as.numeric(Total.Background.Fragment)) %>%
    mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
    mutate(ResultType = "All")%>%
    mutate(Gradient = paste0(Gradient, " minutes"))

  quant_refined <- read.csv(paste0(fileRoot1, "_Refined_Quant.csv"))  %>%
    mutate(Normalized.Area = gsub("#N/A", 0, Normalized.Area)) %>%
    mutate(Library = gsub(".*\\_","", fileRoot1)) %>%
    mutate(Search = paste0(Library, " refined"))%>%
    mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
    mutate(FileRoot = paste0(fileRoot1, "_Refined")) %>%
    mutate(MSMethod = paste0(IsolationWindow, " Th, ",
                             InjectionTime, " ms, ", type))%>%
    mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
    mutate(Total.Background.Fragment = as.numeric(Total.Background.Fragment)) %>%
    mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
    mutate(ResultType = "Refined")%>%
    mutate(Gradient = paste0(Gradient, " minutes"))

  quant_full <- bind_rows(quant_unrefined, quant_refined)

  files <- list(...)

  for (f in files) {
    quant_unrefined <- read.csv(paste0(f, "_Quant.csv"))  %>%
      mutate(Normalized.Area = gsub("#N/A", 0, Normalized.Area)) %>%
      mutate(Library = gsub(".*\\_","", f)) %>%
      mutate(Search = paste0(Library, " unrefined"))  %>%
      mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
      mutate(FileRoot = f) %>%
      mutate(MSMethod = paste0(IsolationWindow, " Th, ",
                               InjectionTime, " ms, ", type))%>%
      mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
      mutate(Total.Background.Fragment = as.numeric(Total.Background.Fragment)) %>%
      mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
      mutate(ResultType = "All")%>%
      mutate(Gradient = paste0(Gradient, " minutes"))

    quant_refined <- read.csv(paste0(f, "_Refined_Quant.csv"))  %>%
      mutate(Normalized.Area = gsub("#N/A", 0, Normalized.Area)) %>%
      mutate(Library = gsub(".*\\_","", f)) %>%
      mutate(Search = paste0(Library, " refined"))%>%
      mutate(type = if_else(str_detect(Description, "dDIA"), "dDIA", "DIA")) %>%
      mutate(FileRoot = paste0(f, "_Refined")) %>%
      mutate(MSMethod = paste0(IsolationWindow, " Th, ",
                               InjectionTime, " ms, ", type))%>%
      mutate(Total.Area.Fragment = as.numeric(Total.Area.Fragment))%>%
      mutate(Total.Background.Fragment = as.numeric(Total.Background.Fragment)) %>%
      mutate(Best.Retention.Time = as.numeric(Best.Retention.Time)) %>%
      mutate(ResultType = "Refined")%>%
      mutate(Gradient = paste0(Gradient, " minutes"))

    quant_full <- bind_rows(quant_full, quant_unrefined, quant_refined)
  }


  return(quant_full)
}


plotCVs <- function(QuantReports, Transitions = "All") {

  peptideSummary <- QuantReports %>%
    filter(ResultType == Transitions) %>%
    group_by(Description, Search, Protein.Accession, Peptide,
             SampleLabel, Analyzer, IsolationWindow, InjectionTime,
             Gradient, MSMethod, Precursor.Mz, Precursor.Charge, type) %>%
    summarise( NormCV = sd(as.numeric(Normalized.Area))/mean(as.numeric(Normalized.Area))*100)%>%
    ungroup() %>%
    arrange(desc(type), IsolationWindow, InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod)))%>%
    filter(SampleLabel == "EV plasma")


  peptideSummary_summary <- peptideSummary %>%
    ungroup() %>%
    group_by(Description, Search, Analyzer, SampleLabel, MSMethod, Gradient) %>%
    summarise(medCV = median(NormCV, na.rm = TRUE))



  ggplot(peptideSummary,
         aes(y = NormCV, x = MSMethod, fill = MSMethod)) +
    #geom_violin(draw_quantiles = c( 0.5), color = "white") +
    geom_violin(color = "black") +
    geom_boxplot(color = "grey", fill = "white", width = 0.1, outlier.size = 0) +
    scale_fill_manual(values = c(cOrange3, cOrange2, cGreen2, cTeal, cYellow)) +
    geom_hline(yintercept = 20, linetype = "dashed") +
    facet_grid(cols = vars(Gradient), scales = "free", space = "free_x") +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
          axis.title.x = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(
            size = 9))+
    ylab("% CV") +
    labs(fill = "Method") +
    guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
    coord_cartesian(ylim = c(0, 50))


}


summarizeSearch <- function(QuantReports, Level) {

  resultsSummary <- QuantReports %>%
    group_by(Analyzer, MSMethod, Search, ResultType, Library,
             Gradient, IsolationWindow,type, InjectionTime) %>%
    summarize(Proteins = n_distinct(Protein.Preferred.Name),
              Peptides = n_distinct(Peptide))  %>%
    ungroup() %>%
    arrange(desc(type), IsolationWindow, InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod)))

  if (Level == "Protein") {
    ggplot(resultsSummary, aes(x = MSMethod, y = Proteins, fill = MSMethod, alpha = ResultType))+
      geom_bar(color = "white", stat = "Identity", position = "identity", fill = "white", alpha = 1) +
      geom_bar(color = "black", stat = "Identity", position = "identity") +
      theme_minimal() +
      facet_grid(cols = vars(Gradient), space = "free_x", scales = "free_x") +
      xlab("") +
      labs(fill = "Samples Searched") +
      scale_fill_manual(values = c(cOrange3, cOrange2, cGreen2, cTeal, cYellow)) +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))+
      ylab("Proteins") +
      scale_alpha_manual(values = c(0.5, 1)) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
            strip.text = element_text(size = 9))
  }

  else {
    ggplot(resultsSummary, aes(x = MSMethod, y = Peptides, fill = MSMethod, alpha = ResultType))+
      geom_bar(color = "white", stat = "Identity", position = "identity", fill = "white", alpha = 1) +
      geom_bar(color = "black", stat = "Identity", position = "identity") +
      theme_minimal() +
      facet_grid(cols = vars(Gradient), space = "free_x", scales = "free_x") +
      xlab("") +
      labs(fill = "Samples Searched") +
      scale_fill_manual(values = c(cOrange3, cOrange2, cGreen2, cTeal, cYellow)) +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))+
      ylab("Peptides") +
      scale_alpha_manual(values = c(0.5, 1)) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
            strip.text = element_text(size = 9))
  }
}



summarizeSearch2 <- function(QuantReports, Level) {

  resultsSummary <- QuantReports %>%
    group_by(Analyzer, MSMethod, Search, ResultType, Library,
             Gradient, IsolationWindow,type, InjectionTime) %>%
    summarize(Proteins = n_distinct(Protein.Preferred.Name),
              Peptides = n_distinct(Peptide))  %>%
    ungroup() %>%
    arrange(desc(type), IsolationWindow, InjectionTime) %>%
    mutate(MSMethod = factor(MSMethod, unique(MSMethod))) %>%
    filter(ResultType == "All")

  if (Level == "Protein") {
    ggplot(resultsSummary, aes(x = MSMethod, y = Proteins, fill = MSMethod))+
      geom_bar(color = "white", stat = "Identity", position = "identity", fill = "white", alpha = 1) +
      geom_bar(color = "black", stat = "Identity", position = "identity", alpha = 1) +
      theme_minimal() +
      facet_grid(cols = vars(Gradient), space = "free_x", scales = "free_x") +
      xlab("") +
      labs(fill = "Samples Searched") +
      scale_fill_manual(values = c(cOrange3, cOrange2, cGreen2, cTeal, cYellow)) +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))+
      ylab("Proteins") +
      #scale_alpha_manual(values = c(0.5, 1)) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
            strip.text = element_text(size = 9))
  }

  else {
    ggplot(resultsSummary, aes(x = MSMethod, y = Peptides, fill = MSMethod))+
      geom_bar(color = "white", stat = "Identity", position = "identity", fill = "white", alpha = 1) +
      geom_bar(color = "black", stat = "Identity", position = "identity", alpha = 1) +
      theme_minimal() +
      facet_grid(cols = vars(Gradient), space = "free_x", scales = "free_x") +
      xlab("") +
      labs(fill = "Samples Searched") +
      scale_fill_manual(values = c(cOrange3, cOrange2, cGreen2, cTeal, cYellow)) +
      guides(fill = guide_legend(title.position = "top", title.hjust = 0.5))+
      ylab("Peptides") +
      #scale_alpha_manual(values = c(0.5, 1)) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
            strip.text = element_text(size = 9))
  }
}


calcFoldChange <- function(FileRoot1, QuantReports, EnrichedList, DepletedList) {

  FCCalc <- filter(QuantReports, FileRoot == FileRoot1)  %>%
    group_by(Replicate, Protein.Accession, Protein.Preferred.Name, SampleLabel, Protein.Gene) %>%
    summarise(NormArea = sum(as.numeric(Normalized.Area)))%>%
    ungroup()  %>%
    group_by(Protein.Gene) %>%
    summarize(PValue = t.test(NormArea[str_detect(SampleLabel, "Total plasma")],
                              NormArea[str_detect(SampleLabel, "EV plasma")])$p.value,
              fold_change = mean(NormArea[1:3]) / mean(NormArea[4:6]),
              ProteinName = first(Protein.Preferred.Name),
              TP = mean(NormArea[str_detect(SampleLabel, "Total plasma")]),
              EV = mean(NormArea[str_detect(SampleLabel, "EV plasma")])) %>%
    mutate(Enriched = if_else(log2(fold_change) > 0.6 & PValue < 0.05,
                              "Enriched in EV prep", "No significant change"))%>%
    mutate(Enriched = if_else(log2(fold_change) < -0.6 & PValue < 0.05,
                              "Depleted in EV prep", Enriched)) %>%
    mutate(ProteinName = gsub("_HUMAN", "", ProteinName)) %>%
    mutate(Label = if_else(ProteinName %in% DepletedList |
                             ProteinName %in% EnrichedList, T, F)) %>%
    mutate(Predicted = if_else(ProteinName %in% EnrichedList,
                               "Enriched in EV prep", "No significant change"))%>%
    mutate(Predicted = if_else(ProteinName %in% DepletedList,
                               "Depleted in EV prep", Predicted)) %>%
    mutate(Rank_EV = rank(-EV))%>%
    mutate(Rank_TP = rank(-TP)) %>%
    mutate(EV = EV/sum(EV)*100) %>%
    mutate(TP = TP/sum(TP)*100)

  return(FCCalc)
}


plotFoldChange <- function(FCCalc, labels) {

  FCPlot <- ggplot(FCCalc, aes(x=log2(fold_change), y=-log10(PValue), col=Enriched)) +
    geom_point(alpha = 0.25, size = 0.5) +
    theme_minimal() +
    scale_color_manual(values=c(cBlue, cRed, "grey")) +
    geom_vline(xintercept=c(-0.6, 0.6), col="black") +
    geom_hline(yintercept=-log10(0.05), col="black") +
    xlab("log2(Fold Change)") +
    ylab("-log10(P value)")+
    theme(legend.position = "right") +
    guides(color = guide_legend(title.position = "top", title.hjust = 0.5,
                                override.aes = list(size=3, alpha = 1))) +
    scale_fill_manual(values=c(cBlue, cRed, "grey"))+
    coord_cartesian(ylim = c(0, 10), xlim = c(-10, 10))

  if (labels == TRUE) {
    FinalPlot <- FCPlot +
      geom_point(data = filter(FCCalc, Label == TRUE), size = 1.5, color = "black", shape = 1) +
      ggrepel::geom_label_repel(data = filter(FCCalc, Label == TRUE),
                                aes(x=log2(fold_change),
                                    y=-log10(PValue),
                                    label = ProteinName, color = Predicted),
                                segment.color = "black", size = 2,
                                min.segment.length = unit(0, 'lines'),
                                nudge_y = 0.5, show.legend = FALSE)

  } else {
    FinalPlot <- FCPlot
  }

  return(FinalPlot)
}


plotRankOrder <- function(FCCalc, labels) {


  EVrank <- ggplot(FCCalc, aes(x = Rank_EV, y = EV)) +
    geom_point(size = 0.7) +
    theme_minimal() +
    ylab("Percent of total signal") +
    xlab("Protein Rank") +
    ggtitle("EV enriched sample") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.8, 0.8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100))+
    scale_color_manual(values=c(cBlue, cRed, "grey"),
                       name = "") +
    coord_cartesian(ylim = c(1e-7, 50))

  TPrank <- ggplot(FCCalc, aes(x = Rank_TP, y = TP)) +
    geom_point(size = 0.7) +
    theme_minimal() +
    ylab("Percent of total signal") +
    xlab("Protein Rank") +
    ggtitle("Total plasma") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.8, 0.8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100))+
    scale_color_manual(values=c(cBlue, cRed, "grey"),
                       name = "") +
    coord_cartesian(ylim = c(1e-7, 50))

  if (labels == TRUE) {
    EVrank <- EVrank +
      geom_point(data = filter(FCCalc, Label == TRUE), size = 1.5, aes(color = Predicted), shape = 1) +
      ggrepel::geom_label_repel(data = filter(FCCalc, Label == TRUE),
                                aes(x=Rank_EV,
                                    y=EV,
                                    label = ProteinName, color = Predicted),
                                max.overlaps = Inf, show.legend = FALSE,
                                min.segment.length = unit(0, 'lines'))
    TPrank <- TPrank +
      geom_point(data = filter(FCCalc, Label == TRUE), size = 1.5, aes(color = Predicted), shape = 1) +
      ggrepel::geom_label_repel(data = filter(FCCalc, Label == TRUE),
                                aes(x=Rank_TP,
                                    y=TP,
                                    label = ProteinName, color = Predicted),
                                max.overlaps = Inf, show.legend = FALSE,
                                min.segment.length = unit(0, 'lines'))

  }
  rankPlot <- ggarrange(TPrank,EVrank, common.legend = TRUE, legend="bottom", nrow = 1)

  return(rankPlot)

}


enrichmentAnalysis <- function(FCCalc, dbs, Direction) {

  enrichedList <- filter(FCCalc, Enriched == Direction) %>%
    mutate(gene = Protein.Gene) %>%
    select(gene) %>%
    mutate(gene = strsplit(as.character(gene), " / ")) %>%
    unnest(gene)

  enrichedList <- enrichr(enrichedList$gene, dbs)

  plotEnrich(enrichedList[[1]], showTerms = 20, numChar = 50, y = "Count", orderBy = "P.value") +
    theme(legend.position = c(0.7, 0.4), legend.title.align=0.5,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ggtitle(paste0("Go Cellular Component Terms ", Direction)) +
    scale_fill_gradient(low = cBlue, high = cRed, name = "P-value",
                        labels=function(x) paste0(format(signif(x, 2)))) +
    labs(fill = "P-value") +
    xlab("")

}


plotRankOrder2 <- function(FCCalc, labels) {


  EVrank <- ggplot(data = FCCalc) +
    geom_point(size = 0.7, color = "#2d3142", aes(x = Rank_EV, y = EV))+
    theme_minimal() +
    ylab("Percent of total signal") +
    xlab("Protein Rank") +
    ggtitle("Enriched plasma") +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.8, 0.8),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100))+
    scale_color_manual(values=c(cRed, "grey"),
                       name = "") +
    coord_cartesian(ylim = c(1e-7, 50), xlim = c(-300, 5000))





  if (labels == TRUE) {

    df <- filter(FCCalc, Label == TRUE)

    labels <- rbind(data.frame(x = df$Rank_EV, y = df$EV, type = "EV Plasma" , ProteinName = df$ProteinName),
                    data.frame(x = df$Rank_TP, y = df$TP, type = "Total plasma", ProteinName = df$ProteinName ))

    EVrank <- EVrank  +
      ggrepel::geom_label_repel(data = filter(FCCalc, Label == TRUE),
                                aes(x=Rank_EV,
                                    y=EV,
                                    label = Protein.Gene, color = cRed),
                                box.padding = 1, force = 0.5, force_pull = 2.1,
                                max.overlaps = Inf, show.legend = FALSE,
                                min.segment.length = unit(0, 'lines')) +
      geom_point(data = filter(FCCalc, Label == TRUE),
                 aes(x=Rank_EV,
                     y=EV), color = cRed)


  }

  return(EVrank)

}




### Processing info from raw files

readRawInfo <- function(fileName, Analyzer, IsolationWindow, SampleType, maxIT, AGCtarget, MethodType) {

  rawFileInfo <- read.csv(fileName) %>%
    mutate(ions = InjectTime*total_S/1000) %>%
    mutate(Analyzer = Analyzer)%>%
    mutate(Ions2 = ifelse(Analyzer == "Orbitrap", rawFv, total_SN))%>%
    mutate(maxIT = maxIT) %>%
    mutate(Target = AGCtarget) %>%
    mutate(IsolationWindow = IsolationWindow) %>%
    mutate(SpeedHz = 1/(60*ScanTime - 60*lag(ScanTime)))%>%
    mutate(Sample = SampleType) %>%
    mutate(MethodType = MethodType) %>%
    mutate(Method = paste0(Analyzer, ", ", IsolationWindow,
                           " Th, ", maxIT, " ms, ",
                           MethodType))

  return(rawFileInfo)
}


plotIons <- function(rawFileInfo1, ...) {
  files <- list(...)
  for (f in files) rawFileInfo1 <- rbind(rawFileInfo1, f)

  rawFileInfo1 <- rawFileInfo1 %>%
    arrange(Analyzer, maxIT, desc(IsolationWindow)) %>%
    mutate(Method=factor(Method, levels = unique(Method)))

  rawFileInfoSummary <- select(rawFileInfo1, Target, Sample, Method) %>%
    unique()

  ggplot(rawFileInfo1, aes(x = Method, y = ions)) +
    geom_violin(aes(fill = Method, color = Method)) +
    theme_minimal() +
    scale_fill_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue)) +
    scale_y_log10() +
    geom_boxplot(color = "grey", fill = "white", width = 0.1, outlier.size = 0) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
    ylab("Total ions per spectrum (calculated)")+
    geom_hline(data = rawFileInfoSummary, aes(yintercept = Target, color = Method))+
    facet_wrap(~Sample, nrow = 2)+
    scale_color_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue))
}


plotIons2 <- function(rawFileInfo1, ...) {
  files <- list(...)
  for (f in files) rawFileInfo1 <- rbind(rawFileInfo1, f)

  rawFileInfo1 <- rawFileInfo1 %>%
    arrange(Analyzer, maxIT, desc(IsolationWindow)) %>%
    mutate(Method=factor(Method, levels = unique(Method)))

  rawFileInfoSummary <- select(rawFileInfo1, Target, Sample, Method) %>%
    unique()

  ggplot(rawFileInfo1, aes(x = Method, y = Ions2)) +
    geom_violin(aes(fill = Method, color = Method)) +
    theme_minimal() +
    scale_fill_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue)) +
    scale_y_log10() +
    geom_boxplot(color = "grey", fill = "white", width = 0.1, outlier.size = 0) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
    ylab("Total ions per spectrum")+
    geom_hline(data = rawFileInfoSummary, aes(yintercept = Target, color = Method)) +
    facet_wrap(~Sample, nrow = 2) +
    scale_color_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue))
}


plotInjectTime <- function(rawFileInfo1, ...) {
  files <- list(...)
  for (f in files) rawFileInfo1 <- rbind(rawFileInfo1, f)

  rawFileInfo1 <- rawFileInfo1 %>%
    arrange(Analyzer, maxIT, desc(IsolationWindow)) %>%
    mutate(Method=factor(Method, levels = unique(Method)))


  ggplot(rawFileInfo1, aes(x = Method, y = InjectTime, color = Method)) +
    geom_violin(aes(fill = Method)) +
    theme_minimal() +
    scale_fill_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue))+
    scale_color_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue)) +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
    ylab("Injection time (ms)")+
    facet_wrap(~Sample, nrow = 2)
}


plotCentroids <- function(rawFileInfo1, ...) {
  files <- list(...)
  for (f in files) rawFileInfo1 <- rbind(rawFileInfo1, f)

  rawFileInfo1 <- rawFileInfo1 %>%
    arrange(Analyzer, maxIT, desc(IsolationWindow)) %>%
    mutate(Method=factor(Method, levels = unique(Method)))


  ggplot(rawFileInfo1, aes(x = Method, y = Centroids, color = Method)) +
    geom_violin(aes(fill = Method)) +
    theme_minimal() +
    scale_fill_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue))+
    scale_color_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen2, cBlue))  +
    theme(legend.position = "none", axis.title.x  = element_blank(),
          axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
    ylab("Centroids per spectrum") +
    coord_cartesian(ylim = c(0, 10000))+
    facet_wrap(~Sample, nrow = 2)
}


plotSpeed <- function(rawFileInfo1, ...) {
  files <- list(...)
  for (f in files) rawFileInfo1 <- rbind(rawFileInfo1, f)

  rawFileInfo1 <- rawFileInfo1 %>%
    arrange(Analyzer, maxIT, desc(IsolationWindow)) %>%
    mutate(Method=factor(Method, levels = unique(Method)))


  ggplot(rawFileInfo1, aes(x = Method, y = SpeedHz)) +
    geom_boxplot(aes(color = Method), outlier.size = 0, outlier.alpha = 0) +
    theme_minimal() +
    scale_color_manual(values = c(cRed, cOrange1, cYellow, cTeal, cGreen1, cBlue))  +
    theme(legend.position = "none", axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 25, vjust = 1, hjust=1)) +
    ylab("Spectrum acquisition Speed (Hz)") +
    coord_cartesian(ylim = c(0, 200))+
    facet_wrap(~Sample, nrow = 2)

}


plotMassError <- function(massErrorFile, ionStatsFiles){
  masserror <- read.delim(massErrorFile,
                          sep = "\t") %>%
    separate_rows(Raw.Spectrum.Ids, Raw.Mass.Errors, Raw.Intensities, Raw.Times, sep = ",") %>%
    mutate(Raw.Mass.Errors = as.numeric(Raw.Mass.Errors))%>%
    mutate(Raw.Intensities = as.numeric(Raw.Intensities)) %>%
    mutate(Raw.Times = as.numeric(Raw.Times))%>%
    mutate(Raw.Spectrum.Ids = as.numeric(gsub("0.1.", "", Raw.Spectrum.Ids))) %>%
    filter(Raw.Times >= Start.Time & Raw.Times <= End.Time)


  ionStatsMassError <- read.csv(ionStatsFiles) %>%
    select(ï..scanNumber, InjectTime) %>%
    set_colnames(c("Raw.Spectrum.Ids", "InjectTime"))


  masserror <- masserror %>%
    left_join(ionStatsMassError, by = "Raw.Spectrum.Ids") %>%
    mutate(ions = Raw.Intensities * InjectTime /1000)

  ggplot(masserror, aes(x = ions, y = Raw.Mass.Errors)) +
    geom_bin2d(bins = 200) +
    theme_minimal() +
    scale_x_log10() +
    ylab("Mass error (ppm)") +
    xlab("Number of ions") +
    scale_fill_gradientn(colors = c("white", cBlue, cRed, cYellow), name = "Fragment ion peaks") +
    theme(legend.position = c(0.85, 0.25))


}


plotIonCountDist <- function(peakInfoFile, threshold = 5) {

  scanIonStats <- read.csv(peakInfoFile) %>%
    select(1:3) %>%
    set_colnames(c("Signal", "Scan", "InjectionTime")) %>%
    mutate(ions = Signal * InjectionTime/1000) %>%
    mutate(InjectionTimeGroup = as.factor(round(InjectionTime/2, 0)*2))


  fractionPeaks <- nrow(filter(scanIonStats, ions <= threshold))/nrow(scanIonStats)

  fractionTIC <- sum(filter(scanIonStats, ions <= threshold)$ions)/sum(scanIonStats$ions)

  print(paste0("Fraction of peaks less than ", threshold, " ions: ", fractionPeaks*100,
               "% \n Fraction of TIC less than", threshold, " ions: ", fractionTIC *100, "%"))


  ggplot(scanIonStats, aes(x = ions, color = InjectionTimeGroup)) +
    geom_density(size =0.5) +
    theme_minimal() +
    scale_x_log10() +
    scale_color_manual(values = c(cRed, cOrange3, cYellow, cGreen2, cTeal, cBlue),
                       labels = c("0-2", "2-4", "4-6", "6-8", "8-10", "10"),
                       name = "Injection time bin (ms)") +
    ylab("Density") +
    xlab("Number of ions") 


}


plotDIAScheme <- function(Width1, Time1, Min1, Max1, Label1 ,
                          Width2, Time2, Min2, Max2, Label2,
                          Width3, Time3, Min3, Max3, Label3,
                          Width4, Time4, Min4, Max4, Label4){

  length1 <- (Max1-Min1)/Width1 * 1
  windowInfo1 <- data.frame(xmin = seq(0, 10^6, by = Time1)[1:length1],
                            xmax = seq(Time1, 10^6, by = Time1)[1:length1],
                            ymin = rep(seq(Min1, Max1, by = Width1), 5)[1:length1],
                            ymax = rep(seq((Min1+Width1), (Max1+Width1), by = Width1), 5)[1:length1]) %>%
    mutate(WindowScheme = Label1)

  length2 <- (Max2-Min2)/Width2 * 1
  windowInfo2 <- data.frame(xmin = seq(0, 10^6, by = Time2)[1:length2],
                            xmax = seq(Time2, 10^6, by = Time2)[1:length2],
                            ymin = rep(seq(Min2, Max2, by = Width2), 6)[1:length2],
                            ymax = rep(seq((Min2+Width2), (Max2+Width2), by = Width2), 6)[1:length2]) %>%
    mutate(WindowScheme = Label2)

  length3 <- (Max3-Min3)/Width3 * 1
  windowInfo3 <- data.frame(xmin = seq(0, 10^6, by = Time3)[1:length3],
                            xmax = seq(Time3, 10^6, by = Time3)[1:length3],
                            ymin = rep(seq(Min3, Max3, by = Width3), 6)[1:length3],
                            ymax = rep(seq((Min3+Width3), (Max3+Width3), by = Width3), 6)[1:length3]) %>%
    mutate(WindowScheme = Label3)

  length4 <- (Max4-Min4)/Width4 * 1
  windowInfo4 <- data.frame(xmin = seq(0, 10^6, by = Time4)[1:length4],
                            xmax = seq(Time4, 10^6, by = Time4)[1:length4],
                            ymin = rep(seq(Min4, Max4, by = Width4), 6)[1:length4],
                            ymax = rep(seq((Min4+Width4), (Max4+Width4), by = Width4), 6)[1:length4]) %>%
    mutate(WindowScheme = Label4)

  windowInfo <- rbind(windowInfo4,
                      windowInfo2,
                      windowInfo3,
                      windowInfo1) %>%
    mutate(WindowScheme = factor(WindowScheme,
                                 levels =c(Label4, Label3, Label2, Label1)))

  ggplot(windowInfo, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax,
                         fill = WindowScheme)) +
    geom_rect(color = "black") +
    theme_minimal()  +
    scale_fill_manual(values = c( cTeal, cYellow, cOrange1, cRed),
                      guide = guide_legend(reverse = TRUE)) +
    labs(fill = "") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size=9),
          legend.text = element_text(size=8),
          legend.key.height = unit(0.35, 'cm'),
          legend.key.width = unit(0.35, 'cm')) +
    ylab("m/z") +
    xlab("Time") +
    coord_cartesian(ylim = c(Min1, Max1))+
    guides(fill=guide_legend(nrow=2, byrow=TRUE)) 


}


plotdDIABounds <- function(MS2Features, bounds) {

  allFeatures <- read.csv(MS2Features)

  bounds <- read.csv(bounds)


  ggplot(allFeatures, aes(x = Average.Measured.Retention.Time, y = Precursor.Mz))+
    geom_hex(binwidth = c(0.7, 20)) +
    theme_minimal() +
    geom_line(data = bounds, aes(x = RT, y = Lowmz), color = "#264653", size = 1) +
    geom_line(data = bounds, aes(x = RT, y = Highmz), color = "#264653", size = 1) +
    scale_fill_gradient2(low = "white", high = "#f94144", mid = "#f9c74f", midpoint = 65, name = "Peptides") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())  +
    ylab ("m/z") +
    xlab("Retention Time (min)") +
    scale_x_continuous(position = "top") +
    coord_cartesian(ylim = c(390, 968), xlim = c(7, 32))


}
