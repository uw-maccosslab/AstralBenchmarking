library(tidyverse, quietly = T)
library(magrittr)
library(ggpubr)
library(enrichR)

setwd("D:/data/Olympus/20230403_Oly_Data/AstralBenchmarking/AstralBenchmarking")



MMCC_AS_3p5msResults <- "data/20230517_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_Ev3_GPFLib"

MMCC_AS_10msResults <- "data/20230517_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_Ev3_GPFLib"

MMCC_AS_10msdDIAResults <- "data/20230517_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev3_GPFLib"

MMCC_OT_90minResults <- "data/20230517_Lumos_MMCC_1ug_90min_OT_22ms_8Th_Ev3_GPFLib"



Plasma_30min_3p5msResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_30min_AS_3p5ms_2Th_GPFLib"

Plasma_30min_10msResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_30min_AS_10ms_4Th_GPFLib"

Plasma_30min_10ms_dDIAResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_30min_AS_dDIA_20ms_4Th_GPFLib"

Plasma_60min_15msResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_60min_AS_15ms_4Th_GPFLib"

Plasma_60min_15ms_dDIAResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_60min_AS_dDIA_15ms_2Th_GPFLib"





#Fig 1

quantMMCC1 <- readQuantReportsMMCC(MMCC_AS_3p5msResults, MMCC_AS_10msResults,
                                   MMCC_AS_10msdDIAResults, MMCC_OT_90minResults)

quantMMCCSelected <- quantMMCC1 %>%
  filter(SampleLabel == "100%" | SampleLabel == "50%" | SampleLabel == "10%") %>%
  mutate(SampleLabel = paste0(10*Analyte.Concentration, " ng"))


DIAbounds <- plotdDIABounds("data/RT_DIA_Scheduling_HeLa.csv", 
                            "data/20230406scheduled_dia_10ms_2Th_75.0Hz_2.00sec_fit.csv")


proteinUpset <- upsetPlot(quantMMCCSelected, "Protein")


exampleScheme <- plotDIAScheme(2, 5, 600, 648, "182 Hz, 2 Th",
                               4, 11.5, 600, 648, "87 Hz, 4 Th",
                               2, 11.5, 612, 632, "87 Hz, 2 Th",
                               8, 37, 600, 648, "27 Hz, 8 Th")



blank <- grid::grid.rect(gp=grid::gpar(col="white"))




ggarrange(blank, ggarrange(exampleScheme, DIAbounds, nrow =1, widths = c(1,1.35),
                    labels = c("B", "C")), proteinUpset, nrow = 3,
          labels = c("A", "", "D"),
          heights = c(0.65, 0.5, 0.7))

ggsave("results/AstralBenchmarking_Figure1.jpg", height = 8, width = 6.66, dpi = 700)


# Fig 2


proteinCV <- plotCVsMMCCProteins(quantMMCCSelected, "All transitions")

peptideCV <-
  plotCVsMMCCAll(quantMMCCSelected, "All transitions")

ggarrange(peptideCV, proteinCV, nrow = 1, labels = c("A", "B"))

ggsave("results/AstralBenchmarking_Figure2.jpg", height = 4, width = 14, dpi = 700)



# Fig 3

LOQhisto <- plotLOQMMCC(quantMMCC1)

LOQsummary <- plotLOQMMCCSummary(quantMMCC1)

pairwiseComp <- plotLOQMMCCHistogram(quantMMCC1, MMCC_AS_3p5msResults,
                     "Astral, 24 minute gradient",
                     MMCC_OT_90minResults,
                     "Orbitrap, 90 minute gradient")

quantRatioPlot <- quantRatioDensity(quantMMCC1, 100, 10)



ggarrange(ggarrange(LOQsummary, pairwiseComp, nrow = 2,
                    labels = c("A", "B"),
                    heights = c(1, 0.75)),LOQhisto, quantRatioPlot, ncol = 3,
          widths = c(0.7, 0.8, 0.8), labels = c("", "C", "D")
)

ggsave("results/AstralBenchmarking_Figure3.jpg", height = 7, width = 14, dpi = 700)


# Figure 4


EnrichedList <- c("CD9", "CD63", "FLOT1" ,
                  "NCAM1",  "PDC6I")

DepletedList <- c("ALBU" ,"TRFE",
                  "APOA1",  "APOB",  "HBA")

QuantReports_All <- readQuantReports(Plasma_30min_3p5msResults, Plasma_30min_10msResults,
                                     Plasma_30min_10ms_dDIAResults, Plasma_60min_15msResults,
                                     Plasma_60min_15ms_dDIAResults)



plotCV <- plotCVs((QuantReports_All), Transitions = "Refined")+
  theme(axis.text.x = element_text(size = 8),
        plot.margin = unit(c(0,1,0,1), 'lines'))


pepSummary <- summarizeSearch(QuantReports_All, Level = "Peptide") +
  theme(axis.text.x = element_blank(),
        plot.margin = unit(c(0,1,0,1), 'lines'))

protSummary <- summarizeSearch(QuantReports_All, Level = "Protein") +
  theme(axis.text.x = element_blank(),
        plot.margin = unit(c(0,1,0,1), 'lines'))


volcano <- plotFoldChange(calcFoldChange(Plasma_60min_15msResults,
                                 QuantReports_All, EnrichedList, DepletedList),
               labels = TRUE)


ggarrange(pepSummary, protSummary, plotCV, volcano,
          labels = c("A", "", "B", "C"),
          ncol = 1, nrow = 4, heights = c(0.65, 0.65, 1.2, 1))

ggsave("results/AstralBenchmarking_Figure4.jpg", height = 8, width = 6.66, dpi = 700)
