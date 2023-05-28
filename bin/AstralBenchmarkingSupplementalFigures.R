library(tidyverse, quietly = T)
library(magrittr)
library(ggpubr)
library(enrichR)

setwd("D:/data/Olympus/20230403_Oly_Data/AstralBenchmarking/AstralBenchmarking")


MMCC_AS_3p5msResults <- "data/20230517_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_Ev3_GPFLib"

MMCC_AS_10msResults <- "data/20230517_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_Ev3_GPFLib"

MMCC_AS_10msdDIAResults <- "data/20230517_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev3_GPFLib"

MMCC_OT_90minResults <- "data/20230517_Lumos_MMCC_1ug_90min_OT_22ms_8Th_Ev3_GPFLib"

MMCC_OT_24minResults <- "data/20230517_OLEP08_MMCC_1ug_MB_24min_OT_23ms_2Th_Ev3_GPFLib"



Plasma_30min_3p5msResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_30min_AS_3p5ms_2Th_GPFLib"

Plasma_30min_10msResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_30min_AS_10ms_4Th_GPFLib"

Plasma_30min_10ms_dDIAResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_30min_AS_dDIA_20ms_4Th_GPFLib"

Plasma_60min_15msResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_60min_AS_15ms_4Th_GPFLib"

Plasma_60min_15ms_dDIAResults <- "data/20230517_OLEP08_EV_TP_1ug_MB_60min_AS_dDIA_15ms_2Th_GPFLib"


QuantReports_All <- readQuantReports(Plasma_30min_3p5msResults, Plasma_30min_10msResults, Plasma_30min_10ms_dDIAResults,
                                  Plasma_60min_15msResults, Plasma_60min_15ms_dDIAResults)



rawInfo_24min_ASChrLib <- readRawInfo(fileName = "data/20230410_OLEP08_HeLa_1ug_MB_24min_AS_23ms_2Th_ChrLib_4_Info_CS.csv",
                                          Analyzer = "Astral",
                                          IsolationWindow = 2,
                                          SampleType = "",
                                          maxIT = 23,
                                          AGCtarget = 50000,
                                          MethodType = "ChrLib")


rawInfo_24min_OTChrLib <- readRawInfo(fileName = "data/20230410_OLEP08_HeLa_1ug_MB_24min_OT_23ms_2Th_ChrLib_4_Info_CS.csv",
                                          Analyzer = "Orbitrap",
                                          IsolationWindow = 2,
                                          SampleType = "",
                                          maxIT = 23,
                                          AGCtarget = 1000000,
                                          MethodType = "ChrLib")

chrLibIons <-
  plotIons2(rawInfo_24min_ASChrLib, rawInfo_24min_OTChrLib) +
  scale_color_manual(values = c(cRed, cTeal, cOrange1, cYellow, cTeal, cGreen1, cBlue))+
  scale_fill_manual(values = c(cRed, cTeal, cOrange1, cYellow, cTeal, cGreen1, cBlue))

chrLibIT <-
  plotInjectTime(rawInfo_24min_ASChrLib, rawInfo_24min_OTChrLib) +
  scale_fill_manual(values = c(cRed, cTeal, cOrange1, cYellow, cTeal, cGreen1, cBlue))+
  scale_color_manual(values = c(cRed, cTeal, cOrange1, cYellow, cTeal, cGreen1, cBlue))

chrLibCentroids <-
  plotCentroids(rawInfo_24min_ASChrLib, rawInfo_24min_OTChrLib) +
  scale_fill_manual(values = c(cRed, cTeal, cOrange1, cYellow, cTeal, cGreen1, cBlue))+
  scale_color_manual(values = c(cRed, cTeal, cOrange1, cYellow, cTeal, cGreen1, cBlue))

chrLibSpeed <-
  plotSpeed(rawInfo_24min_ASChrLib, rawInfo_24min_OTChrLib) +
  scale_color_manual(values = c(cRed, cTeal, cOrange1, cYellow, cTeal, cGreen1, cBlue))

ggarrange(chrLibIons, chrLibIT,
          chrLibCentroids, chrLibSpeed,
          labels = c("A", "B", "C", "D"))

ggsave("results/AstralBenchmarking_FigureS1.jpg", height = 8, width = 9, dpi = 700)



rawInfo_24min_AS_3p5ms_2Th <- readRawInfo(fileName = "data/20230406_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_I_1_Info_CS.csv",
                                          Analyzer = "Astral",
                                          IsolationWindow = 2,
                                          SampleType = "",
                                          maxIT = 3.5,
                                          AGCtarget = 50000,
                                          MethodType = "DIA")


rawInfo_24min_AS_10ms_4Th <- readRawInfo(fileName = "data/20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1_Info_CS.csv",
                                         Analyzer = "Astral",
                                         IsolationWindow = 4,
                                         SampleType = "",
                                         maxIT = 10,
                                         AGCtarget = 50000,
                                         MethodType = "DIA")

rawInfo_24min_AS_10ms_2Th <- readRawInfo(fileName = "data/20230406_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_I_1_Info_CS.csv",
                                         Analyzer = "Astral",
                                         IsolationWindow = 2,
                                         SampleType = "",
                                         maxIT = 10,
                                         AGCtarget = 50000,
                                         MethodType = "dDIA")

rawInfo_24min_OT_23ms_24Th <- readRawInfo(fileName = "data/20230406_OLEP08_MMCC_1ug_MB_24min_OT_23ms_24Th_I_1_Info_CS.csv",
                                          Analyzer = "Orbitrap",
                                          IsolationWindow = 24,
                                          SampleType = "",
                                          maxIT = 23,
                                          AGCtarget = 1000000,
                                          MethodType = "DIA")


rawInfo_24min_OT_23ms_8Th <- readRawInfo(fileName = "data/20210906_LRH_MMCC_Static_K_2_Info_CS.csv",
                                         Analyzer = "Orbitrap",
                                         IsolationWindow = 8,
                                         SampleType = "",
                                         maxIT = 23,
                                         AGCtarget = 1000000,
                                         MethodType = "DIA")

rawInfo_24min_OT_23ms_2Th <- readRawInfo(fileName = "data/20230406_OLEP08_MMCC_1ug_MB_24min_OT_75Th_23ms_2Th_I_1_Info_CS.csv",
                                         Analyzer = "Orbitrap",
                                         IsolationWindow = 2,
                                         SampleType = "",
                                         maxIT = 23,
                                         AGCtarget = 1000000,
                                         MethodType = "DIA*")

quantIons <-
  plotIons2(rawInfo_24min_AS_3p5ms_2Th, rawInfo_24min_AS_10ms_4Th,
            rawInfo_24min_AS_10ms_2Th, rawInfo_24min_OT_23ms_8Th, rawInfo_24min_OT_23ms_2Th)

quantIT <-
  plotInjectTime(rawInfo_24min_AS_3p5ms_2Th, rawInfo_24min_AS_10ms_4Th,
                 rawInfo_24min_AS_10ms_2Th, rawInfo_24min_OT_23ms_8Th, rawInfo_24min_OT_23ms_2Th)

quantCentroids <-
  plotCentroids(rawInfo_24min_AS_3p5ms_2Th, rawInfo_24min_AS_10ms_4Th,
                rawInfo_24min_AS_10ms_2Th, rawInfo_24min_OT_23ms_8Th, rawInfo_24min_OT_23ms_2Th)

quantSpeed <-
  plotSpeed(rawInfo_24min_AS_3p5ms_2Th, rawInfo_24min_AS_10ms_4Th,
            rawInfo_24min_AS_10ms_2Th, rawInfo_24min_OT_23ms_8Th, rawInfo_24min_OT_23ms_2Th)


ggarrange(quantIons, quantIT,
          quantCentroids, quantSpeed,
          labels = c("A", "B", "C", "D"))

ggsave("results/AstralBenchmarking_FigureS2.jpg", height = 8, width = 9.5, dpi = 700)



quantMMCC1 <- readQuantReportsMMCC(MMCC_AS_3p5msResults, MMCC_AS_10msResults,
                                   MMCC_AS_10msdDIAResults, MMCC_OT_90minResults)

quantMMCCSelected <- quantMMCC1 %>%
  filter(SampleLabel == "100%" | SampleLabel == "50%" | SampleLabel == "10%") %>%
  mutate(SampleLabel = paste0(10*Analyte.Concentration, " ng"))


peptideUpset <- upsetPlot(quantMMCCSelected, "Peptide")

peptideUpset

ggsave("results/AstralBenchmarking_FigureS3.jpg", height = 6, width = 9, dpi = 700)




quantMMCC_Shortmassrange <- readQuantReportsMMCC(MMCC_AS_3p5msResults, MMCC_AS_10msResults,
                                                 MMCC_AS_10msdDIAResults, MMCC_OT_24minResults) %>%
  filter(Precursor.Mz <= 576 & Precursor.Mz >= 500.5)


peptideUpsetShort <- upsetPlot(quantMMCC_Shortmassrange, "Peptide")

peptideUpsetShort

ggsave("results/AstralBenchmarking_FigureS4.jpg", height = 6, width = 9, dpi = 700)



quantRatioDensityProtein(quantMMCC1, 100, 10)

ggsave("results/AstralBenchmarking_FigureS7.jpg", height = 9, width = 9, dpi = 700)



LOQhisto_short <- plotLOQMMCC(quantMMCC_Shortmassrange)

LOQsummary_short <- plotLOQMMCCSummary(quantMMCC_Shortmassrange)

pairwiseComp_short <- plotLOQMMCCHistogram(quantMMCC_Shortmassrange, MMCC_AS_10msResults,
                                           "Astral",
                                           MMCC_OT_24minResults,
                                           "Orbitrap")

quantRatioPlot_short <- quantRatioDensity(quantMMCC_Shortmassrange, 100, 10)


ggarrange(ggarrange(LOQsummary_short, pairwiseComp_short, nrow = 2,
                    labels = c("A", "B"),
                    heights = c(1, 0.75)),LOQhisto_short, quantRatioPlot_short, ncol = 3,
          widths = c(0.8, 0.8, 0.7), labels = c("", "C", "D")
)

ggsave("results/AstralBenchmarking_FigureS8.jpg", height = 7, width = 14, dpi = 700)


quantRatioDensity(quantMMCC1, 100, 1, Alpha = 0.1, multiplier = 250)

ggsave("results/AstralBenchmarking_FigureS9.jpg", height = 9, width = 9, dpi = 700)


LOQ_Histo_dDIA_3p5 <- plotLOQMMCCHistogram(quantMMCC1, MMCC_AS_10msdDIAResults, 
                                           "Astral, 2 Th, 10 ms, dDIA",
                                           MMCC_AS_3p5msResults,
                                           "Astral, 2 Th, 3.5 ms, DIA") +
  scale_fill_manual(values = c(cYellow, "grey", cRed),
                    name = "Better LLOQ in: ")




LOQ_Histo_dDIA_10 <- plotLOQMMCCHistogram(quantMMCC1, MMCC_AS_10msdDIAResults, 
                     "Astral, 2 Th, 10 ms, dDIA",
                     MMCC_AS_10msResults,
                     "Astral, 4 Th, 10 ms, DIA") +
  scale_fill_manual(values = c(cYellow, "grey", cOrange1),
                    name = "Better LLOQ in: ")


ggarrange(LOQ_Histo_dDIA_3p5, LOQ_Histo_dDIA_10, nrow = 1, labels = c("A", "B"))

ggsave("results/AstralBenchmarking_FigureS11.jpg", height = 5, width = 9, dpi = 700)







EnrichedList <- c("CD9", "CD63", "FLOT1" ,
                  "NCAM1",  "PDC6I")

DepletedList <- c("ALBU" ,"TRFE",
                  "APOA1",  "APOB",  "HBA")

plotRankOrder(calcFoldChange(Plasma_60min_15msResults,
                             QuantReports_All, EnrichedList, DepletedList),
              labels = TRUE)


ggsave("results/AstralBenchmarking_FigureS15.jpg", height = 5, width = 9, dpi = 700)



plasmaCV_All <- plotCVs(QuantReports_All, Transitions = "All")

plasmaCV_Refined <- plotCVs(QuantReports_All, Transitions = "Refined")

ggarrange(plasmaCV_All, plasmaCV_Refined, labels = c("A", "B"),
          nrow = 2)

ggsave("results/AstralBenchmarking_FigureSK.jpg", height = 9, width = 9, dpi = 700)




quantMMCC_Shortmassrange2 <- quantMMCC_Shortmassrange %>%
  filter(SampleLabel == "100%" | SampleLabel == "50%" | SampleLabel == "10%") %>%
  mutate(SampleLabel = paste0(10*Analyte.Concentration, " ng"))

proteinCVShort <- plotCVsMMCCProteins(quantMMCC_Shortmassrange2, "All transitions")

peptideCVShort <-
  plotCVsMMCCAll(quantMMCC_Shortmassrange2, "All transitions")

ggarrange(peptideCVShort, proteinCVShort, nrow = 1, labels = c("A", "B"))

ggsave("results/AstralBenchmarking_FigureS5.jpg", height = 7, width = 14, dpi = 700)



plotIonCountDist("data/20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_3_PeakInfo_120000to130000.csv")

ggsave("results/AstralBenchmarking_FigureS14.jpg", height = 7, width = 14, dpi = 700)




quantRatioDensity(quantMMCC_Shortmassrange, 100, 1, Alpha = 0.8, multiplier = 250)

ggsave("results/AstralBenchmarking_FigureS10.jpg", height = 9, width = 9, dpi = 700)




plotMassError("data/20230403_OLEP08_EV_1ug_MB_30min_AS_10ms_4Th_1_MassErrorChromatograms.tsv", 
              "data/20230403_OLEP08_EV_1ug_MB_30min_AS_10ms_4Th_1_Info_CS.csv")


ggsave("results/AstralBenchmarking_FigureS12.jpg", height = 6, width = 8, dpi = 700)


ResultsSummary <- writeLOQMMCCSummary(quantMMCC1)

write.table(writeLOQMMCCSummary(quantMMCC1), "results/SupplementalTable2.tsv",
            row.names = FALSE, quote = FALSE, sep = "\t")

resultsShort <- writeLOQMMCCSummary(quantMMCC_Shortmassrange)



plotLOQMMCCSummary(quantMMCC1,
                   OTmultiplier = 90,
                   ASmultiplier = 24,
                   AxisLabel = "Peptides per hour")

ggsave("results/AstralBenchmarking_FigureS6.jpg", height = 6, width = 8, dpi = 700)

