library(tidyverse, quietly = T)
library(magrittr)
library(ggpubr)
library(enrichR)

setwd("D:/data/Olympus/20230403_Oly_Data/bin")



rawInfo_24min_ASChrLib <- readRawInfo(fileName = "../MMCC/20230410_OLEP08_HeLa_1ug_MB_24min_AS_23ms_2Th_ChrLib_4_Info_CS.csv",
                                          Analyzer = "Astral",
                                          IsolationWindow = 2,
                                          SampleType = "",
                                          maxIT = 23,
                                          AGCtarget = 50000,
                                          MethodType = "ChrLib")


rawInfo_24min_OTChrLib <- readRawInfo(fileName = "../MMCC/20230410_OLEP08_HeLa_1ug_MB_24min_OT_23ms_2Th_ChrLib_4_Info_CS.csv",
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

ggsave("SupplementalFigures/HeLaChrLibIonStats.jpg", height = 8, width = 9, dpi = 700)




peptideUpset <- upsetPlot(quantMMCCSelected, "Peptide")

peptideUpset

ggsave("SupplementalFigures/HeLaPeptideUpset.jpg", height = 6, width = 9, dpi = 700)



rawInfo_24min_AS_3p5ms_2Th <- readRawInfo(fileName = "../MMCC/20230406_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_I_1_Info_CS.csv",
                                      Analyzer = "Astral",
                                      IsolationWindow = 2,
                                      SampleType = "",
                                      maxIT = 3.5,
                                      AGCtarget = 50000,
                                      MethodType = "DIA")


rawInfo_24min_AS_10ms_4Th <- readRawInfo(fileName = "../MMCC/20230406_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_I_1_Info_CS.csv",
                                      Analyzer = "Astral",
                                      IsolationWindow = 4,
                                      SampleType = "",
                                      maxIT = 10,
                                      AGCtarget = 50000,
                                      MethodType = "DIA")

rawInfo_24min_AS_10ms_2Th <- readRawInfo(fileName = "../MMCC/20230406_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_I_1_Info_CS.csv",
                                         Analyzer = "Astral",
                                         IsolationWindow = 2,
                                         SampleType = "",
                                         maxIT = 10,
                                         AGCtarget = 50000,
                                         MethodType = "dDIA")

rawInfo_24min_OT_23ms_24Th <- readRawInfo(fileName = "../MMCC/20230406_OLEP08_MMCC_1ug_MB_24min_OT_23ms_24Th_I_1_Info_CS.csv",
                                         Analyzer = "Orbitrap",
                                         IsolationWindow = 24,
                                         SampleType = "",
                                         maxIT = 23,
                                         AGCtarget = 1000000,
                                         MethodType = "DIA")


rawInfo_24min_OT_23ms_8Th <- readRawInfo(fileName = "../MMCC/20210906_LRH_MMCC_Static_K_2_Info_CS.csv",
                                          Analyzer = "Orbitrap",
                                          IsolationWindow = 8,
                                          SampleType = "",
                                          maxIT = 23,
                                          AGCtarget = 1000000,
                                          MethodType = "DIA")

rawInfo_24min_OT_23ms_2Th <- readRawInfo(fileName = "../MMCC/20230406_OLEP08_MMCC_1ug_MB_24min_OT_75Th_23ms_2Th_I_1_Info_CS.csv",
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

ggsave("SupplementalFigures/HeLaQuantIonStats.jpg", height = 8, width = 9.5, dpi = 700)





quantMMCC_Shortmassrange <- readQuantReportsMMCC("../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_Ev1p12_GPFLib",
                                                 "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_Ev1p12_GPFLib",
                                                 "../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_OT_23ms_2Th_Ev2_GPFLib",
                                                 "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev1p12_GPFLib") %>%
  filter(Precursor.Mz <= 576 & Precursor.Mz >= 500.5)

quantRatioDensity(quantMMCC_Shortmassrange, 100, 1, Alpha = 0.8, multiplier = 250)

ggsave("SupplementalFigures/20230501_Ratio_shortMassRange_100to1.jpg", height = 9, width = 9, dpi = 700)



peptideUpsetShort <- upsetPlot(quantMMCC_Shortmassrange, "Peptide")

peptideUpsetShort

ggsave("SupplementalFigures/HeLaPeptideUpsetShortMassRange.jpg", height = 6, width = 9, dpi = 700)



quantMMCC_Shortmassrange <- quantMMCC_Shortmassrange %>%
  filter(SampleLabel == "100%" | SampleLabel == "50%" | SampleLabel == "10%") %>%
  mutate(SampleLabel = paste0(10*Analyte.Concentration, " ng"))

proteinCVShort <- plotCVsMMCCProteins(quantMMCC_Shortmassrange, "All transitions")

peptideCVShort <-
  plotCVsMMCCAll(quantMMCC_Shortmassrange, "All transitions")

ggarrange(peptideCVShort, proteinCVShort, nrow = 1, labels = c("A", "B"))

ggsave("SupplementalFigures/20230501_CVSummary_Short_v2.jpg", height = 7, width = 14, dpi = 700)




#writeDiffactoInputMMCC("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_Ev1_GPFLib", quantMMCC1)

#writeDiffactoInputMMCC("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_Ev1_GPFLib", quantMMCC1)

#writeDiffactoInputMMCC("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_OT_23ms_24Th_Ev2_GPFLib", quantMMCC1)

#writeDiffactoInputMMCC("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev1_GPFLib", quantMMCC1)


#plotDiffactoRatio(diffactoRatioReadIn("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_Ev1_GPFLib",
#                                      "Astral, 4 Th, 10 ms, DIA"),
#              diffactoRatioReadIn("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_Ev1_GPFLib",
  #                                    "Astral, 2 Th, 3.5 ms, DIA"),
 #                 diffactoRatioReadIn("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev1_GPFLib",
 #                                     "Astral, 2 Th, 10 ms, dDIA"),
 #             diffactoRatioReadIn("../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_OT_23ms_24Th_Ev2_GPFLib",
  #                                "Orbitrap, 24 Th, 23 ms, DIA"))


quantRatioDensityProtein(quantMMCC1, 100, 10)


ggsave("SupplementalFigures/HeLaProteinRawQuant_10percent.jpg", height = 9, width = 9, dpi = 700)


LOQ_Histo_dDIA_3p5 <- plotLOQMMCCHistogram(quantMMCC1, "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_Ev1p12_GPFLib",
                     "2 Th, 3.5 ms, DIA",
                     "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev1p12_GPFLib",
                     "2 Th, 10 ms, dDIA")





LOQ_Histo_dDIA_10 <- plotLOQMMCCHistogram(quantMMCC1, "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_Ev1p12_GPFLib",
                     "4 Th, 10 ms, DIA",
                     "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev1p12_GPFLib",
                     "2 Th, 10 ms, dDIA")
ggarrange(LOQ_Histo_dDIA_3p5, LOQ_Histo_dDIA_10, nrow = 1, labels = c("A", "B"))

ggsave("SupplementalFigures/HeLaLOQComp_dDIAvsDIA.jpg", height = 5, width = 9, dpi = 700)



plotLOQRatioScatter(quantMMCC1, "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_3p5ms_2Th_Ev1p12_GPFLib",
                     "2 Th, 3.5 ms, DIA",
                     "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_dDIA_10ms_2Th_Ev1p12_GPFLib",
                     "2 Th, 10 ms, dDIA") +
  coord_cartesian(xlim = c(0, 5e5))


plotRankOrder(calcFoldChange("../Plasma_60min/20230404_OLEP08_EV_TP_1ug_MB_60min_AS_dDIA_15ms_2Th_GPFLib",
                             QuantReports_All, EnrichedList, DepletedList),
              labels = TRUE)


ggsave("SupplementalFigures/PlasmaRankOrderPlot_v1.jpg", height = 5, width = 9, dpi = 700)



DepletedList2 <- c()

EnrichedList2 <- c("CRP", "PLMN", "PCSK9" ,
                  "PLF4",  "CXCL7", "CFAD",  "INHBC",
                   "APOE" ,
                  "APOA1",  "APOB", "IL18", "PARK7", "A4", "SYUA")

plotRankOrder2(calcFoldChange("../Plasma_60min/20230404_OLEP08_EV_TP_1ug_MB_60min_AS_dDIA_15ms_2Th_GPFLib",
                              QuantReports_All, EnrichedList2, DepletedList2),
               labels = TRUE)

ggsave("SupplementalFigures/PlasmaRankOrderPlot_v2.jpg", height = 4.5, width = 8, dpi = 700)

plotRankOrder2(calcFoldChange("../Plasma_60min/20230404_OLEP08_EV_TP_1ug_MB_60min_AS_15ms_4Th_GPFLib",
                              QuantReports_All, EnrichedList2, DepletedList2),
               labels = TRUE)

ggsave("SupplementalFigures/PlasmaRankOrderPlot_v3.jpg", height = 5, width = 8, dpi = 700)




LOQhisto_short <- plotLOQMMCC(quantMMCC_Shortmassrange)

LOQsummary_short <- plotLOQMMCCSummary(quantMMCC_Shortmassrange)

pairwiseComp_short <- plotLOQMMCCHistogram(quantMMCC_Shortmassrange, "../MMCC/20230428_OLEP08_MMCC_1ug_MB_24min_AS_10ms_4Th_Ev1p12_GPFLib",
                                     "Astral",
                                     "../MMCC/20230425_OLEP08_MMCC_1ug_MB_24min_OT_23ms_2Th_Ev2_GPFLib",
                                     "Orbitrap")

quantRatioPlot_short <- quantRatioDensity(quantMMCC_Shortmassrange, 100, 10)


ggarrange(ggarrange(LOQsummary_short, pairwiseComp_short, nrow = 2,
                    labels = c("A", "B"),
                    heights = c(1, 0.75)),LOQhisto_short, quantRatioPlot_short, ncol = 3,
          widths = c(0.8, 0.8, 0.7), labels = c("", "C", "D")
)

ggsave("SupplementalFigures/20230501_LOQSummary_ShortMassRange.jpg", height = 7, width = 14, dpi = 700)



quantRatioDensity(quantMMCC1, 100, 1, Alpha = 0.1, multiplier = 250)

ggsave("SupplementalFigures/20230501_Ratio_100to1.jpg", height = 9, width = 9, dpi = 700)





 plotCVs(QuantReports_All)

 ggsave("SupplementalFigures/20230504_PlasmaCV.jpg", height = 7, width = 9, dpi = 700)
