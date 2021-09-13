

rm(list = ls())

library(ggplot2)
library(patchwork)

######################
## Load the curves of the severe malaria signature performance in other infections 

# DF
load("./Objs/Dengue1_Curves.rda")
load("./Objs/Dengue2_Curves.rda")
load("./Objs/Dengue3_Curves.rda")
load("./Objs/Dengue4_Curves.rda")
load("./Objs/Dengue5_Curves.rda")
load("./Objs/Dengue6_Curves.rda")

## Many infections
load("./Objs/Adenovirus1_Curves.rda")
load("./Objs/ManyInfections1_Curves.rda")
load("./Objs/ManyInfections2_Curves.rda")
load("./Objs/ManyInfections3_Curves.rda")

## TB
load("./Objs/TB1_Curves.rda")
load("./Objs/TB2_Curves.rda")
load("./Objs/TB3_Curves.rda")
load("./Objs/TB4_Curves.rda")

# TB and HIV
load("./Objs/HIVandTB_Curves.rda")

# West Nile virus
load("./Objs/WestNile1_Curves.rda")

# Meningitis
load("./Objs/Meningitis_Curves_severe.rda")

#############################################
ROC_Dengue$labels$title <- "GSE51808 (DF vs healthy)"
ROC_Dengue$theme$plot.title$size <- 6

ROC_Dengue2$labels$title <- "GSE96656 (DF vs healthy)"
ROC_Dengue2$theme$plot.title$size <- 6

ROC_Dengue3$labels$title <- "GSE25001 (un-complicated vs complicated DF)"  
ROC_Dengue3$theme$plot.title$size <- 5

ROC_Dengue4$labels$title <- "GSE18090 (DF vs healthy)"
ROC_Dengue4$theme$plot.title$size <- 6

ROC_Dengue5$labels$title <- "GSE17924 (un-complicated vs complicated DF)"
ROC_Dengue5$theme$plot.title$size <- 5

ROC_Dengue6$labels$title <- "GSE13052 (un-complicated vs complicated DF)"
ROC_Dengue6$theme$plot.title$size <- 5


ROC_Adenovirus$labels$title <- "GSE40396 (Bacterial and viral infections vs healthy)"
ROC_Adenovirus$theme$plot.title$size <- 6
ROC_ManyInfections1$labels$title <- "GSE42026 (Bacterial and viral infections vs healthy)"
ROC_ManyInfections1$theme$plot.title$size <- 6
ROC_ManyInfections2$labels$title <- "GSE6269-GPL96 (Bacterial and viral infections vs healthy)"
ROC_ManyInfections2$theme$plot.title$size <- 6
ROC_ManyInfections3$labels$title <- "GSE63990 (Bacterial and viral infections vs healthy)"
ROC_ManyInfections3$theme$plot.title$size <- 6

ROC_TB1$labels$title <- "GSE19444"
ROC_TB2$labels$title <- "GSE73408"
ROC_TB3$labels$title <- "GSE62525"
ROC_TB4$labels$title <- "GSE83456"

ROC_HIVandTB$labels$title <- "GSE39940 (HIV/TB vs healthy)"
ROC_HIVandTB$theme$plot.title$size <- 6
ROC_WestNile$labels$title <- "GSE46681 (Asymptomatic vs severe West Nile viral infection)"
ROC_WestNile$theme$plot.title$size <- 6


ROC_Meningitis2_severe$labels$title <- "GSE80496"
ROC_Meningitis4_severe$labels$title <- "GSE40586"
ROC_Meningitis2_severe$theme$plot.title$size <- 9
ROC_Meningitis4_severe$theme$plot.title$size <- 9

##############################################
## Dengue Fever
tiff(filename = "./Figs/CompMalariaSigPerformance_Dengue.tiff", width = 2500, height = 2000, res = 350)
((ROC_Dengue / ROC_Dengue2 & theme(plot.tag = element_text(size = 8))) | (ROC_Dengue4 / ROC_Dengue3 & theme(plot.tag = element_text(size = 8))) | (ROC_Dengue5 / ROC_Dengue6 & theme(plot.tag = element_text(size = 8))
)) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the severe malaria signature in dengue fever datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

##############################################
## TB 
png(filename = "./Figs/CompMalariaSigPerformance_TB.png", width = 2500, height = 2000, res = 300)
(ROC_TB1 / ROC_TB2 | ROC_TB3 / ROC_TB4 & theme(plot.tag = element_text(size = 12))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the severe malaria signature in TB datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )
dev.off()


##############################################
## Many infections
png(filename = "./Figs/CompMalariaSigPerformance_OtherInfections.png", width = 2500, height = 2000, res = 300)
(ROC_Adenovirus / ROC_ManyInfections1 | ROC_ManyInfections2 / ROC_ManyInfections3 | ROC_HIVandTB / ROC_WestNile & theme(plot.tag = element_text(size = 12))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the severe malaria signature in other viral and bacterial infections datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

##############################################
## Meningitis
tiff(filename = "./Figs/new/SevereMalariaSigPerformance_Meningitis.tiff", width = 2500, height = 2000, res = 300)
(ROC_Meningitis2_severe | ROC_Meningitis4_severe & theme(plot.tag = element_text(size = 12))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the severe malaria signature in meningitis datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

