
rm(list = ls())

library(ggplot2)
library(patchwork)

######################
## Load the curves of the cerebral malaria signature performance in other infections 

# DF
load("./Objs/Dengue1_Curves_cerebral.rda")
load("./Objs/Dengue2_Curves_cerebral.rda")
load("./Objs/Dengue3_Curves_cerebral.rda")
load("./Objs/Dengue4_Curves_cerebral.rda")
load("./Objs/Dengue5_Curves_cerebral.rda")
load("./Objs/Dengue6_Curves_cerebral.rda")

## Many infections
load("./Objs/Adenovirus1_Curves_cerebral.rda")
load("./Objs/ManyInfections1_Curves_cerebral.rda")
load("./Objs/ManyInfections2_Curves_cerebral.rda")
load("./Objs/ManyInfections3_Curves_cerebral.rda")

## TB
load("./Objs/TB1_Curves_cerebral.rda")
load("./Objs/TB2_Curves_cerebral.rda")
load("./Objs/TB3_Curves_cerebral.rda")
load("./Objs/TB4_Curves_cerebral.rda")

# TB and HIV
load("./Objs/HIVandTB_Curves_cerebral.rda")

# West Nile virus
load("./Objs/WestNile1_Curves_cerebral.rda")

#############################################
ROC_Dengue_cerebral$labels$title <- "GSE51808 (DF vs healthy)"
ROC_Dengue_cerebral$theme$plot.title$size <- 6

ROC_Dengue2_cerebral$labels$title <- "GSE96656 (DF vs healthy)"
ROC_Dengue2_cerebral$theme$plot.title$size <- 6

ROC_Dengue3_cerebral$labels$title <- "GSE25001 (un-complicated vs complicated DF)"  
ROC_Dengue3_cerebral$theme$plot.title$size <- 5

ROC_Dengue4_cerebral$labels$title <- "GSE18090 (DF vs healthy)"
ROC_Dengue4_cerebral$theme$plot.title$size <- 6

ROC_Dengue5_cerebral$labels$title <- "GSE17924 (un-complicated vs complicated DF)"
ROC_Dengue5_cerebral$theme$plot.title$size <- 5

ROC_Dengue6_cerebral$labels$title <- "GSE13052 (un-complicated vs complicated DF)"
ROC_Dengue6_cerebral$theme$plot.title$size <- 5


ROC_Adenovirus_cerebral$labels$title <- "GSE40396 (Bacterial and viral infections vs healthy)"
ROC_Adenovirus_cerebral$theme$plot.title$size <- 6
ROC_ManyInfections1_cerebral$labels$title <- "GSE42026 (Bacterial and viral infections vs healthy)"
ROC_ManyInfections1_cerebral$theme$plot.title$size <- 6
ROC_ManyInfections2_cerebral$labels$title <- "GSE6269-GPL96 (Bacterial and viral infections vs healthy)"
ROC_ManyInfections2_cerebral$theme$plot.title$size <- 6
ROC_ManyInfections3_cerebral$labels$title <- "GSE63990 (Bacterial and viral infections vs healthy)"
ROC_ManyInfections3_cerebral$theme$plot.title$size <- 6

ROC_TB1_cerebral$labels$title <- "GSE19444"
ROC_TB2_cerebral$labels$title <- "GSE73408"
ROC_TB3_cerebral$labels$title <- "GSE62525"
ROC_TB4_cerebral$labels$title <- "GSE83456"

ROC_HIVandTB_cerebral$labels$title <- "GSE39940 (HIV/TB vs healthy)"
ROC_HIVandTB_cerebral$theme$plot.title$size <- 6
ROC_WestNile_cerebral$labels$title <- "GSE46681 (Asymptomatic vs severe West Nile viral infection)"
ROC_WestNile_cerebral$theme$plot.title$size <- 6

##############################################
## Dengue Fever
tiff(filename = "./Figs/CerebralMalariaSigPerformance_Dengue.tiff", width = 2500, height = 2000, res = 350)
((ROC_Dengue_cerebral / ROC_Dengue2_cerebral & theme(plot.tag = element_text(size = 8))) | (ROC_Dengue4_cerebral / ROC_Dengue3_cerebral & theme(plot.tag = element_text(size = 8))) | (ROC_Dengue5_cerebral / ROC_Dengue6_cerebral & theme(plot.tag = element_text(size = 8))
)) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the cerebral malaria signatures in dengue fever datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

##############################################
## TB 
png(filename = "./Figs/CerebralMalariaSigPerformance_TB.png", width = 2500, height = 2000, res = 300)
(ROC_TB1_cerebral / ROC_TB2_cerebral | ROC_TB3_cerebral / ROC_TB4_cerebral & theme(plot.tag = element_text(size = 12))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the cerebral malaria signatures in TB datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 15, face = "bold"))
  )
dev.off()


##############################################
## Many infections
png(filename = "./Figs/CerebralMalariaSigPerformance_OtherInfections.png", width = 2500, height = 2000, res = 300)
(ROC_Adenovirus_cerebral / ROC_ManyInfections1_cerebral | ROC_ManyInfections2_cerebral / ROC_ManyInfections3_cerebral | ROC_HIVandTB_cerebral / ROC_WestNile_cerebral & theme(plot.tag = element_text(size = 12))
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The performance of the cerebral malaria signatures in other viral and bacterial infections datasets',
    tag_levels = c('A', '1'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

##############################################
