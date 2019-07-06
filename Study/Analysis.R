#Real Data Analysis
#Load Data
setwd("C:/Users/jeff/Desktop/Study/Data/Final")
#Breast Cancer
BreastSelectedA_PI3K_Akt_CN = read.csv("Breast_GPL96_PI3K-Akt_CN.csv", header = TRUE)[ ,-1]
BreastSelectedA_Wnt_CN = read.csv("Breast_GPL96_Wnt_CN.csv", header = TRUE)[ ,-1]
BreastSelectedA_PI3K_Akt_unique_CN = read.csv("Breast_GPL96_PI3K-Akt_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
BreastSelectedA_Wnt_unique_CN = read.csv("Breast_GPL96_Wnt_uniqueprobe_CN.csv", header = TRUE)[ ,-1]

BreastSelectedA_PI3K_Akt_NC = read.csv("Breast_GPL96_PI3K-Akt_NC.csv", header = TRUE)[ ,-1]
BreastSelectedA_Wnt_NC = read.csv("Breast_GPL96_Wnt_NC.csv")[ ,-1]
BreastSelectedA_PI3K_Akt_unique_NC = read.csv("Breast_GPL96_PI3K-Akt_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
BreastSelectedA_Wnt_unique_NC = read.csv("Breast_GPL96_Wnt_uniqueprobe_NC.csv", header = TRUE)[ ,-1]


#Colorectal Cancer
ColorectalSelectedB_p53_CN	 = read.csv("Colorectal_GPL570_p53_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_PI3K_Akt_CN	 = read.csv("Colorectal_GPL570_PI3K-Akt_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Ras_CN	 = read.csv("Colorectal_GPL570_Ras_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_TGF_beta_CN	 = read.csv("Colorectal_GPL570_TGF-beta_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Wnt_CN	 = read.csv("Colorectal_GPL570_Wnt_CN.csv", header = TRUE)[ ,-1]

ColorectalSelectedB_p53_unique_CN	 = read.csv("Colorectal_GPL570_p53_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_PI3K_Akt_unique_CN	 = read.csv("Colorectal_GPL570_PI3K-Akt_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Ras_unique_CN	 = read.csv("Colorectal_GPL570_Ras_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_TGF_beta_unique_CN	 = read.csv("Colorectal_GPL570_TGF-beta_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Wnt_unique_CN	 = read.csv("Colorectal_GPL570_Wnt_uniqueprobe_CN.csv", header = TRUE)[ ,-1]


ColorectalSelectedB_p53_NC	 = read.csv("Colorectal_GPL570_p53_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_PI3K_Akt_NC	 = read.csv("Colorectal_GPL570_PI3K-Akt_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Ras_NC	 = read.csv("Colorectal_GPL570_Ras_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_TGF_beta_NC	 = read.csv("Colorectal_GPL570_TGF-beta_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Wnt_NC	 = read.csv("Colorectal_GPL570_Wnt_NC.csv", header = TRUE)[ ,-1]

ColorectalSelectedB_p53_unique_NC	 = read.csv("Colorectal_GPL570_p53_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_PI3K_Akt_unique_NC	 = read.csv("Colorectal_GPL570_PI3K-Akt_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Ras_unique_NC	 = read.csv("Colorectal_GPL570_Ras_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_TGF_beta_unique_NC	 = read.csv("Colorectal_GPL570_TGF-beta_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
ColorectalSelectedB_Wnt_unique_NC	 = read.csv("Colorectal_GPL570_Wnt_uniqueprobe_NC.csv", header = TRUE)[ ,-1]


#Leukemia	
LeukemiaSelectedA_Jak_STAT_CN	 = read.csv("Leukemia_GPL96_Jak-STAT_CN.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_mTOR_CN	 = read.csv("Leukemia_GPL96_mTOR_CN.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_PI3K_Akt_CN	 = read.csv("Leukemia_GPL96_PI3K-Akt_CN.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_Wnt_CN	 = read.csv("Leukemia_GPL96_Wnt_CN.csv", header = TRUE)[ ,-1]

LeukemiaSelectedA_Jak_STAT_unique_CN	 = read.csv("Leukemia_GPL96_Jak-STAT_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_mTOR_unique_CN	 = read.csv("Leukemia_GPL96_mTOR_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_PI3K_Akt_unique_CN	 = read.csv("Leukemia_GPL96_PI3K-Akt_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_Wnt_unique_CN	 = read.csv("Leukemia_GPL96_Wnt_uniqueprobe_CN.csv", header = TRUE)[ ,-1]

LeukemiaSelectedA_Jak_STAT_NC	 = read.csv("Leukemia_GPL96_Jak-STAT_NC.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_mTOR_NC	 = read.csv("Leukemia_GPL96_mTOR_NC.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_PI3K_Akt_NC	 = read.csv("Leukemia_GPL96_PI3K-Akt_NC.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_Wnt_NC	 = read.csv("Leukemia_GPL96_Wnt_NC.csv", header = TRUE)[ ,-1]

LeukemiaSelectedA_Jak_STAT_unique_NC	 = read.csv("Leukemia_GPL96_Jak-STAT_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_mTOR_unique_NC	 = read.csv("Leukemia_GPL96_mTOR_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_PI3K_Akt_unique_NC	 = read.csv("Leukemia_GPL96_PI3K-Akt_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
LeukemiaSelectedA_Wnt_unique_NC	 = read.csv("Leukemia_GPL96_Wnt_uniqueprobe_NC.csv", header = TRUE)[ ,-1]


#Lung Cancer(GPL96)	
LungSelectedA_Cell_cycle_CN 	 = read.csv("Lung_GPL96_Cell cycle_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_EGFR_CN 	 = read.csv("Lung_GPL96_EGFR_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_Jak_STAT_CN 	 = read.csv("Lung_GPL96_Jak-STAT_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_mapk_CN 	 = read.csv("Lung_GPL96_mapk_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_mTOR_CN 	 = read.csv("Lung_GPL96_mTOR_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_p53_CN 	 = read.csv("Lung_GPL96_p53_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_Wnt_CN 	 = read.csv("Lung_GPL96_Wnt_CN.csv", header = TRUE)[ ,-1] 

LungSelectedA_Cell_cycle_unique_CN 	 = read.csv("Lung_GPL96_Cell cycle_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_EGFR_unique_CN 	 = read.csv("Lung_GPL96_EGFR_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_Jak_STAT_unique_CN 	 = read.csv("Lung_GPL96_Jak-STAT_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_mapk_unique_CN 	 = read.csv("Lung_GPL96_mapk_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_mTOR_unique_CN 	 = read.csv("Lung_GPL96_mTOR_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_p53_unique_CN 	 = read.csv("Lung_GPL96_p53_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedA_Wnt_unique_CN 	 = read.csv("Lung_GPL96_Wnt_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 

LungSelectedA_Cell_cycle_NC 	 = read.csv("Lung_GPL96_Cell cycle_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_EGFR_NC 	 = read.csv("Lung_GPL96_EGFR_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_Jak_STAT_NC 	 = read.csv("Lung_GPL96_Jak-STAT_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_mapk_NC 	 = read.csv("Lung_GPL96_mapk_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_mTOR_NC 	 = read.csv("Lung_GPL96_mTOR_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_p53_NC 	 = read.csv("Lung_GPL96_p53_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_Wnt_NC 	 = read.csv("Lung_GPL96_Wnt_NC.csv", header = TRUE)[ ,-1] 

LungSelectedA_Cell_cycle_unique_NC 	 = read.csv("Lung_GPL96_Cell cycle_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_EGFR_unique_NC 	 = read.csv("Lung_GPL96_EGFR_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_Jak_STAT_unique_NC 	 = read.csv("Lung_GPL96_Jak-STAT_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_mapk_unique_NC 	 = read.csv("Lung_GPL96_mapk_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_mTOR_unique_NC 	 = read.csv("Lung_GPL96_mTOR_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_p53_unique_NC 	 = read.csv("Lung_GPL96_p53_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedA_Wnt_unique_NC 	 = read.csv("Lung_GPL96_Wnt_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 


#Lung Cancer(GPL570)	
LungSelectedB_Cell_cycle_CN 	 = read.csv("Lung_GPL570_Cell cycle_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_EGFR_CN 	 = read.csv("Lung_GPL570_EGFR_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_Jak_STAT_CN 	 = read.csv("Lung_GPL570_Jak-STAT_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_mapk_CN 	 = read.csv("Lung_GPL570_mapk_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_mTOR_CN 	 = read.csv("Lung_GPL570_mTOR_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_p53_CN 	 = read.csv("Lung_GPL570_p53_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_Wnt_CN 	 = read.csv("Lung_GPL570_Wnt_CN.csv", header = TRUE)[ ,-1] 

LungSelectedB_Cell_cycle_unique_CN 	 = read.csv("Lung_GPL570_Cell cycle_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_EGFR_unique_CN 	 = read.csv("Lung_GPL570_EGFR_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_Jak_STAT_unique_CN 	 = read.csv("Lung_GPL570_Jak-STAT_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_mapk_unique_CN 	 = read.csv("Lung_GPL570_mapk_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_mTOR_unique_CN 	 = read.csv("Lung_GPL570_mTOR_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_p53_unique_CN 	 = read.csv("Lung_GPL570_p53_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 
LungSelectedB_Wnt_unique_CN 	 = read.csv("Lung_GPL570_Wnt_uniqueprobe_CN.csv", header = TRUE)[ ,-1] 

LungSelectedB_Cell_cycle_NC 	 = read.csv("Lung_GPL570_Cell cycle_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_EGFR_NC 	 = read.csv("Lung_GPL570_EGFR_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_Jak_STAT_NC 	 = read.csv("Lung_GPL570_Jak-STAT_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_mapk_NC 	 = read.csv("Lung_GPL570_mapk_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_mTOR_NC 	 = read.csv("Lung_GPL570_mTOR_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_p53_NC 	 = read.csv("Lung_GPL570_p53_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_Wnt_NC 	 = read.csv("Lung_GPL570_Wnt_NC.csv", header = TRUE)[ ,-1] 

LungSelectedB_Cell_cycle_unique_NC 	 = read.csv("Lung_GPL570_Cell cycle_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_EGFR_unique_NC 	 = read.csv("Lung_GPL570_EGFR_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_Jak_STAT_unique_NC 	 = read.csv("Lung_GPL570_Jak-STAT_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_mapk_unique_NC 	 = read.csv("Lung_GPL570_mapk_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_mTOR_unique_NC 	 = read.csv("Lung_GPL570_mTOR_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_p53_unique_NC 	 = read.csv("Lung_GPL570_p53_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 
LungSelectedB_Wnt_unique_NC 	 = read.csv("Lung_GPL570_Wnt_uniqueprobe_NC.csv", header = TRUE)[ ,-1] 


#Ovarian Cancer	
OvarianSelectedA_mTOR_mean_unique_CN	 = read.csv("Ovarian_GPL96_mTOR_mean_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
OvarianSelectedA_Jak_STAT_mean_unique_CN	 = read.csv("Ovarian_GPL96_Jak-STAT_mean_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
OvarianSelectedA_p53_mean_unique_CN	 = read.csv("Ovarian_GPL96_p53_mean_uniqueprobe_CN.csv", header = TRUE)[ ,-1]
OvarianSelectedA_PI3K_Akt_mean_unique_CN	 = read.csv("Ovarian_GPL96_PI3K-Akt_mean_uniqueprobe_CN.csv", header = TRUE)[ ,-1]

OvarianSelectedA_mTOR_mean_unique_NC	 = read.csv("Ovarian_GPL96_mTOR_mean_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
OvarianSelectedA_Jak_STAT_mean_unique_NC	 = read.csv("Ovarian_GPL96_Jak-STAT_mean_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
OvarianSelectedA_p53_mean_unique_NC	 = read.csv("Ovarian_GPL96_p53_mean_uniqueprobe_NC.csv", header = TRUE)[ ,-1]
OvarianSelectedA_PI3K_Akt_mean_unique_NC	 = read.csv("Ovarian_GPL96_PI3K-Akt_mean_uniqueprobe_NC.csv", header = TRUE)[ ,-1]


###############################################################################################################
#Classified by different Diseases
Breast = list()
#Cancer Normal
Breast[["BreastA_PI3K_Akt_CN"]] = t(BreastSelectedA_PI3K_Akt_CN)
Breast[["BreastA_Wnt_CN"]] = t(BreastSelectedA_Wnt_CN)
#unique probe
Breast[["BreastA_PI3K_Akt_unique_CN"]] = t(BreastSelectedA_PI3K_Akt_unique_CN)
Breast[["BreastA_Wnt_unique_CN"]] = t(BreastSelectedA_Wnt_unique_CN)

#Normal Cancer
Breast[["BreastA_PI3K_Akt_NC"]] = t(BreastSelectedA_PI3K_Akt_NC)
Breast[["BreastA_Wnt_NC"]] = t(BreastSelectedA_Wnt_NC)
#unique probe
Breast[["BreastA_PI3K_Akt_unique_NC"]] = t(BreastSelectedA_PI3K_Akt_unique_NC)
Breast[["BreastA_Wnt_unique_NC"]] = t(BreastSelectedA_Wnt_unique_NC)



Colorectal = list()
#Cancer Normal
Colorectal[["ColorectalB_p53_CN"]] = t(ColorectalSelectedB_p53_CN)
Colorectal[["ColorectalB_PI3K_Akt_CN"]] = t(ColorectalSelectedB_PI3K_Akt_CN)
Colorectal[["ColorectalB_Ras_CN"]] = t(ColorectalSelectedB_Ras_CN)
Colorectal[["ColorectalB_TGF_beta_CN"]] = t(ColorectalSelectedB_TGF_beta_CN)
Colorectal[["ColorectalB_Wnt_CN"]] = t(ColorectalSelectedB_Wnt_CN)
#unique probe
Colorectal[["ColorectalB_p53_unique_CN"]] = t(ColorectalSelectedB_p53_unique_CN)
Colorectal[["ColorectalB_PI3K_Akt_unique_CN"]] = t(ColorectalSelectedB_PI3K_Akt_unique_CN)
Colorectal[["ColorectalB_Ras_unique_CN"]] = t(ColorectalSelectedB_Ras_unique_CN)
Colorectal[["ColorectalB_TGF_beta_unique_CN"]] = t(ColorectalSelectedB_TGF_beta_unique_CN)
Colorectal[["ColorectalB_Wnt_unique_CN"]] = t(ColorectalSelectedB_Wnt_unique_CN)

#Normal Cancer
Colorectal[["ColorectalB_p53_NC"]] = t(ColorectalSelectedB_p53_NC)
Colorectal[["ColorectalB_PI3K_Akt_NC"]] = t(ColorectalSelectedB_PI3K_Akt_NC)
Colorectal[["ColorectalB_Ras_NC"]] = t(ColorectalSelectedB_Ras_NC)
Colorectal[["ColorectalB_TGF_beta_NC"]] = t(ColorectalSelectedB_TGF_beta_NC)
Colorectal[["ColorectalB_Wnt_NC"]] = t(ColorectalSelectedB_Wnt_NC)
#unique probe
Colorectal[["ColorectalB_p53_unique_NC"]] = t(ColorectalSelectedB_p53_unique_NC)
Colorectal[["ColorectalB_PI3K_Akt_unique_NC"]] = t(ColorectalSelectedB_PI3K_Akt_unique_NC)
Colorectal[["ColorectalB_Ras_unique_NC"]] = t(ColorectalSelectedB_Ras_unique_NC)
Colorectal[["ColorectalB_TGF_beta_unique_NC"]] = t(ColorectalSelectedB_TGF_beta_unique_NC)
Colorectal[["ColorectalB_Wnt_unique_NC"]] = t(ColorectalSelectedB_Wnt_unique_NC)



Leukemia_CN = list()
#Cancer Normal
Leukemia_CN[["LeukemiaA_Jak_STAT_CN"]] = t(LeukemiaSelectedA_Jak_STAT_CN)
Leukemia_CN[["LeukemiaA_mTOR_CN"]] = t(LeukemiaSelectedA_mTOR_CN)
Leukemia_CN[["LeukemiaA_PI3K_Akt_CN"]] = t(LeukemiaSelectedA_PI3K_Akt_CN)
Leukemia_CN[["LeukemiaA_Wnt_CN"]] = t(LeukemiaSelectedA_Wnt_CN)
#unique probe
Leukemia_CN[["LeukemiaA_Jak_STAT_unique_CN"]] = t(LeukemiaSelectedA_Jak_STAT_unique_CN)
Leukemia_CN[["LeukemiaA_mTOR_unique_CN"]] = t(LeukemiaSelectedA_mTOR_unique_CN)
Leukemia_CN[["LeukemiaA_PI3K_Akt_unique_CN"]] = t(LeukemiaSelectedA_PI3K_Akt_unique_CN)
Leukemia_CN[["LeukemiaA_Wnt_unique_CN"]] = t(LeukemiaSelectedA_Wnt_unique_CN)

Leukemia_NC = list()
#Normal Cancer
Leukemia_NC[["LeukemiaA_Jak_STAT_NC"]] = t(LeukemiaSelectedA_Jak_STAT_NC)
Leukemia_NC[["LeukemiaA_mTOR_NC"]] = t(LeukemiaSelectedA_mTOR_NC)
Leukemia_NC[["LeukemiaA_PI3K_Akt_NC"]] = t(LeukemiaSelectedA_PI3K_Akt_NC)
Leukemia_NC[["LeukemiaA_Wnt_NC"]] = t(LeukemiaSelectedA_Wnt_NC)
#unique probe
Leukemia_NC[["LeukemiaA_Jak_STAT_unique_NC"]] = t(LeukemiaSelectedA_Jak_STAT_unique_NC)
Leukemia_NC[["LeukemiaA_mTOR_unique_NC"]] = t(LeukemiaSelectedA_mTOR_unique_NC)
Leukemia_NC[["LeukemiaA_PI3K_Akt_unique_NC"]] = t(LeukemiaSelectedA_PI3K_Akt_unique_NC)
Leukemia_NC[["LeukemiaA_Wnt_unique_NC"]] = t(LeukemiaSelectedA_Wnt_unique_NC)



LungA = list()
#Cancer Normal
LungA[["LungA_Cell_cycle_CN"]] = t(LungSelectedA_Cell_cycle_CN) 
LungA[["LungA_EGFR_CN"]] = t(LungSelectedA_EGFR_CN) 
LungA[["LungA_Jak_STAT_CN"]] = t(LungSelectedA_Jak_STAT_CN) 
LungA[["LungA_mapk_CN"]] = t(LungSelectedA_mapk_CN) 
LungA[["LungA_mTOR_CN"]] = t(LungSelectedA_mTOR_CN) 
LungA[["LungA_p53_CN"]] = t(LungSelectedA_p53_CN) 
LungA[["LungA_Wnt_CN"]] = t(LungSelectedA_Wnt_CN) 
#unique probe
LungA[["LungA_Cell_cycle_unique_CN"]] = t(LungSelectedA_Cell_cycle_unique_CN) 
LungA[["LungA_EGFR_unique_CN"]] = t(LungSelectedA_EGFR_unique_CN) 
LungA[["LungA_Jak_STAT_unique_CN"]] = t(LungSelectedA_Jak_STAT_unique_CN) 
LungA[["LungA_mapk_unique_CN"]] = t(LungSelectedA_mapk_unique_CN) 
LungA[["LungA_mTOR_unique_CN"]] = t(LungSelectedA_mTOR_unique_CN) 
LungA[["LungA_p53_unique_CN"]] = t(LungSelectedA_p53_unique_CN) 
LungA[["LungA_Wnt_unique_CN"]] = t(LungSelectedA_Wnt_unique_CN) 

#Normal Cancer
LungA[["LungA_Cell_cycle_NC"]] = t(LungSelectedA_Cell_cycle_NC) 
LungA[["LungA_EGFR_NC"]] = t(LungSelectedA_EGFR_NC) 
LungA[["LungA_Jak_STAT_NC"]] = t(LungSelectedA_Jak_STAT_NC) 
LungA[["LungA_mapk_NC"]] = t(LungSelectedA_mapk_NC) 
LungA[["LungA_mTOR_NC"]] = t(LungSelectedA_mTOR_NC) 
LungA[["LungA_p53_NC"]] = t(LungSelectedA_p53_NC) 
LungA[["LungA_Wnt_NC"]] = t(LungSelectedA_Wnt_NC) 
#unique probe
LungA[["LungA_Cell_cycle_unique_NC"]] = t(LungSelectedA_Cell_cycle_unique_NC) 
LungA[["LungA_EGFR_unique_NC"]] = t(LungSelectedA_EGFR_unique_NC) 
LungA[["LungA_Jak_STAT_unique_NC"]] = t(LungSelectedA_Jak_STAT_unique_NC) 
LungA[["LungA_mapk_unique_NC"]] = t(LungSelectedA_mapk_unique_NC) 
LungA[["LungA_mTOR_unique_NC"]] = t(LungSelectedA_mTOR_unique_NC) 
LungA[["LungA_p53_unique_NC"]] = t(LungSelectedA_p53_unique_NC) 
LungA[["LungA_Wnt_unique_NC"]] = t(LungSelectedA_Wnt_unique_NC) 



LungB = list()
#Cancer Normal
LungB[["LungB_Cell_cycle_CN"]] = t(LungSelectedB_Cell_cycle_CN) 
LungB[["LungB_EGFR_CN"]] = t(LungSelectedB_EGFR_CN) 
LungB[["LungB_Jak_STAT_CN"]] = t(LungSelectedB_Jak_STAT_CN) 
LungB[["LungB_mapk_CN"]] = t(LungSelectedB_mapk_CN) 
LungB[["LungB_mTOR_CN"]] = t(LungSelectedB_mTOR_CN) 
LungB[["LungB_p53_CN"]] = t(LungSelectedB_p53_CN) 
LungB[["LungB_Wnt_CN"]] = t(LungSelectedB_Wnt_CN) 
#unique probe
LungB[["LungB_Cell_cycle_unique_CN"]] = t(LungSelectedB_Cell_cycle_unique_CN) 
LungB[["LungB_EGFR_unique_CN"]] = t(LungSelectedB_EGFR_unique_CN) 
LungB[["LungB_Jak_STAT_unique_CN"]] = t(LungSelectedB_Jak_STAT_unique_CN) 
LungB[["LungB_mapk_unique_CN"]] = t(LungSelectedB_mapk_unique_CN) 
LungB[["LungB_mTOR_unique_CN"]] = t(LungSelectedB_mTOR_unique_CN) 
LungB[["LungB_p53_unique_CN"]] = t(LungSelectedB_p53_unique_CN) 
LungB[["LungB_Wnt_unique_CN"]] = t(LungSelectedB_Wnt_unique_CN) 

#Normal Cancer
LungB[["LungB_Cell_cycle_NC"]] = t(LungSelectedB_Cell_cycle_NC) 
LungB[["LungB_EGFR_NC"]] = t(LungSelectedB_EGFR_NC) 
LungB[["LungB_Jak_STAT_NC"]] = t(LungSelectedB_Jak_STAT_NC) 
LungB[["LungB_mapk_NC"]] = t(LungSelectedB_mapk_NC) 
LungB[["LungB_mTOR_NC"]] = t(LungSelectedB_mTOR_NC) 
LungB[["LungB_p53_NC"]] = t(LungSelectedB_p53_NC) 
LungB[["LungB_Wnt_NC"]] = t(LungSelectedB_Wnt_NC) 
#unique probe
LungB[["LungB_Cell_cycle_unique_NC"]] = t(LungSelectedB_Cell_cycle_unique_NC) 
LungB[["LungB_EGFR_unique_NC"]] = t(LungSelectedB_EGFR_unique_NC) 
LungB[["LungB_Jak_STAT_unique_NC"]] = t(LungSelectedB_Jak_STAT_unique_NC) 
LungB[["LungB_mapk_unique_NC"]] = t(LungSelectedB_mapk_unique_NC) 
LungB[["LungB_mTOR_unique_NC"]] = t(LungSelectedB_mTOR_unique_NC) 
LungB[["LungB_p53_unique_NC"]] = t(LungSelectedB_p53_unique_NC) 
LungB[["LungB_Wnt_unique_NC"]] = t(LungSelectedB_Wnt_unique_NC) 



Ovarian_CN = list()
#Cancer Normal, unique mean probe
Ovarian_CN[["OvarianA_mTOR_mean_unique_CN"]] = t(OvarianSelectedA_mTOR_mean_unique_CN)
Ovarian_CN[["OvarianA_Jak_STAT_mean_unique_CN"]] = t(OvarianSelectedA_Jak_STAT_mean_unique_CN)
Ovarian_CN[["OvarianA_p53_mean_unique_CN"]] = t(OvarianSelectedA_p53_mean_unique_CN)
Ovarian_CN[["OvarianA_PI3K_Akt_mean_unique_CN"]] = t(OvarianSelectedA_PI3K_Akt_mean_unique_CN)

Ovarian_NC = list()
#Normal Cancer, unique mean probe
Ovarian_NC[["OvarianA_mTOR_mean_unique_NC"]] = t(OvarianSelectedA_mTOR_mean_unique_NC)
Ovarian_NC[["OvarianA_Jak_STAT_mean_unique_NC"]] = t(OvarianSelectedA_Jak_STAT_mean_unique_NC)
Ovarian_NC[["OvarianA_p53_mean_unique_NC"]] = t(OvarianSelectedA_p53_mean_unique_NC)
Ovarian_NC[["OvarianA_PI3K_Akt_mean_unique_NC"]] = t(OvarianSelectedA_PI3K_Akt_mean_unique_NC)



################################################################################################################################################
#Classified by unique probe
NonUniqueProbe = list()
#Cancer Normal
NonUniqueProbe[["BreastA_PI3K_Akt_CN"]] = t(BreastSelectedA_PI3K_Akt_CN)
NonUniqueProbe[["BreastA_Wnt_CN"]] = t(BreastSelectedA_Wnt_CN)
NonUniqueProbe[["ColorectalB_p53_CN"]] = t(ColorectalSelectedB_p53_CN)
NonUniqueProbe[["ColorectalB_PI3K_Akt_CN"]] = t(ColorectalSelectedB_PI3K_Akt_CN)
NonUniqueProbe[["ColorectalB_Ras_CN"]] = t(ColorectalSelectedB_Ras_CN)
NonUniqueProbe[["ColorectalB_TGF_beta_CN"]] = t(ColorectalSelectedB_TGF_beta_CN)
NonUniqueProbe[["ColorectalB_Wnt_CN"]] = t(ColorectalSelectedB_Wnt_CN)
NonUniqueProbe[["LungA_Cell_cycle_CN"]] = t(LungSelectedA_Cell_cycle_CN) 
NonUniqueProbe[["LungA_EGFR_CN"]] = t(LungSelectedA_EGFR_CN) 
NonUniqueProbe[["LungA_Jak_STAT_CN"]] = t(LungSelectedA_Jak_STAT_CN) 
NonUniqueProbe[["LungA_mapk_CN"]] = t(LungSelectedA_mapk_CN) 
NonUniqueProbe[["LungA_mTOR_CN"]] = t(LungSelectedA_mTOR_CN) 
NonUniqueProbe[["LungA_p53_CN"]] = t(LungSelectedA_p53_CN) 
NonUniqueProbe[["LungA_Wnt_CN"]] = t(LungSelectedA_Wnt_CN) 
NonUniqueProbe[["LungB_Cell_cycle_CN"]] = t(LungSelectedB_Cell_cycle_CN) 
NonUniqueProbe[["LungB_EGFR_CN"]] = t(LungSelectedB_EGFR_CN) 
NonUniqueProbe[["LungB_Jak_STAT_CN"]] = t(LungSelectedB_Jak_STAT_CN) 
NonUniqueProbe[["LungB_mapk_CN"]] = t(LungSelectedB_mapk_CN) 
NonUniqueProbe[["LungB_mTOR_CN"]] = t(LungSelectedB_mTOR_CN) 
NonUniqueProbe[["LungB_p53_CN"]] = t(LungSelectedB_p53_CN) 
NonUniqueProbe[["LungB_Wnt_CN"]] = t(LungSelectedB_Wnt_CN) 
#=====================================================================================
NonUniqueProbe[["BreastA_PI3K_Akt_NC"]] = t(BreastSelectedA_PI3K_Akt_NC)
NonUniqueProbe[["BreastA_Wnt_NC"]] = t(BreastSelectedA_Wnt_NC)
NonUniqueProbe[["ColorectalB_p53_NC"]] = t(ColorectalSelectedB_p53_NC)
NonUniqueProbe[["ColorectalB_PI3K_Akt_NC"]] = t(ColorectalSelectedB_PI3K_Akt_NC)
NonUniqueProbe[["ColorectalB_Ras_NC"]] = t(ColorectalSelectedB_Ras_NC)
NonUniqueProbe[["ColorectalB_TGF_beta_NC"]] = t(ColorectalSelectedB_TGF_beta_NC)
NonUniqueProbe[["ColorectalB_Wnt_NC"]] = t(ColorectalSelectedB_Wnt_NC)
NonUniqueProbe[["LungA_Cell_cycle_NC"]] = t(LungSelectedA_Cell_cycle_NC) 
NonUniqueProbe[["LungA_EGFR_NC"]] = t(LungSelectedA_EGFR_NC) 
NonUniqueProbe[["LungA_Jak_STAT_NC"]] = t(LungSelectedA_Jak_STAT_NC) 
NonUniqueProbe[["LungA_mapk_NC"]] = t(LungSelectedA_mapk_NC) 
NonUniqueProbe[["LungA_mTOR_NC"]] = t(LungSelectedA_mTOR_NC) 
NonUniqueProbe[["LungA_p53_NC"]] = t(LungSelectedA_p53_NC) 
NonUniqueProbe[["LungA_Wnt_NC"]] = t(LungSelectedA_Wnt_NC) 
NonUniqueProbe[["LungB_Cell_cycle_NC"]] = t(LungSelectedB_Cell_cycle_NC) 
NonUniqueProbe[["LungB_EGFR_NC"]] = t(LungSelectedB_EGFR_NC) 
NonUniqueProbe[["LungB_Jak_STAT_NC"]] = t(LungSelectedB_Jak_STAT_NC) 
NonUniqueProbe[["LungB_mapk_NC"]] = t(LungSelectedB_mapk_NC) 
NonUniqueProbe[["LungB_mTOR_NC"]] = t(LungSelectedB_mTOR_NC) 
NonUniqueProbe[["LungB_p53_NC"]] = t(LungSelectedB_p53_NC) 
NonUniqueProbe[["LungB_Wnt_NC"]] = t(LungSelectedB_Wnt_NC)

#***************************************************************************************
UniqueProbe = list()
UniqueProbe[["BreastA_PI3K_Akt_unique_CN"]] = t(BreastSelectedA_PI3K_Akt_unique_CN)
UniqueProbe[["BreastA_Wnt_unique_CN"]] = t(BreastSelectedA_Wnt_unique_CN)
UniqueProbe[["ColorectalB_p53_unique_CN"]] = t(ColorectalSelectedB_p53_unique_CN)
UniqueProbe[["ColorectalB_PI3K_Akt_unique_CN"]] = t(ColorectalSelectedB_PI3K_Akt_unique_CN)
UniqueProbe[["ColorectalB_Ras_unique_CN"]] = t(ColorectalSelectedB_Ras_unique_CN)
UniqueProbe[["ColorectalB_TGF_beta_unique_CN"]] = t(ColorectalSelectedB_TGF_beta_unique_CN)
UniqueProbe[["ColorectalB_Wnt_unique_CN"]] = t(ColorectalSelectedB_Wnt_unique_CN)
UniqueProbe[["LungA_Cell_cycle_unique_CN"]] = t(LungSelectedA_Cell_cycle_unique_CN) 
UniqueProbe[["LungA_EGFR_unique_CN"]] = t(LungSelectedA_EGFR_unique_CN) 
UniqueProbe[["LungA_Jak_STAT_unique_CN"]] = t(LungSelectedA_Jak_STAT_unique_CN) 
UniqueProbe[["LungA_mapk_unique_CN"]] = t(LungSelectedA_mapk_unique_CN) 
UniqueProbe[["LungA_mTOR_unique_CN"]] = t(LungSelectedA_mTOR_unique_CN) 
UniqueProbe[["LungA_p53_unique_CN"]] = t(LungSelectedA_p53_unique_CN) 
UniqueProbe[["LungA_Wnt_unique_CN"]] = t(LungSelectedA_Wnt_unique_CN) 
UniqueProbe[["LungB_Cell_cycle_unique_CN"]] = t(LungSelectedB_Cell_cycle_unique_CN) 
UniqueProbe[["LungB_EGFR_unique_CN"]] = t(LungSelectedB_EGFR_unique_CN) 
UniqueProbe[["LungB_Jak_STAT_unique_CN"]] = t(LungSelectedB_Jak_STAT_unique_CN) 
UniqueProbe[["LungB_mapk_unique_CN"]] = t(LungSelectedB_mapk_unique_CN) 
UniqueProbe[["LungB_mTOR_unique_CN"]] = t(LungSelectedB_mTOR_unique_CN) 
UniqueProbe[["LungB_p53_unique_CN"]] = t(LungSelectedB_p53_unique_CN) 
UniqueProbe[["LungB_Wnt_unique_CN"]] = t(LungSelectedB_Wnt_unique_CN) 
#=====================================================================================
UniqueProbe[["BreastA_PI3K_Akt_unique_NC"]] = t(BreastSelectedA_PI3K_Akt_unique_NC)
UniqueProbe[["BreastA_Wnt_unique_NC"]] = t(BreastSelectedA_Wnt_unique_NC)
UniqueProbe[["ColorectalB_p53_unique_NC"]] = t(ColorectalSelectedB_p53_unique_NC)
UniqueProbe[["ColorectalB_PI3K_Akt_unique_NC"]] = t(ColorectalSelectedB_PI3K_Akt_unique_NC)
UniqueProbe[["ColorectalB_Ras_unique_NC"]] = t(ColorectalSelectedB_Ras_unique_NC)
UniqueProbe[["ColorectalB_TGF_beta_unique_NC"]] = t(ColorectalSelectedB_TGF_beta_unique_NC)
UniqueProbe[["ColorectalB_Wnt_unique_NC"]] = t(ColorectalSelectedB_Wnt_unique_NC)
UniqueProbe[["LungA_Cell_cycle_unique_NC"]] = t(LungSelectedA_Cell_cycle_unique_NC) 
UniqueProbe[["LungA_EGFR_unique_NC"]] = t(LungSelectedA_EGFR_unique_NC) 
UniqueProbe[["LungA_Jak_STAT_unique_NC"]] = t(LungSelectedA_Jak_STAT_unique_NC) 
UniqueProbe[["LungA_mapk_unique_NC"]] = t(LungSelectedA_mapk_unique_NC) 
UniqueProbe[["LungA_mTOR_unique_NC"]] = t(LungSelectedA_mTOR_unique_NC) 
UniqueProbe[["LungA_p53_unique_NC"]] = t(LungSelectedA_p53_unique_NC) 
UniqueProbe[["LungA_Wnt_unique_NC"]] = t(LungSelectedA_Wnt_unique_NC) 
UniqueProbe[["LungB_Cell_cycle_unique_NC"]] = t(LungSelectedB_Cell_cycle_unique_NC) 
UniqueProbe[["LungB_EGFR_unique_NC"]] = t(LungSelectedB_EGFR_unique_NC) 
UniqueProbe[["LungB_Jak_STAT_unique_NC"]] = t(LungSelectedB_Jak_STAT_unique_NC) 
UniqueProbe[["LungB_mapk_unique_NC"]] = t(LungSelectedB_mapk_unique_NC) 
UniqueProbe[["LungB_mTOR_unique_NC"]] = t(LungSelectedB_mTOR_unique_NC) 
UniqueProbe[["LungB_p53_unique_NC"]] = t(LungSelectedB_p53_unique_NC) 
UniqueProbe[["LungB_Wnt_unique_NC"]] = t(LungSelectedB_Wnt_unique_NC) 


################################################################################################################################################
#Normality test
nonnormal = function(AAA, genenumber = 4, times = 100, adjacent = TRUE, disease = "Leukemia", log2 = FALSE){
  library(energy)
  library(MVN)
  sig <<- 0.05
  names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
  number = function(X){
    for (i in 1:length(names)) assign(names[i], 0)
    
    if(disease == "Leukemia"){
      selected = read.csv("C:/Users/jeff/Desktop/Study/Data/Leukemia/Platform_HG-U133A(GPL96)/Probes and samples/GSE9476_selected.csv", header = FALSE)
      tumor = "tumor"
      normal = "normal"
    }else if(disease == "Ovarian"){
      selected = read.csv("C:/Users/jeff/Desktop/Study/Data/Ovarian cancer/Platform_HG-U133A(GPL96)/Probes and samples/tcga_ov_preselected.csv", header = FALSE)
      tumor = 1
      normal = 0
    }
    
    ntum = sum(selected[,1] == tumor)
    nnor = sum(selected[,1] == normal)
    if(adjacent == TRUE){
      AA <<- X[1:(nrow(X)/2),]  #cancer 
      BB <<- X[((nrow(X)/2)+1):nrow(X),]   #normal
    }else if(adjacent == "CN"){
      AA <<- X[1:ntum,]
      BB <<- X[ntum + 1: ntum + nnor]
    }else if(adjacent == "NC"){
      AA <<- X[1:nnor,]
      BB <<- X[nnor + 1: ntum + nnor]
    }else{
      AA <<- c()
      BB <<- c()
    }
    
    if(log2){
      AA = log2(AA)
    }

    for(j in 1:times){
      if(j == 11){
        j = 1001
      }
      print(j)
      set.seed(j)
      
      probe = sample(1:ncol(AA), genenumber)
      resultM <<- mvn(AA[, probe], mvnTest = "mardia")$multivariateNormality
      resultH <<- mvn(AA[, probe], mvnTest = "hz")$multivariateNormality
      resultR <<- mvn(AA[, probe], mvnTest = "royston")$multivariateNormality
      resultN <<- mvnorm.etest(AA[, probe], R = 999)
      if(sum(c(resultH$`p value` < sig, resultR$`p value` < sig, resultN$p.value < sig)) >= 2){
        count2 = count2 + 1
      }
      
      if(as.numeric(as.character(resultM$`p value`[which.min(as.numeric(as.character(resultM$`p value`[-3])))])) < sig){
        countM = countM + 1
      }

      if(resultH$`p value` < sig){
        countH = countH + 1
        if(resultR$`p value` < sig){
          countHR = countHR + 1
          if(resultN$p.value < sig){
            countHRN = countHRN + 1
          }
        }
        if(resultN$p.value < sig){
          countHN = countHN + 1
        }
      }
      
      if(resultR$`p value` < sig){
        countR = countR + 1
        if(resultN$p.value < sig){
          countRN = countRN + 1
        }
      }
      
      if(resultN$p.value < sig){
        countN = countN + 1
      }
    }
    count = matrix(c(countM, countH, countR, countN, countHR, countHN, countRN, count2, countHRN), 1, length(names))
    return(count/times)
  }
  result = t(sapply(AAA, function(x) number(x)))
  colnames(result) = names
  rownames(result) = names(AAA)
  if(log2){
    write.csv(result, 
              file = (paste("log2_", as.character(deparse(substitute(AAA))), "_", genenumber, "probes_", times, "times_result", ".csv", sep = "")), row.names = TRUE)
    cat("\nFile saved as: ")
    cat(paste(getwd(), "/", "log2_" , as.character(deparse(substitute(AAA))), "_", genenumber, "probes_", times, "times_result", ".csv", sep = ""),"\n")
  }else{
    write.csv(result, 
              file = (paste(as.character(deparse(substitute(AAA))), "_", genenumber, "probes_", times, "times_result", ".csv", sep = "")), row.names = TRUE)
    cat("\nFile saved as: ")
    cat(paste(getwd(),"/",as.character(deparse(substitute(AAA))), "_", genenumber, "probes_", times, "times_result", ".csv", sep = ""),"\n")
  }
  
  return(result)
}
#############################################################################################################

setwd("C:/Users/jeff/Desktop/Study/Results")
#Classified by different Diseases
nonnormal(Breast, genenumber = 4, times = 500)
nonnormal(Colorectal, genenumber = 4, times = 500)
nonnormal(LungA, genenumber = 4, times = 500)
nonnormal(LungB, genenumber = 4, times = 500)
summary(Leukemia_CN)
nonnormal(Leukemia_CN, genenumber = 4, times = 500, adjacent = "CN", disease = "Leukemia")
nonnormal(Leukemia_CN, genenumber = 4, times = 500, adjacent = "CN", disease = "Leukemia", log2 = TRUE)
#nonnormal(Leukemia_CN, genenumber = 2, times = 500, adjacent = "CN", disease = "Leukemia")
summary(Leukemia_NC)
nonnormal(Leukemia_NC, genenumber = 4, times = 500, adjacent = "NC", disease = "Leukemia")
nonnormal(Leukemia_NC, genenumber = 4, times = 500, adjacent = "NC", disease = "Leukemia", log2 = TRUE)
#nonnormal(Leukemia_NC, genenumber = 2, times = 500, adjacent = "NC", disease = "Leukemia")
summary(Ovarian_CN)
nonnormal(Ovarian_CN, genenumber = 4, times = 500, adjacent = "CN", disease = "Ovarian")
nonnormal(Ovarian_CN, genenumber = 4, times = 500, adjacent = "CN", disease = "Ovarian", log2 = TRUE)
#nonnormal(Ovarian_CN, genenumber = 2, times = 500, adjacent = "CN", disease = "Ovarian")
summary(Ovarian_NC)
nonnormal(Ovarian_NC, genenumber = 4, times = 500, adjacent = "NC", disease = "Ovarian")
nonnormal(Ovarian_NC, genenumber = 4, times = 500, adjacent = "NC", disease = "Ovarian", log2 = TRUE)
#nonnormal(Ovarian_NC, genenumber = 2, times = 500, adjacent = "NC", disease = "Ovarian")


#Classified by unique probes(Genes)
nonnormal(NonUniqueProbe, genenumber = 4, times = 500)
nonnormal(UniqueProbe, genenumber = 4, times = 500)
#nonnormal(NonUniqueProbe, genenumber = 2, times = 500)
#nonnormal(UniqueProbe, genenumber = 2, times = 500)



#log
nonnormal(NonUniqueProbe, genenumber = 4, times = 500, log2 = TRUE)
nonnormal(UniqueProbe, genenumber = 4, times = 500, log2 = TRUE)




#All Real Data Analysis finish
#=================================================================================================
