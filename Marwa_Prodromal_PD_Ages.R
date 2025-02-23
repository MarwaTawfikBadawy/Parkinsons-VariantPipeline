---
title: "gVCF AGES visualization"
author: "Marwa Tawfik Badawy"
date: "2023-07-12"
output: html_document
---

$$Marwa-PhD-Bioinformatics-Analysis-Workflow$$ <br>
**1- Downloading the data from PPMI server with my approved access** <br>
**2- Storing them on 3 different external hard drives to be used later on** <br>
**3- Creating excel sheet containing all of the information for each group and subgroup of the patients and healthy individuals** <br>
**4- Looking for the codes as represented in this markdown to get overall data visualization** <br>
**5- Starting unzipping the compressed files in terminal by using these two codes** <br>
**5-A- sudo mount and the Path-name of the used external hard** <br>
**5-B- gunzip followed by the Path-name of each file (patient ID)** <br>
*P.S. If the "gunzip" code is not working, we can use the code "open" followed by the Path-name of the folder that containing the required files*


**Loading Libraries**
```{r}
library(tidyverse)
library(rmarkdown)
library(knitr)
library(fastqcr)
library(fastmap)
library(dplyr)
library(plyr)
library(here)
library(ggplot2)
```

*Loading HC ages*
```{r}
library(readxl)
Marwa_HC <- read_excel("/Users/mahamarwatawfik/WGS_gVCF/PPMI_Study/Data_Ages/Marwa_HC_Data.xlsx")
```

```{r}
names(Marwa_HC)
```

```{r}
group_by(Marwa_HC)
```

```{r}
table(Marwa_HC$Subgroup)
```


```{r}
table(Marwa_HC$Sex)
table(Marwa_HC$Age)
```

```{r}
library(dplyr)
data.frame(Marwa_HC) %>% group_by(Age)
```

*Classifying HC participants based on their age groups*

```{r}
Marwa_HC$HC_Age_Group <- cut(Marwa_HC$Age,
                         breaks = c(-Inf
                                    ,5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60
                                    , 65, 70, 75, 80, 85, Inf),
                         labels = c("0-4 years", "5-9 years", "10-14 years"
                                    , "15-19 years","20-24 years", "25-29 years"
                                    ,"30-34 years", "35-39 years", "40-44 years"
                                    , "45-49 years", "50-54 years", "55-59 years"
                                    , "60-64 years", "65-69 years", "70-74 years"
                                    , "75-79 years", "80-84 years", "85+ years"),
                         right = FALSE)

```

```{r}
View(Marwa_HC)
```

*Plotting HC age groups*
```{r}
P1 <- ggplot(data = na.omit(Marwa_HC),aes(x = HC_Age_Group, colour = Sex,
                                               fill = Sex)) + ggtitle("HC_Age_Groups") +
  theme(plot.title = element_text(color = "red", size = 18, face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5),  # Center the title horizontally
        plot.margin = unit(c(1, 1, 2, 1), "lines"))  # Adjust margins
P1 <- P1+geom_bar()
P1 
```

*Loading Pro Genetic Patients ages*

```{r}
library(readxl)
Marwa_Pro_Gene <- read_excel("/Users/mahamarwatawfik/WGS_gVCF/PPMI_Study/Data_Ages/Marwa_Pro_Gene.xlsx")
```

```{r}
names(Marwa_Pro_Gene)
```


```{r}
group_by(Marwa_Pro_Gene)
```


```{r}
table(Marwa_Pro_Gene$Sex)
table(Marwa_Pro_Gene$Age)
```

*Classifying Pro Genetic patients based on their age groups*

```{r}
Marwa_Pro_Gene$Pro_Gene_Age_Group <- cut(Marwa_Pro_Gene$Age,
                         breaks = c(-Inf
                                    ,5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60
                                    , 65, 70, 75, 80, 85, Inf),
                         labels = c("0-4 years", "5-9 years", "10-14 years"
                                    , "15-19 years","20-24 years", "25-29 years"
                                    ,"30-34 years", "35-39 years", "40-44 years"
                                    , "45-49 years", "50-54 years", "55-59 years"
                                    , "60-64 years", "65-69 years", "70-74 years"
                                    , "75-79 years", "80-84 years", "85+ years"),
                         right = FALSE)

```

```{r}
view(Marwa_Pro_Gene)
```


*Plotting Pro Genetic age groups*
```{r}
P2 <- ggplot(data = na.omit(Marwa_Pro_Gene),aes(x = Pro_Gene_Age_Group, colour = Sex,
                                              fill = Sex)) + ggtitle("Pro_Gen_Age_Groups") +
  theme(plot.title = element_text(color = "red", size = 20, face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5),  # Center the title horizontally
        plot.margin = unit(c(1, 1, 2, 1), "lines"))  # Adjust margins
P2 <- P2+geom_bar()
P2
```
*Loading Healthy Control ages*

```{r}
library(readxl)
Marwa_Pro_Hypo <- read_excel("/Users/mahamarwatawfik/WGS_gVCF/PPMI_Study/Data_Ages/Marwa_Pro_Hypo.xlsx")
```

```{r}
names(Marwa_Pro_Hypo)
```

```{r}
table(Marwa_Pro_Hypo$Sex)
table(Marwa_Pro_Hypo$Age)
```

*Classifying Healthy Control based on their age groups*

```{r}
Marwa_Pro_Hypo$Pro_Hypo_Age_Group <- cut(Marwa_Pro_Hypo$Age,
                         breaks = c(-Inf
                                    ,5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60
                                    , 65, 70, 75, 80, 85, Inf),
                         labels = c("0-4 years", "5-9 years", "10-14 years"
                                    , "15-19 years","20-24 years", "25-29 years"
                                    ,"30-34 years", "35-39 years", "40-44 years"
                                    , "45-49 years", "50-54 years", "55-59 years"
                                    , "60-64 years", "65-69 years", "70-74 years"
                                    , "75-79 years", "80-84 years", "85+ years"),
                         right = FALSE)

```

```{r}
view(Marwa_Pro_Hypo)
```


*Plotting Healthy Control age groups*
```{r}
P3 <- ggplot(data = na.omit(Marwa_Pro_Hypo),aes(x = Pro_Hypo_Age_Group, colour = Sex,
                                                fill = Sex)) + ggtitle("Pro_Hypo_Age_Groups") +
  theme(plot.title = element_text(color = "red", size = 20, face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5),  # Center the title horizontally
        plot.margin = unit(c(1, 1, 2, 1), "lines"))  # Adjust margins
P3 <- P3+geom_bar()
P3
```

*Loading Prodromal Genetic Patients ages*

```{r}
library(readxl)
Marwa_Pro_RBD <- read_excel("/Users/mahamarwatawfik/WGS_gVCF/PPMI_Study/Data_Ages/Marwa_Pro_RBD.xlsx")
```

```{r}
names(Marwa_Pro_RBD)
```

```{r}
table(Marwa_Pro_RBD$Sex)
table(Marwa_Pro_RBD$Age)
```

*Classifying Prodromal Genetic patients based on their age groups*

```{r}
Marwa_Pro_RBD$Pro_RBD_Age_Group <- cut(Marwa_Pro_RBD$Age,
                         breaks = c(-Inf
                                    ,5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60
                                    , 65, 70, 75, 80, Inf),
                         labels = c("0-4 years", "5-9 years", "10-14 years"
                                    , "15-19 years","20-24 years", "25-29 years"
                                    ,"30-34 years", "35-39 years", "40-44 years"
                                    , "45-49 years", "50-54 years", "55-59 years"
                                    , "60-64 years", "65-69 years", "70-74 years"
                                    , "75-79 years", "80-84 years"),
                         right = FALSE)

```

```{r}
view(Marwa_Pro_RBD)
```

*Plotting Prodromal genetic age groups*
```{r}
P4 <- ggplot(data = na.omit(Marwa_Pro_RBD),aes(x = Pro_RBD_Age_Group, colour = Sex,
                                                fill = Sex)) + ggtitle("Pro_RBD_Age_Groups") +
  theme(plot.title = element_text(color = "red", size = 20, face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.5),  # Center the title horizontally
      plot.margin = unit(c(1, 1, 2, 1), "lines"))  # Adjust margins

P4 <- P4+geom_bar()
P4

```





