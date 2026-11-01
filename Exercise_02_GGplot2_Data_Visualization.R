#===============================================================
# EXERCISE 02 - DATA VISUALIZATION WITH GGPLOT2 — Grammar of Graphs
# From Genome to Drug | 2025
# Università degli Studi di Perugia
#===============================================================
# Dataset: TCGA-LAML (The Cancer Genome Atlas — Acute Myeloid
#          Leukemia) clinical and survival data
# Goal:    Explore patient clinical variables using ggplot2
#          and build 11 different types of plots
#===============================================================

#---------------------------------------------------------------
# INSTALL AND LOAD REQUIRED LIBRARIES
#---------------------------------------------------------------
# The easiest way to get all tidyverse packages at once:
# install.packages("tidyverse")

# Alternatively, install individually:
# install.packages("dplyr")
# install.packages("ggplot2")

library(dplyr)
library(tidyverse)
library(ggplot2)

#---------------------------------------------------------------
# LOAD DATASETS
#---------------------------------------------------------------
# You will need two files:
#   - TCGA_LAML_clinical.txt   (clinical data)
#   - TCGA_LAML_survival.txt   (survival data)
# A file chooser window will open for each one.

clinical_data <- read.delim(file.choose(), header = TRUE, row.names = 1)
head(clinical_data)

survival_data <- read.delim(file.choose(), header = TRUE, row.names = 1)
head(survival_data)

# Remove the last column from survival data (contains redundant info)
survival_data <- survival_data[, -ncol(survival_data)]
head(survival_data)

#---------------------------------------------------------------
# PLOT 1 — Histogram: Age distribution at diagnosis
#---------------------------------------------------------------
# Question: How are patient ages distributed at the time
#           of initial diagnosis?

plot1 <- ggplot(clinical_data, aes(x = age_at_initial_pathologic_diagnosis)) +
  geom_histogram(binwidth = 5, fill = "#00A8FF", color = "black") +
  labs(title = "Patient Age Distribution at Diagnosis",
       x = "Age at Initial Diagnosis",
       y = "Count") +
  scale_x_continuous(
    breaks = seq(
      min(clinical_data$age_at_initial_pathologic_diagnosis, na.rm = TRUE),
      max(clinical_data$age_at_initial_pathologic_diagnosis, na.rm = TRUE),
      by = 10)) +
  theme_minimal()

plot1

#---------------------------------------------------------------
# PLOT 2 — Bar plot: Survival by tumor status
#---------------------------------------------------------------
# Question: Among patients with active vs. tumor-free status,
#           how many are alive vs. deceased?

plot2 <- ggplot(clinical_data, aes(x = person_neoplasm_cancer_status, fill = vital_status)) +
  geom_bar(position = "dodge2") +
  labs(title = "Survival by Tumor Status",
       x = "Tumor Status",
       y = "Number of Patients") +
  theme_minimal()

plot2

#---------------------------------------------------------------
# PLOT 3 — Boxplot: Tumor length by age group
#---------------------------------------------------------------
# Question: Does tumor size differ across age groups?
# We first create age group categories.

clinical_data <- clinical_data %>%
  mutate(age_group = cut(age_at_initial_pathologic_diagnosis,
                         breaks = c(0, 30, 60, Inf),
                         labels = c("Under 30", "31-60", "Above 60"),
                         right = FALSE))

plot3 <- ggplot(clinical_data, aes(x = age_group, y = pathologic_tumor_length, fill = age_group)) +
  geom_boxplot() +
  labs(title = "Tumor Length Distribution by Age Group",
       x = "Age Group",
       y = "Tumor Length") +
  scale_fill_manual(values = c("#E74C3C", "#3498DB", "#2ECC71")) +
  theme_minimal()

plot3

#---------------------------------------------------------------
# PLOT 4.1 — Scatter plot: Age vs. pathologic tumor length
#---------------------------------------------------------------
# Question: Is there a linear relationship between patient age
#           and tumor size at diagnosis?

plot4_1 <- ggplot(clinical_data, aes(x = age_at_initial_pathologic_diagnosis,
                                     y = pathologic_tumor_length)) +
  geom_point(alpha = 1, color = "blue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Age vs. Pathologic Tumor Length",
       x = "Age at Diagnosis",
       y = "Pathologic Tumor Length") +
  theme_minimal()

plot4_1

#---------------------------------------------------------------
# PLOT 4.2 — Scatter plot: Radiologic vs. pathologic tumor length
#---------------------------------------------------------------
# Question: How well does radiologic measurement predict
#           the actual pathologic tumor size?

# Ensure columns are numeric before plotting
clinical_data$radiologic_tumor_length <- as.numeric(as.character(clinical_data$radiologic_tumor_length))
clinical_data$pathologic_tumor_length <- as.numeric(as.character(clinical_data$pathologic_tumor_length))

plot4_2 <- ggplot(clinical_data, aes(x = radiologic_tumor_length,
                                     y = pathologic_tumor_length)) +
  geom_point(alpha = 1, color = "blue") +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Radiologic vs. Pathologic Tumor Length",
       x = "Radiologic Tumor Length",
       y = "Pathologic Tumor Length") +
  theme_minimal()

plot4_2

#---------------------------------------------------------------
# PLOT 5 — Bar plot: Number of patients by vital status
#---------------------------------------------------------------
# Question: How many patients are alive vs. deceased
#           in the TCGA-LAML cohort?

plot5 <- ggplot(clinical_data, aes(x = vital_status, fill = vital_status)) +
  geom_bar() +
  labs(title = "Number of Patients by Vital Status",
       x = "Vital Status",
       y = "Number of Patients") +
  theme_minimal()

plot5

#---------------------------------------------------------------
# PLOT 6 — Density plot: Age distribution by tumor status
#---------------------------------------------------------------
# Question: Do patients with active tumor differ in age
#           from those who are tumor-free?

plot6 <- ggplot(clinical_data, aes(x = age_at_initial_pathologic_diagnosis,
                                   fill = person_neoplasm_cancer_status)) +
  geom_density(alpha = 0.5) +
  labs(title = "Age Distribution by Tumor Status",
       x = "Age at Diagnosis",
       y = "Density") +
  theme_minimal()

plot6

#---------------------------------------------------------------
# PLOT 7 — Violin + Boxplot: Age distribution by tumor status
#---------------------------------------------------------------
# Question: How does the full distribution of age compare
#           between tumor status groups?
# Note: Violin plots show the full distribution shape,
#       while the inner boxplot shows median and IQR.

plot7 <- ggplot(clinical_data, aes(x = person_neoplasm_cancer_status,
                                   y = age_at_initial_pathologic_diagnosis,
                                   fill = person_neoplasm_cancer_status)) +
  geom_violin(alpha = 0.5, width = 1.0) +
  geom_boxplot(width = 0.1, outlier.color = "black", fill = "white") +
  labs(title = "Age Distribution by Tumor Status",
       x = "Tumor Status",
       y = "Age at Diagnosis") +
  theme_minimal()

plot7

#---------------------------------------------------------------
# PLOT 8 — Facet wrap: Age vs. tumor length by vital status
#---------------------------------------------------------------
# Question: Does the relationship between age and tumor size
#           differ between patients who are alive vs. deceased?

plot8 <- ggplot(clinical_data, aes(x = age_at_initial_pathologic_diagnosis,
                                   y = pathologic_tumor_length,
                                   color = vital_status)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~ vital_status) +
  labs(title = "Age vs. Tumor Length by Vital Status",
       x = "Age at Diagnosis",
       y = "Tumor Length") +
  theme_minimal()

plot8

#---------------------------------------------------------------
# PLOT 9 — Density plot: Survival time distribution
#---------------------------------------------------------------
# Question: How is overall survival time distributed
#           between patients who are alive vs. deceased?

plot9 <- ggplot(survival_data, aes(x = OS.time, color = OS)) +
  geom_density() +
  labs(title = "Overall Survival Time Distribution",
       x = "Days to Event",
       y = "Density") +
  theme_minimal()

plot9

#---------------------------------------------------------------
# PLOT 10 — Pie chart: Gender distribution
#---------------------------------------------------------------
# Question: What is the male-to-female ratio in the cohort?

gender_counts <- clinical_data %>%
  group_by(gender) %>%
  summarise(n = n())

plot10 <- ggplot(gender_counts, aes(x = "", y = n, fill = gender)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(title = "Gender Distribution") +
  theme_minimal()

plot10

#---------------------------------------------------------------
# PLOT 11 — Boxplot: Tumor length by treatment outcome
#---------------------------------------------------------------
# Question: Does tumor size differ based on whether
#           additional treatment was completed successfully?
# Note: Make sure clinical_data_clean is defined before
#       running this plot (remove NA rows if needed):
# clinical_data_clean <- clinical_data %>%
#   filter(!is.na(additional_treatment_completion_success_outcome))

plot11 <- ggplot(clinical_data_clean,
                 aes(x = additional_treatment_completion_success_outcome,
                     y = pathologic_tumor_length,
                     fill = additional_treatment_completion_success_outcome)) +
  geom_boxplot() +
  labs(title = "Pathologic Tumor Length by Treatment Outcome",
       x = "Additional Treatment Completion Outcome",
       y = "Pathologic Tumor Length") +
  scale_fill_manual(values = c("#E74C3C", "#3498DB", "#2ECC71", "#F39C12", "#8E44AD")) +
  theme_minimal()

plot11

#===============================================================
# END OF SCRIPT
#===============================================================
