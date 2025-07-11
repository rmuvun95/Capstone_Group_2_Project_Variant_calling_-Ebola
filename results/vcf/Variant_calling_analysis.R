#==================================================================#
# Description: R analysis for variant calling analysis
# Date: 2025-July-10
# Author: Raissa, John, Rennatha, Liya 
#==================================================================#

# Installing packages - only done once

install.packages(c("openxlsx", "dplyr", "janitor", "stringi", "stringr", "ggplot2", "data.table"))

# Loading all the required libraries

library(openxlsx)
library(dplyr)
library(janitor)
library(stringi)
library(stringr)
library(ggplot2)
library(data.table)

#====================================================================#

# Loading the variants file
df <- fread("hq_allvariants.tsv")

# Getting to know about your data table
dim (df) # prints out the number of rows and columns
colnames(df) # prints out columns name
unique(df$SAMPLE_ID) # Prints out unique sample IDs

#====================================================================#
# Calculating the total number of high-quality variants per sample
number_hq_per_sample <- df %>%
  clean_names() %>%ß
  group_by(sample_id) %>%
  summarise(n_hq = n(), .groups = "drop")
  

View(number_hq_per_sample)

# Find out which kind of variants are in each sample
variants_type <- df %>%
  clean_names() %>%
  mutate(variant_type = case_when(
    (ref == "A" & alt == "G") | (ref == "G" & alt == "A") ~ "Transition",
    (ref == "C" & alt == "T") | (ref == "T" & alt == "C") ~ "Transition",
    (ref %in% c("A", "G", "C", "T") & alt %in% c("A", "G", "C", "T")) ~ "Transversion",
    TRUE ~ "Other"
  ))
        
# Finding out the number of variant types and proportions in each sample

prop_variants_type <- variants_type %>%
  group_by(sample_id, variant_type) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(prop = count / sum(count) * 100) %>%
  ungroup()

# Finding out the most common variant types

variants_type <- df %>%
  clean_names() %>%
  mutate(
    substitution = paste0(ref, ">", alt),  # Create substitution string like "A>G"
    variant_type = case_when(
      (ref == "A" & alt == "G") | (ref == "G" & alt == "A") ~ "Transition",
      (ref == "C" & alt == "T") | (ref == "T" & alt == "C") ~ "Transition",
      (ref %in% c("A", "T", "C", "G") & alt %in% c("A", "T", "C", "G")) ~ "Transversion",
      TRUE ~ "Other"
    )
  )

# Count how often each substitution occurs by variant type
common_variants <- variants_type %>%
  count(variant_type, substitution, sort = TRUE)

# Display unique substitutions

unique_variants <- variants_type %>%
  select(substitution) %>%
  distinct()

#====================================================================#

# Visualization
# Plotting the Number of High-Quality Variants in each sample

plot_numbers <- ggplot(number_hq_per_sample, aes(x = sample_id, y = n_hq, fill = sample_id)) +
  geom_bar(stat = "identity") +
  xlab("Sample ID") +
  ylab("Number of High-Quality Variants") +
  ggtitle("HQ Variants per Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = FALSE)  # Hide legend if not needed


plot_numbers

# Plotting the proportion of Ts vs. Tv variants per sample

plot_proportions <- ggplot(prop_variants_type, aes(x = sample_id, y = prop, fill = variant_type)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(prop, 1), "%")), 
            position = position_stack(vjust = 0.5), size = 3, color = "black") +
  xlab("Sample ID") +
  ylab("Proportion (%)") +
  ggtitle("Transitions vs Transversions per Sample") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Transition" = "darkblue", "Transversion" = "orange")) +
  guides(fill = guide_legend(title = "Variant Type"))

plot_proportions

# Optional: Looking at the susbtitutions
plot_substitutions <- ggplot(common_variants, aes(x = reorder(substitution, n), y = n, fill = variant_type)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # horizontal bars for better readability
  xlab("Substitution") +
  ylab("Count") +
  ggtitle("Counts of Substitutions by Variant Type") +
  theme_minimal() +
  scale_fill_manual(values = c("Transition" = "skyblue", "Transversion" = "orange")) +
  guides(fill = guide_legend(title = "Variant Type"))

plot_substitutions



