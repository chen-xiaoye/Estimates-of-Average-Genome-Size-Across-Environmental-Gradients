library(vegan)
library(ggplot2)
library(scales)
library(agricolae)
library(ggtext)
library(dplyr)

####Fig5a####
otu_all <- read_excel("Supplcompst.xlsx", sheet = "otutable")
gsize <- read_excel("Supplenv.xlsx", sheet = "genomesize")
total <- read_excel("Supplenv.xlsx", sheet = "group")

otu_all <- merge(otu_all, gsize, by = "OTUID", all = TRUE)
otu_all <- as.data.frame(otu_all)

df1 <- otu_all %>%
  mutate(count = if_else(!is.na(genomesize), 1, 0))

sample_cols <- names(df1)[2:41]
result_list <- list()

for (sample in sample_cols) {
  temp <- df1 %>%
    select(OTUID, all_of(sample), count) %>%
    rename(abundance = all_of(sample)) %>%
    filter(abundance > 0) %>%  # 只考虑在该样品中出现的OTU
    group_by(count) %>%
    summarise(otu_number = n()) %>%
    mutate(
      perctge = otu_number / sum(otu_number),
      mod = sample,
      covery = ifelse(count == 0, "unmatched", "matched")
    ) %>%
    select(mod, covery, perctge)
  
  result_list[[sample]] <- temp
}

otu_sample_summary <- bind_rows(result_list)

unmatched <- otu_sample_summary[otu_sample_summary$covery == "unmatched",]

total <- merge(total, unmatched, by.x = "sampleid", by.y = "mod")

cor <-  cor.test(total$pH, total$perctge, method = "spearman")
p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
unmatched_text <- ifelse(p_value < 0.001,
                         paste(
                           sprintf("R = %.3f", R),
                           sprintf("P = %.3e", p_value),  # 使用 sprintf 格式化 P 值
                           sep = " "
                         ),
                         paste(
                           sprintf("R = %.3f", R),
                           sprintf("P = %.3f", p_value),
                           sep = " "
                         ))


ggplot(total, aes(x = pH,y = perctge,color = "pH")) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", color = "white", size =3 ,se = TRUE, group = 1) +
  scale_color_manual(values = "#8B008B",name = "pH") +
  scale_y_continuous(
    limits = c(0, 0.5),
    name = "Unmatched proportion (%)",
    labels = scales::percent_format(accuracy = 1))+
  labs(
    x = "pH",
    subtitle = unmatched_text
  ) +
  theme_classic() +
  theme(
    plot.subtitle = element_markdown(size = 15, hjust = 0.5, color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 17, face = "bold"),
    legend.position = "none")

####Fig5b####
rm(list = ls())
otu_all <- read_excel("Supplcompst.xlsx", sheet = "Feng_otutable")
gsize <- read_excel("Supplenv.xlsx", sheet = "Feng_genomesize")
total <- read_excel("Supplenv.xlsx", sheet = "Feng_group")
otu_all <- merge(otu_all, gsize, by = "OTUID", all = TRUE)

otu_all <- as.data.frame(otu_all)

df1 <- otu_all %>%
  mutate(count = if_else(!is.na(genomesize), 1, 0))

sample_cols <- names(df1)[2:38]
result_list <- list()

for (sample in sample_cols) {
  temp <- df1 %>%
    select(OTUID, all_of(sample), count) %>%
    rename(abundance = all_of(sample)) %>%
    filter(abundance > 0) %>%  # 只考虑在该样品中出现的OTU
    group_by(count) %>%
    summarise(otu_number = n()) %>%
    mutate(
      perctge = otu_number / sum(otu_number),
      mod = sample,
      covery = ifelse(count == 0, "unmatched", "matched")
    ) %>%
    select(mod, covery, perctge)
  
  result_list[[sample]] <- temp
}

otu_sample_summary <- bind_rows(result_list)

unmatched <- otu_sample_summary[otu_sample_summary$covery == "unmatched",]

total <- merge(total, unmatched, by.x = "samplename", by.y = "mod")

total$Ec_log2 <- log2(total$Ec)

cor <- cor.test(total$Ec_log2, total$perctge, method = "spearman")

p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
unmatched_text <- ifelse(p_value < 0.001,
                         paste(
                           sprintf("R = %.3f", R),
                           sprintf("P = %.3e", p_value),  # 使用 sprintf 格式化 P 值
                           sep = " "
                         ),
                         paste(
                           sprintf("R = %.3f", R),
                           sprintf("P = %.3f", p_value),
                           sep = " "
                         ))


ggplot(total, aes(x = Ec_log2,y = perctge,color = "log2(EC)")) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", color = "white", size =3 ,se = TRUE, group = 1) +
  scale_color_manual(values = "#8B008B",name = "Log2(EC)") +
  scale_y_continuous(
    limits = c(0, 0.5),
    name = "Unmatched proportion (%)",
    labels = scales::percent_format(accuracy = 1))+
  labs(
    x = "log2(EC)",
    subtitle = unmatched_text
  ) +
  theme_classic() +
  theme(
    plot.subtitle = element_markdown(size = 15, hjust = 0.5, color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 17, face = "bold"),
    legend.position = "none")
