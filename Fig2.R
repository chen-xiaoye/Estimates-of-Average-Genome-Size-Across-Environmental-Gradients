library(vegan)
library(ggplot2)
library(scales)
library(agricolae)
library(ggtext)
library(readxl)
####Fig2-a####
total <- read_excel("Supplenv.xlsx", sheet = "group")


cor <-  cor.test(total$pH, total$metagenome_AGS, method = "spearman")
p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
label_text <- ifelse(p_value < 0.001,
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3e", p_value),  # 使用 sprintf 格式化 P 值
                       sep = ", "
                     ),
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3f", p_value),
                       sep = ", "
                     ))


ggplot(total, aes(x = pH, y = metagenome_AGS, color = pH, size = metagenome_AGS)) +
  geom_point() +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#2e6cae", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "Metagenomes",x = "pH", y = "Average genome size (Mbp)", subtitle = label_text) +
  scale_y_continuous(labels = label_number(scale = 1e-6)) + 
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5),
        legend.position = "none")


####Fig2-b####

cor <-  cor.test(total$pH, total$gtdb_AGS, method = "spearman")
p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
label_text <- ifelse(p_value < 0.001,
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3e", p_value),  # 使用 sprintf 格式化 P 值
                       sep = ","
                     ),
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3f", p_value),
                       sep = ","
                     ))


ggplot(total, aes(x = pH, y = gtdb_AGS, color = pH, size = gtdb_AGS)) +
  geom_point() +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#2e6cae", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(title = "16S metabarcoding", x = "pH", y = "Average Genome size (Mbp)", subtitle = label_text) +
  scale_y_continuous(
    limits = c(1500000, 5000000),          # 轴范围
    breaks = seq(1500000, 5000000, 500000),
    labels = label_number(scale = 1e-6) # 刻度间隔（0, 2, 4, 6, 8, 10）
  ) +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5),
        legend.position = "none")

####Fig2-c####
g_total <- read_excel("Supplenv.xlsx", sheet = "Feng_group")

g_total$Ec_log2 <- log2(g_total$Ec)
cor <- cor.test(g_total$Ec_log2, g_total$Bac_GS, method = "spearman")
p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
label_text <- ifelse(p_value < 0.001,
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3e", p_value),  # 使用 sprintf 格式化 P 值
                       sep = ", "
                     ),
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3f", p_value),
                       sep = ", "
                     ))

ggplot(g_total, aes(x = Ec_log2, y = Bac_GS, color = Ec_log2, size = Bac_GS)) +
  geom_point() +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#4a9e48", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "log2(EC)") +
  labs(title = "Metagenomes",x = "log2(EC)", y = "Average genome size (Mbp)", subtitle = label_text) +
  scale_y_continuous(labels = label_number(scale = 1e-6)) + 
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5),
        legend.position = "none")


####Fig2-d####
cor <- cor.test(g_total$Ec_log2, g_total$gtdb_genomesize, method = "spearman")
p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
label_text <- ifelse(p_value < 0.001,
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3e", p_value),  # 使用 sprintf 格式化 P 值
                       sep = ", "
                     ),
                     paste(
                       sprintf("R = %.3f", R),
                       sprintf("P = %.3f", p_value),
                       sep = ", "
                     ))

ggplot(g_total, aes(x = Ec_log2, y = gtdb_genomesize, color = Ec_log2, size = gtdb_genomesize)) +
  geom_point() +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(method = "lm", color = "white", se = TRUE, aes(group = 1)) +  # 线性拟合
  scale_color_gradientn(colors = "#4a9e48", values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "log2(EC)") +
  labs(title = "16S metabarcoding",x = "log2(EC)", y = "Average genome size (Mbp)", subtitle = label_text) +
  scale_y_continuous(labels = label_number(scale = 1e-6)) + 
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5),
        legend.position = "none")
