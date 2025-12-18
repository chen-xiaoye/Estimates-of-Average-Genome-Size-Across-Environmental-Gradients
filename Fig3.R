library(vegan)
library(ggplot2)
library(dplyr)
library(ggtext)
library(scales)
library(tidyr)
library(tibble)

KO <- read_excel("Supplcompst.xlsx", sheet = "KO_table") %>%
  column_to_rownames(var = colnames(.)[1])
grp <- read_excel("Supplenv.xlsx", sheet = "group")
bac_otu <- read_excel("Supplcompst.xlsx", sheet = "otutable") %>%
  column_to_rownames(var = colnames(.)[1])
bac_otu <- t(bac_otu)

KO_r <- as.data.frame(t(KO))
KO_richness <- round(KO_r)


alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- vegan::diversity(x, index = 'shannon', base = base)
  Simpson <- vegan::diversity(x, index = 'simpson', base = base)
  Pielou <- Shannon / log(ACE, base)
  observed_OTU <- rowSums(x>0)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, observed_OTU, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}


ko_alpha <- alpha(KO_richness, base = 2)

ko_alpha$samplename <- rownames(ko_alpha)
ko_alpha <- ko_alpha[order(ko_alpha$samplename),]

ko_comb <- merge(ko_alpha, grp, by.x = "samplename", by.y = "sampleid")


bac_alpha <- alpha(bac_otu, base = 2)
bac_alpha$samplename <- rownames(bac_alpha)

colnames(bac_alpha)[2] <- "amplicon_shannon"
ko_comb <- merge(ko_comb, bac_alpha[,c(2,9)], by = "samplename")


####Fig3a####

ko_comb$Shannon_norm <- rescale(ko_comb$Shannon, to = c(0, 1))
ko_comb$amplic_shannon_norm <- rescale(ko_comb$amplicon_shannon, to = c(0, 1))
norm <- ko_comb[,c(1,19,26,27)]


norm1 <- norm %>% mutate(
  grp = if_else(Shannon_norm > amplic_shannon_norm, "FA", "AF")
  
)
df <- norm1 %>%
  pivot_longer(
    cols = c(Shannon_norm, amplic_shannon_norm), 
    names_to = "types",
    values_to = "values"
  )


model <- lm(amplic_shannon_norm ~ pH + I(pH^2), data = ko_comb)
r_square <- summary(model)$r.squared
p <- summary(model)$coefficients[2,4]
apl_label_text <- ifelse(p < 0.001,
                         paste(
                           sprintf("R² = %.3f", r_square),
                           sprintf("P = %.3e", p),  # 使用 sprintf 格式化 P 值
                           sep = ", "
                         ),
                         paste(
                           sprintf("R² = %.3f", r_square),
                           sprintf("P = %.3f", p),
                           sep = ", "
                         ))

model <- lm(Shannon_norm ~ pH + I(pH^2), data = ko_comb)
r_square <- summary(model)$r.squared
p <- summary(model)$coefficients[2,4]
meta_label_text <- ifelse(p < 0.001,
                          paste(
                            sprintf("R² = %.3f", r_square),
                            sprintf("P = %.3e", p),  # 使用 sprintf 格式化 P 值
                            sep = ", "
                          ),
                          paste(
                            sprintf("R² = %.3f", r_square),
                            sprintf("P = %.3f", p),
                            sep = ", "
                          ))


library(ggnewscale)
ggplot(df, aes(x = pH, y = values, color = types)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "Shannon_norm" = "#E41A1C",         
      "amplic_shannon_norm" = "#377EB8" ),
    labels = c("Shannon_norm" = "H'.fun.scaled", "amplic_shannon_norm" = "H'.tax.scaled"))+
  geom_smooth(se = F, method = "lm",formula = y ~ poly(x, 2, raw = TRUE),size = 2) +  
  new_scale_color()+
  geom_path(aes(group = samplename, color = grp),linewidth = 0.8) +
  scale_color_manual(values = c("#377EB8","#E41A1C" )) +
  guides(color = "none")  +
  labs(x = "pH", y = "Shannon diversity (scaled)", subtitle = paste0(apl_label_text," ",meta_label_text)) +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5),
        legeng.position = "right")

####Fig3b####
col1 <- c("blue","#336dff","#1E90FF","#80B1DA","#ffcc33","#ffb133","lightcoral","red")
ko_comb$dis <- ko_comb$Shannon_norm - ko_comb$amplic_shannon_norm

cor <-  cor.test(ko_comb$pH, ko_comb$dis, method = "spearman")
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


ggplot(ko_comb, aes(x = pH, y = dis, color = pH)) +
  geom_point(size = 5) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "pH", y = "H'.fun.scaled - H'.tax.scaled (ΔH')", subtitle = label_text) +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))
####Fig3c####
cor <-  cor.test(ko_comb$dis, ko_comb$metagenome_AGS, method = "spearman")
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
ggplot(ko_comb, aes(x = dis, y = metagenome_AGS, color = pH)) +
  geom_point(size = 5) +
  scale_y_continuous(labels = label_number(scale = 1e-6)) + 
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "ΔH'", y = "Metagenomic AGS", subtitle = label_text) +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))
####Fig3d####
cor <-  cor.test(ko_comb$dis, ko_comb$gtdb_AGS, method = "spearman")
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

ggplot(ko_comb, aes(x = dis, y = gtdb_AGS, color = pH)) +
  geom_point(size = 5) +
  scale_y_continuous(labels = label_number(scale = 1e-6)) + 
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, values = scales::rescale(c(8, 8.5, 9.0, 9.5,10,10.5)),name = "pH") +
  labs(x = "ΔH'", y = "16S metabarcoding AGS", subtitle = label_text) +
  #ggtitle("KO richness") +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))

####Fig3e####
comb <- read_excel("Supplenv.xlsx", sheet = "Feng_group")

comb$Shannon_meta_norm <- rescale(comb$ko_shannon, to = c(0, 1))
comb$Shannon_ampl_norm <- rescale(comb$amplic_shannon, to = c(0, 1))
comb$dis <- comb$Shannon_meta_norm - comb$Shannon_ampl_norm
comb$log2_Ec <- log2(comb$Ec)

norm <- comb[,c(1,17,18,20)]

library(tidyr)

norm1 <- norm %>% mutate(
  grp = if_else(Shannon_meta_norm > Shannon_ampl_norm, "FA", "AF")
  
)
df <- norm1 %>%
  pivot_longer(
    cols = c(Shannon_meta_norm, Shannon_ampl_norm), 
    names_to = "types",
    values_to = "values"
  )

cor <-  cor.test(comb$log2_Ec, comb$Shannon_ampl_norm, method = "spearman")
p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
ampl_label_text <- ifelse(p_value < 0.001,
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

cor <-  cor.test(comb$log2_Ec, comb$Shannon_meta_norm, method = "spearman")
p_value <- cor$p.value
p_value <- p.adjust(p_value, method = "fdr")
R <- cor$estimate
meta_label_text <- ifelse(p_value < 0.001,
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


ggplot(df, aes(x = log2_Ec, y = values, color = types)) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      "Shannon_meta_norm" = "#E41A1C",        
      "Shannon_ampl_norm" = "#377EB8" ),
    labels = c("Shannon_meta_norm" = "H'.fun.scaled", "Shannon_ampl_norm" = "H'.tax.scaled"))+
  geom_smooth(se = F, method = "lm",size = 2) +  
  new_scale_color()+
  geom_path(aes(group = samplename, color = grp),linewidth = 0.8) +
  scale_color_manual(values = c("#377EB8","#E41A1C" )) +
  guides(color = "none")  +
  labs(x = "log2(EC)", y = "Shannon diversity (scaled)", subtitle = paste0(ampl_label_text," ",meta_label_text)) +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))

####Fig3f####
col1 <- c("#3f5efb","#5b5ae5","#6b59da","#9653ba","#b150a5","#cb4c90","#dd4a83","#fc466b")
cor <-  cor.test(comb$log2_Ec, comb$dis, method = "spearman")
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


ggplot(comb, aes(x = log2_Ec, y = dis, color = log2_Ec)) +
  geom_point(size = 5) +
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, name = "log2(EC)") +
  labs(x = "log2(EC)", y = "ΔH'", subtitle = label_text) +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))

####Fig3g####
cor <-  cor.test(comb$dis, comb$Bac_GS, method = "spearman")
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


ggplot(comb, aes(x = dis, y = Bac_GS, color = log2_Ec)) +
  geom_point(size = 5) +
  scale_y_continuous(labels = label_number(scale = 1e-6)) + 
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, name = "log2(EC)") +
  labs(x = "ΔH'", y = "Metagenomic AGS", subtitle = label_text) +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))


####Fig3h####
cor <-  cor.test(comb$dis, comb$gtdb_genomesize, method = "spearman")
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


ggplot(comb, aes(x = dis, y = gtdb_genomesize, color = log2_Ec)) +
  geom_point(size = 5) +
  scale_y_continuous(labels = label_number(scale = 1e-6)) + 
  scale_size_continuous(range = c(3,10),guide = "none") +
  geom_smooth(color = "white", se = TRUE, aes(group = 1), method = "lm") +  # 线性拟合
  scale_color_gradientn(colors = col1, name = "log2(EC)") +
  labs(x = "ΔH'", y = "16S metabarcoding based AGS", subtitle = label_text) +
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 18, hjust = 0.5, color = "black"),
        axis.title.y = element_text(colour = "black", size = 20, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 20, face = "bold"),
        axis.text.x = element_text(size = 16, color = "black",face = "bold"),
        axis.text.y = element_text(size = 16, color = "black",face = "bold"),
        legend.title = element_text(size = 20, face = "bold", color = "black", hjust=0),
        legend.text = element_text(size = 18, face = "bold", color = "black"),
        plot.title = element_text(size = 25, face = "bold",hjust = 0.5))
