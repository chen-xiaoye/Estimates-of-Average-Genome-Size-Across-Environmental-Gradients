library(vegan)
library(ggplot2)
library(scales)
library(agricolae)
library(ggtext)
library(dplyr)

all0 <- read_excel("Supplenv.xlsx", sheet = "AGS_both")



plotcol <- c("#3286BC","red")

all0$Type <- factor(all0$Type, levels = c("Wang et al.","this study"))


this_stdy <- subset(all0, Type == "this study")

model <- lm(genomesize ~ pH, data = this_stdy)
r_square <- summary(model)$r.squared
p <- summary(model)$coefficients[2,4]
label_text_1 <- ifelse(p < 0.001,
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
Wang <- subset(all0, Type == "Wang et al.")
model <- lm(genomesize ~ pH, data = Wang)
r_square <- summary(model)$r.squared
p <- summary(model)$coefficients[2,4]
label_text_2 <- ifelse(p < 0.001,
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


ggplot(all0, aes(x = pH, y = genomesize, color = Type, shape = Type)) +
  geom_point(size = 5) +
  geom_smooth(method = "lm", color = "white", se = TRUE,aes(group = Type)) +  # 线性拟合
  scale_color_manual(name = "Type", values = plotcol)+
  scale_shape_manual(name = "Type", values=c(15,16)) +
  labs(x = "pH", y = "Average genome size (Mbp)",subtitle = paste0(label_text_1," ",label_text_2)) +
  scale_x_continuous(limits = c(3,11), breaks = seq(3,11,1)) +
  guides(color = guide_legend(order = 1),  # Habitat 图例的顺序为 1
         shape = "none")+ # Location 图例的顺序为 2
  theme_bw()+
  theme(plot.subtitle = element_markdown(size = 15, hjust = 0.5, color = "black"),
        plot.title = element_text(size = 20, face = "bold",hjust = 0.5),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 17, face = "bold"))
