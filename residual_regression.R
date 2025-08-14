# Load libraries
library(ggplot2)
library(dplyr)
library(broom)

############################################ATAC-seq residuals regression analysis####################################################################
# Load ATAC data
df1 <- read.table("ATAC_OCT4_data.tsv", header = TRUE, sep = "\t")

# Create long-format data
long_df1 <- bind_rows(
  df1 %>% filter(RT == "ES") %>%
    transmute(RT = "ES",
              log2_off = log2(ATAC_OCT4.OFF),
              log2_on  = log2(ATAC_OCT4.ON)),
  df1 %>% filter(RT == "MS") %>%
    transmute(RT = "MS",
              log2_off = log2(ATAC_OCT4.OFF),
              log2_on  = log2(ATAC_OCT4.ON)),
  df1 %>% filter(RT == "LS") %>%
    transmute(RT = "LS",
              log2_off = log2(ATAC_OCT4.OFF),
              log2_on  = log2(ATAC_OCT4.ON))
) %>%
  mutate(residual = log2_on - log2_off)

# Linear model summaries
model_stats <- long_df1 %>%
  group_by(RT) %>%
  do(tidy(lm(residual ~ log2_on, data = .))) %>%
  filter(term == "log2_on") %>%
  mutate(label = sprintf("%s: Slope = %.3f, P = %.3g", RT, estimate, p.value))

# Define group colors (same as plot)
group_colors <- c("ES" = "#3664CC", "MS" = "#32AA10", "LS" = "#FF9933")

# Get axis limits for positioning
x_max <- max(long_df1$log2_on, na.rm = TRUE)
y_max <- max(long_df1$residual, na.rm = TRUE)

# Label placement settings
x_label <- x_max - 1.5   # shift left from edge
line_spacing <- 0.3      # more vertical space
y_positions <- y_max - (0:(nrow(model_stats)-1)) * line_spacing

# Plot base
p1 <- ggplot(long_df1, aes(x = log2_on, y = residual, color = RT)) +
  geom_point(alpha = 0.75, size = 2.5) +
  scale_color_manual(values = group_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  labs(title = "Residual regression analysis for ATAC signal",
       x = "log2(ON)",
       y = "Residuals: log2(ON) - log2(OFF)",
       color = "") +
  xlim(-1.4,4) +
  ylim(-1.4,4) +
  theme_test() +
  theme(axis.text.x = element_text(size=16, color="black"),
        axis.text = element_text(size=16, color="black"),
        axis.title = element_text(size=16, color="black"),
        title = element_text(size=14, color="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Add annotations — left aligned with extra spacing
for (i in seq_len(nrow(model_stats))) {
  p1 <- p1 + annotate("text",
                    x = x_label, 
                    y = y_positions[i],
                    label = model_stats$label[i],
                    hjust = 0.25, 
                    vjust = -4,
                    color = group_colors[model_stats$RT[i]],
                    size = 5)
}

# Print
print(p1)


############################################EdU 4h OCT4 12h OFF residuals regression analysis####################################################################

# Load EdU 4h OCT4 OFF 12h data
df2 <- read.table("EdU_4h_OCT4_12h_data.tsv", header = TRUE, sep = "\t")

# Create long-format data
long_df2 <- bind_rows(
  df2 %>% filter(RT == "ES") %>%
    transmute(RT = "ES",
              log2_off = log2(EdU_4h_OCT4_12h_OFF_ES),
              log2_on  = log2(EdU_4h_OCT4_ON_ES)),
  df2 %>% filter(RT == "MS") %>%
    transmute(RT = "MS",
              log2_off = log2(EdU_4h_OCT4_12h_OFF_MS),
              log2_on  = log2(EdU_4h_OCT4_ON_MS)),
  df2 %>% filter(RT == "LS") %>%
    transmute(RT = "LS",
              log2_off = log2(EdU_4h_OCT4_12h_OFF_LS),
              log2_on  = log2(EdU_4h_OCT4_ON_LS))
) %>%
  mutate(residual = log2_on - log2_off)

# Linear model summaries
model_stats <- long_df2 %>%
  group_by(RT) %>%
  do(tidy(lm(residual ~ log2_on, data = .))) %>%
  filter(term == "log2_on") %>%
  mutate(label = sprintf("%s: Slope = %.3f, P = %.3g", RT, estimate, p.value))

# Define group colors (same as plot)
group_colors <- c("ES" = "#3664CC", "MS" = "#32AA10", "LS" = "#FF9933")

# Get axis limits for positioning
x_max <- max(long_df2$log2_on, na.rm = TRUE)
y_max <- max(long_df2$residual, na.rm = TRUE)

# Label placement settings
x_label <- x_max - 1.5   # shift left from edge
line_spacing <- 0.3      # more vertical space
y_positions <- y_max - (0:(nrow(model_stats)-1)) * line_spacing

# Plot base
p2 <- ggplot(long_df2, aes(x = log2_on, y = residual, color = RT)) +
  geom_point(alpha = 0.75, size = 2.5) +
  scale_color_manual(values = group_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  labs(title = "Residual regression analysis for EdU signal 12h",
       x = "log2(ON)",
       y = "Residuals: log2(ON) - log2(OFF)",
       color = "") +
  xlim(-7.6,6) +
  ylim(-3,3.5) +
  theme_test() +
  theme(axis.text.x = element_text(size=16, color="black"),
        axis.text = element_text(size=16, color="black"),
        axis.title = element_text(size=16, color="black"),
        title = element_text(size=14, color="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Add annotations — left aligned with extra spacing
for (i in seq_len(nrow(model_stats))) {
  p2 <- p2 + annotate("text",
                    x = x_label, 
                    y = y_positions[i],
                    label = model_stats$label[i],
                    hjust = 0.3, 
                    vjust = -0.5,
                    color = group_colors[model_stats$RT[i]],
                    size = 5)
}

# Print
print(p2)

############################################EdU 4h OCT4 18h OFF residuals regression analysis####################################################################

# Load EdU 4h OCT4 OFF 18h data
df3 <- read.table("EdU_4h_OCT4_18h_data.tsv", header = TRUE, sep = "\t")

# Create long-format data
long_df3 <- bind_rows(
  df3 %>% filter(RT == "ES") %>%
    transmute(RT = "ES",
              log2_off = log2(EdU_4h_OCT4_18h_OFF_ES),
              log2_on  = log2(EdU_4h_OCT4_ON_ES)),
  df3 %>% filter(RT == "MS") %>%
    transmute(RT = "MS",
              log2_off = log2(EdU_4h_OCT4_18h_OFF_MS),
              log2_on  = log2(EdU_4h_OCT4_ON_MS)),
  df3 %>% filter(RT == "LS") %>%
    transmute(RT = "LS",
              log2_off = log2(EdU_4h_OCT4_18h_OFF_LS),
              log2_on  = log2(EdU_4h_OCT4_ON_LS))
) %>%
  mutate(residual = log2_on - log2_off)

# Linear model summaries
model_stats <- long_df3 %>%
  group_by(RT) %>%
  do(tidy(lm(residual ~ log2_off, data = .))) %>%
  filter(term == "log2_off") %>%
  mutate(label = sprintf("%s: Slope = %.3f, P = %.3g", RT, estimate, p.value))

# Define group colors (same as plot)
group_colors <- c("ES" = "#3664CC", "MS" = "#32AA10", "LS" = "#FF9933")

# Get axis limits for positioning
x_max <- max(long_df3$log2_off, na.rm = TRUE)
y_max <- max(long_df3$residual, na.rm = TRUE)

# Label placement settings
x_label <- x_max - 1.5   # shift left from edge
line_spacing <- 0.3      # more vertical space
y_positions <- y_max - (0:(nrow(model_stats)-1)) * line_spacing

# Plot base
p3 <- ggplot(long_df3, aes(x = log2_off, y = residual, color = RT)) +
  geom_point(alpha = 0.75, size = 2.5) +
  scale_color_manual(values = group_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2) +
  labs(title = "Residual regression analysis for EdU signal 18h",
       x = "log2(OFF)",
       y = "Residuals: log2(ON) - log2(OFF)",
       color = "") +
  xlim(-7.6,3.5) +
  ylim(-3,3.5) +
  theme_test() +
  theme(axis.text.x = element_text(size=16, color="black"),
        axis.text = element_text(size=16, color="black"),
        axis.title = element_text(size=16, color="black"),
        title = element_text(size=14, color="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

# Add annotations — left aligned with extra spacing
for (i in seq_len(nrow(model_stats))) {
  p3 <- p3 + annotate("text",
                    x = x_label, 
                    y = y_positions[i],
                    label = model_stats$label[i],
                    hjust = 0.75, 
                    vjust = 1,
                    color = group_colors[model_stats$RT[i]],
                    size = 5)
}

# Print
print(p3)


