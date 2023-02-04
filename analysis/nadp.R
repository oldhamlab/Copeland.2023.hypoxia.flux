library(tidyverse)

df <-
  tar_read(nad_data)[["nad"]] |>
  filter(is.na(conc)) |>
  pivot_longer(cols = c("a", "b", "c"), names_to = "rep", values_to = "value") |>
  group_by(across(experiment:replicate)) |>
  wmo::remove_nested_outliers(value, remove = TRUE) |>
  summarize(value = mean(value)) |>
  group_by(across(c(experiment, oxygen:treatment, replicate))) |>
  summarize(ratio = value[nucleotide == "NADPH"] / value[nucleotide == "NADP"]) |>
  group_by(experiment, oxygen, treatment) |>
  summarize(ratio = mean(ratio)) |>
  ungroup() |>
  mutate(
    grand_mean = mean(ratio[!(oxygen == "21%" & treatment == "BAY")]),
    oxygen = factor(oxygen, levels = c("21%", "0.5%")),
    treatment = factor(treatment, levels = c("none", "DMSO", "BAY"))
  ) |>
  group_by(experiment) |>
  mutate(
    exp_mean = mean(ratio[!(oxygen == "21%" & treatment == "BAY")]),
    corr_fac = grand_mean - exp_mean,
    corr_ratio = ratio + corr_fac
  ) |>
  filter(!(experiment == "nadp_2023-01-11" & oxygen == "21%" & treatment == "BAY"))

lmerTest::lmer(ratio ~ oxygen * treatment + (1 | experiment), data = df) |>
  emmeans::emmeans(~ oxygen * treatment) |>
  pairs(simple = "each", combine = TRUE) |>
  broom::tidy()

ggplot(df) +
  aes(
    x = treatment,
    y = corr_ratio,
    fill = oxygen
  ) +
  stat_summary(
    geom = "col",
    fun = "mean",
    position = position_dodge()
  ) +
  stat_summary(
    geom = "errorbar",
    fun.data = "mean_se",
    position = position_dodge(width = 0.9),
    width = 0.2,
    linewidth = 0.25
  ) +
  geom_point(
    aes(group = oxygen, shape = experiment),
    position = position_dodge(width = 0.9)
  )

