This package can be used for mediation analysis with repeatedly measured mediators and survival outcomes. 
The jm_mediate() function computes total, direct, and indirect effects over time from a joint model object (fit in JMbayes2). This is a g-computation approach, broadly following [this approach by Zheng & Liu (2021)](https://pubmed.ncbi.nlm.nih.gov/33871871/).

---

```r

# Install from GitHub:
pak::pak("linakramerr/jmed")

library(jmed)

# Example use:
# First fit a joint model:
library(JMbayes2)
data("prothro")
data("prothros")

# longitudinal submodel
lmefit <- lme(pro~ ns(time, 3):treat + ns(time,3), # flexible for other specifications
              random = ~time | id,
              data = prothro,
              control = lmeControl(opt = "optim"))
# survival submodel
coxfit <- coxph(Surv(Time, death)~ treat,
                data = prothros, x = TRUE, na.action = na.omit)

# jointmodel
jointfit <- JMbayes2::jm(coxfit, lmefit, time_var = "time", n_chains =2L,
                         n_iter = 10000L, n_burnin = 1000L, n_thin = 3L,
                         data_Surv = prothros)

summary(jointfit)

# Mediation analysis using jm_mediate():
med_prothro <- jm_mediate(
  jointfit = jointfit,
  ds_surv = prothros,
  ds_long = prothro,
  time_eval = 0:12,
  trt_name = "treat",
  time_var = "time",
  n_mcmc = 500
)

## Visualizing results:
library(ggplot2)
p_colors <- c("Direct" = "#1B9E77", "Indirect" = "#D95F02", "Total" = "#7570B3")

med_plot <- ggplot(data = med_prothro$summary,
                   aes(x = time, y = estimate, color = effect, fill = effect))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15,
              linewidth = 0, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(name = "Years after randomization") +
  scale_y_continuous(name = "Cumulative risk difference")+
  scale_color_manual(values = p_colors) +
  scale_fill_manual(values = p_colors) +
  facet_wrap(~effect)+
  labs(title = "Prednisone vs placebo, N = 488; mediator: prothrombin")+
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 11)
  )

med_plot


```
![Output med_plot.](/man/000003.png)
