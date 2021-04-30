# global mean and variances : all samples by old and young
df_resu.mat <- data.frame("mean.old"=rowMeans(
  fmat[,str_detect(colnames(fmat),"Old")]))
df_resu.mat$variances.old <- rowVars(
  fmat[,str_detect(colnames(fmat),"Old")]
)
df_resu.mat$mean.young = rowMeans(
  fmat[,str_detect(colnames(fmat),"Young")]
)
df_resu.mat$variances.young <- rowVars(
  fmat[,str_detect(colnames(fmat),"Young")]
)
# visualize mean count
df_resu.mat$id <- rownames(df_resu.mat)
df_resuMelt <- reshape2::melt(df_resu.mat,id=c("id"))
df_resuMelt <- df_resuMelt %>% mutate(metric= case_when(
  str_detect(variable,"mean") ~ "mean",
  str_detect(variable,"variance") ~ "variance"
) )
df_resuMelt <- df_resuMelt %>% mutate(age= case_when(
  str_detect(variable,"old") ~ "Old",
  str_detect(variable,"young") ~ "Young"
) )
ggplot(df_resuMelt, aes(age,value)) +
  geom_boxplot() + 
  geom_jitter(size=.1, alpha=.1, color="lightgray") + 
  scale_y_log10() +
  facet_wrap(~metric) + labs(title="Full dataset" ,
                             subtitle="Mean and Variance (logarithmic)")

ggplot(df_resuMelt, aes(age,value)) +
  geom_boxplot() + 
  geom_jitter(size=.1, alpha=.1, color="lightgray") + 
  scale_y_log10() +
  labs(title="Mean only")
