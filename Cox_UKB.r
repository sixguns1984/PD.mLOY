
#############  UKB  PD  #############

###  Cox  
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
library(openxlsx)
library(survival)
library(survminer)
library(MatchIt)


clinical_male <- fread("c:/Users/bioin/Desktop/ppn/Project/AMP_PD/UKB/clinical_male.csv",data.table = F)

clinical_male <- left_join(clinical_male,clinical_AD)

clinical_male_LOY <- clinical_male[,c("mLRRY","LOY","state","survival_time",
                                      "baseline_state","baseline_ADstate","baseline_age" ,"smoke_now",
                                      # "baseline_NSD","coffee_3YN","bmi", "Diabetes", "hypertension","pesticide",
                                      "PRS",
                                      "coffee_YN",  "tea_YN", 
                                      "alcohol")]
table(clinical_male_LOY$alcohol)

clinical_male_LOY[clinical_male_LOY == ""] <- NA

### filter samples baseline non-AGE
filtered_data <- clinical_male_LOY[!is.na(clinical_male_LOY$baseline_age),]

### filter samples baseline non-PD
filtered_data <- filtered_data[filtered_data$baseline_state == "HC",]


variables <- c("baseline_age" ,"smoke_now",
               "PRS",
               "coffee_YN",  "tea_YN", 
               "alcohol")


filtered_data <- filtered_data %>%
  mutate(survival_time = ifelse(state == "PD", survival_time, 15))
  
filtered_data <- filtered_data %>%
  mutate(status = ifelse(state == "PD", 1, 0))

# for mLRRY
cox_model <- coxph(Surv(survival_time, status) ~ mLRRY+baseline_age+smoke_now+PRS+coffee_YN+
                     tea_YN+alcohol, data = filtered_data)
summary(cox_model)


# for LOY
filtered_data <- filtered_data[filtered_data$baseline_age >= 50,]
filtered_data <- filtered_data[filtered_data$survival_time <= 15,]

cox_model <- coxph(Surv(survival_time, status) ~ LOY+baseline_age+smoke_now+PRS+coffee_YN+
                     tea_YN+alcohol, data = filtered_data)
summary(cox_model)


cox.zph(cox_model)


## plot
plot <- ggsurvplot(
  survfit(Surv(survival_time, status) ~ LOY, data = filtered_data),  # 将LOY作为分组
  data = filtered_data,
  pval = TRUE,         # 显示P值
  risk.table = TRUE,   # 显示风险表
  # conf.int = TRUE,     # 显示置信区间
  xlim = c(0,15),
  ylim = c(0.97,1),
  xlab = "Years after blood sampling",  # X轴标签
  ylab = "Probability of PD-free follow-up",  # Y轴标签
  # legend.title = "",   # 图例标题
  # legend.labs = c("without LOY", "with LOY"),  # 图例标签
  palette = c("black", "red"),  # 设置线条颜色
  ggtheme = theme_minimal()     # 主题样式
)
plot


# average in groups plot
adj1<-ggadjustedcurves(cox_model, 
                       data = filtered_data, 
                       onf.int=T,
                       ylim = c(0.98,1),
                       xlim = c(0,15), 
                       xlab = "Years after blood sampling",  # X轴标签
                       ylab = "Probability of PD-free follow-up",  # Y轴标签
                       pval = TRUE,         # 显示P值
                       risk.table = TRUE,   # 显示风险表
                       # palette = c("black", "red"),  # 设置线条颜色
                       # ggtheme = theme_minimal(),
                       title = paste0("HR = ", HR[1], " 95% CI = ", CI_lower[1], "-", CI_upper[1], " P = ", p_value[1]),
                       variable = "LOY",
                       method = "average",
                       palette = "aaas")
adj1


pdf("c:/Users/bioin/Desktop/ppn/Project/AMP_PD/UKB/cox.pdf")
print(p1)
print(adj1)
print(plot)
dev.off()


