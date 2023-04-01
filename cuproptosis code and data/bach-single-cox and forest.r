## 安装成熟的版本
install.packages("ezcox")
## 安装开发版本
devtools::install_github("ShixiangWang/ezcox")

library(ezcox)

input_path <- "yourpath"
out_path <- "yourpath"

## 读取临床信息
STAD_clinial <- read.csv(paste0(input_path,'data_clinical_patient.txt'),sep = "\t")
## 性别改成因子
STAD_clinial$SEX <-  factor(STAD_clinial$SEX)
## 修改生存状态为0，1
STAD_clinial$OS_STATUS <- ifelse(STAD_clinial$OS_STATUS=="1:DECEASED",1,0)

res <-  ezcox(STAD_clinial, 
              covariates = c("SEX","AGE","AJCC_PATHOLOGIC_TUMOR_STAGE"), ## 变量
              time = "OS_MONTHS",
              status = "OS_STATUS",
              return_models = TRUE
              )

library(forestmodel)
mds <- get_models(res)  ## 提取模型
show_models(mds)  ## 森林图
show_models(mds, merge_models = TRUE, drop_controls = TRUE) ## 森林图

## 直接画森林图
show_forest(STAD_clinial, 
            covariates = c("SEX","AGE","AJCC_PATHOLOGIC_TUMOR_STAGE"),
            time = "OS_MONTHS",
            status = "OS_STATUS")
