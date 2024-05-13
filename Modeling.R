
## Install packages
setRepositories(ind = 1:7)

## Load libraries
library(data.table)
library(fst)
library(feather)
library(dplyr)
library(caret)
library(parallel)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(tictoc)
library(car)

## Set working dir.
WORK_DIR <- "D:/2023/서버/IU/bilje/DM/06.13 HW3/"
setwd(WORK_DIR)
getwd()

## Data Load
Hospital1 <- as.data.frame(fread("./GivenMicrobialData/MicrobiomeData_Hospital1.tsv"))
Hospital2 <- as.data.frame(read_feather("./GivenMicrobialData/MicrobiomeData_Hospital2.feather"))
Hospital3 <- as.data.frame(read_fst("./GivenMicrobialData/MicrobiomeData_Hospital3.fst"))
dim(Hospital1)
dim(Hospital2)
dim(Hospital3)

originalData <- rbind(Hospital1, Hospital2, Hospital3)
dim(originalData)

sum(grepl("microbiome", colnames(originalData), ignore.case = TRUE))

# feature 몇 개의 분포 확인
sf <- sample(colnames(originalData), 9)
qqnorm_plots <- lapply(sf, function(column) {
  ggplot(data = originalData, aes(sample = get(column))) +
    stat_qq() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = paste0("QQ-Normal Plot of ", column))
})
grid.arrange(grobs = qqnorm_plots, ncol = 3)

# 각 feature의 p-value 계산
# # anova test
# result <- c()
# for(i in 2:(ncol(originalData))){
#   result[i-1]<-summary(aov(originalData[,i]~originalData$Disease))[[1]][1,5]
# }
# dim(originalData)
# length(result)

#등분산성 확인
var <- c()
for(i in 2:(ncol(originalData))){
  var[i-1] <- leveneTest(originalData[,i] ~ originalData$Disease)$`Pr(>F)`
}
length(var)
#분산이 동일하지 않은 열 수
sum(var < 0.05, na.rm = TRUE) + sum(is.na(var))
#분산이 동일할 수 있는 열 수
sum(var >= 0.05, na.rm = TRUE)

# kruskal.test
result <- c()
for(i in 2:(ncol(originalData))){
  result[i-1] <- kruskal.test(originalData[,i]~originalData$Disease)$p.value
}
dim(originalData)
length(result)

## p값 table 만들기
pValTable <- data.frame(feature = colnames(originalData)[2:ncol(originalData)], pVal = result)
dim(pValTable)

## p 값이 0
pValTable %>%
  filter(pVal == 0) %>% 
  nrow()
## p 값이 0.05이하인 feature selection
pValTable_Rel <- pValTable %>%
  filter(pVal <= 0.05) %>% 
  arrange(pVal)
dim(pValTable_Rel)
head(pValTable_Rel)

selected_features <- pValTable_Rel$feature[1:100]
Data_selected <- originalData[, c(selected_features, "Disease")]
dim(Data_selected)

# 선택한 feature 분포 확인
sf1 <- colnames(Data_selected)[1:5]
qqnorm_plots1 <- lapply(sf1, function(column) {
  ggplot(data = Data_selected, aes(sample = get(column))) +
    stat_qq() +
    geom_abline(intercept = 0, slope = 1, color = "red") +
    labs(title = paste0("QQ-Normal Plot of ", column))
})
grid.arrange(grobs = qqnorm_plots1, ncol = 3)

### Modeling part

# train, test로 나눔
index <- createDataPartition(Data_selected$Disease, p = 0.8, list = FALSE)
trainData <- Data_selected[index, ]
testData <- Data_selected[-index, ]

## Define training control 
trainControl <- trainControl(method = "cv", number = 10)

## For parallel processing
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

## knn
tic("knn Modelling")
model_knn <- train(Disease~., data = trainData, method = "knn", trControl = trainControl)
toc()
model_knn
#예측 결과
predictionResult_knn <- predict(model_knn, newdata = testData)
table_knn <- table(predictionResult_knn, testData$Disease)
print(table_knn)
# Accuracy
accuracy_knn <- sum(diag(table_knn)) / sum(table_knn)
print(accuracy_knn)
accTable <- data.frame(model = c("knn"), accuracy = accuracy_knn)


## rf
tic("rf Modelling")
model_rf <- train(Disease~., data = trainData, method = "rf", trControl = trainControl)
toc()
model_rf
#예측 결과
predictionResult_rf <- predict(model_rf, newdata = testData)
table_rf <- table(predictionResult_rf, testData$Disease)
print(table_rf)
# Accuracy
accuracy_rf <- sum(diag(table_rf)) / sum(table_rf)
print(accuracy_rf)
accTable <- rbind(accTable, data.frame(model = c("rf"), accuracy = accuracy_rf))


## rpart
tic("rpart Modelling")
model_rpart <- train(Disease~., data = trainData, method = "rpart", trControl = trainControl)
toc()
model_rpart
#예측 결과
predictionResult_rpart <- predict(model_rpart, newdata = testData)
table_rpart <- table(predictionResult_rpart, testData$Disease)
print(table_rpart)
# Accuracy
accuracy_rpart <- sum(diag(table_rpart)) / sum(table_rpart)
print(accuracy_rpart)
accTable <- rbind(accTable, data.frame(model = c("rpart"), accuracy = accuracy_rpart))


## svmRadial
tic("svmRadial Modelling")
model_svmRadial <- train(Disease~., data = trainData, method = "svmRadial", trControl = trainControl)
toc()
model_svmRadial
#예측 결과
predictionResult_svmRadial <- predict(model_svmRadial, newdata = testData)
table_svmRadial <- table(predictionResult_svmRadial, testData$Disease)
print(table_svmRadial)
# Accuracy
accuracy_svmRadial <- sum(diag(table_svmRadial)) / sum(table_svmRadial)
print(accuracy_svmRadial)
accTable <- rbind(accTable, data.frame(model = c("svmRadial"), accuracy = accuracy_svmRadial))

## svmLinear3
tic("svmLinear3 Modelling")
model_svmLinear3 <- train(Disease~., data = trainData, method = "svmLinear3", trControl = trainControl)
toc()
model_svmLinear3
#예측 결과
predictionResult_svmLinear3 <- predict(model_svmLinear3, newdata = testData)
table_svmLinear3 <- table(predictionResult_svmLinear3, testData$Disease)
print(table_svmLinear3)
# Accuracy
accuracy_svmLinear3 <- sum(diag(table_svmLinear3)) / sum(table_svmLinear3)
print(accuracy_svmLinear3)
accTable <- rbind(accTable, data.frame(model = c("svmLinear3"), accuracy = accuracy_svmLinear3))


## C5.0
tic("C5.0 Modelling")
model_C5 <- train(Disease~., data = trainData, method = "C5.0", trControl = trainControl)
toc()
model_C5
#예측 결과
predictionResult_C5 <- predict(model_C5, newdata = testData)
table_C5 <- table(predictionResult_C5, testData$Disease)
print(table_C5)
# Accuracy
accuracy_C5 <- sum(diag(table_C5)) / sum(table_C5)
print(accuracy_C5)
accTable <- rbind(accTable, data.frame(model = c("C5.0"), accuracy = accuracy_C5))


## For stop parallel processing
stopCluster(cl)


## test data에 대한 결과
# Accuracy
accTable %>% arrange(-accuracy)


### rf 모델 기반
# 질병별 정확도(Accuracy) 계산
accuracy <- diag(table_rf) / rowSums(table_rf)

# 질병별 민감도(Sensitivity) (True Positive Rate) 계산
sensitivity <- diag(table_rf) / colSums(table_rf)

# 질병별 특이도(Specificity) (True Negative Rate) 계산
total <- sum(table_rf)
true_negatives <- total - rowSums(table_rf) - colSums(table_rf) + diag(table_rf)
specificity <- true_negatives / (total - colSums(table_rf))

# 결과 출력
result <- data.frame(Disease = colnames(table_rf), Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity)
rownames(result) <- 1:nrow(result)
print(result)

########################## Boxplot 그리기 ######################################
# Accuracy table
acc_table <- data.frame(
  knn = model_knn$resample$Accuracy,
  rf = model_rf$resample$Accuracy,
  svmRadial = model_svmRadial$resample$Accuracy,
  svmLinear3 = model_svmLinear3$resample$Accuracy,
  C5.0 = model_C5$resample$Accuracy)
acc_table


# 열 이름을 열로 변환하여 데이터 프레임 생성
acc_long <- pivot_longer(acc_table, cols = everything(), names_to = "model", values_to = "accuracy")
acc_long <- acc_long %>% 
  arrange(model) %>% 
  mutate(model = reorder(model, -accuracy))

# Boxplot 그리기
pastel_colors <- c("#FFB7C5", "#FFE8B7", "#B7FFD5", "#B7E5FF", "#FFB7F5")

acc_plot <- ggplot(acc_long, aes(x = model, y = accuracy, fill = model)) +
  geom_boxplot(color = "black") +
  geom_point(size = 2, color = "black") +
  xlab("Model") +
  ylab("Accuracy") +
  ylim(0.7,0.9) +
  ggtitle('Accuracy for modeling') +
  theme_classic() +
  scale_fill_manual(values = pastel_colors)

acc_plot
ggsave(filename = "/disk7/bilje/kiost/Accuracy.png", plot = acc_plot, width = 10, height = 6)

unique(Data_selected$Disease) %>% 
  as.data.frame()


# 호동행렬 그리기
predictionResult_rf <- predict(model_rf, newdata = testData)
# Create a confusion matrix
confusionMatrix <- confusionMatrix(predictionResult_rf, as.factor(testData$Disease))
# Extract the confusion matrix table
confusionMatrixTable <- as.data.frame(confusionMatrix$table)
# Plot the confusion matrix
confusion_plot <- ggplot(confusionMatrixTable, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Freq), color = "black", size = 3) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Summary of Classfication", x = "", y = "Prediction") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = "black"))

confusion_plot
ggsave(filename = "/disk7/bilje/kiost/Confusion.png", plot = confusion_plot, width = 10, height = 6)

