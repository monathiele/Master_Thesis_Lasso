
###################################################################################
###################################################################################
############################### 1- DATA PREPARATION ###############################
###################################################################################
###################################################################################

###################################################################################
## Read Data ######################################################################
###################################################################################

library(readxl)

# Read Construction data
cons_data <- as.data.frame(read_excel("./DATA/CONSTRUCTION_SAMPLE.xlsx"))
# Remove lession labels (third colum with lession name etc...)
cons_data$`Cause of Lameness` <- NULL
# Remove Severity (a posteriori information)
cons_data$Severity <- NULL
# Remove N:L (not used, correlation)
cons_data$`N:L` <- NULL
# Move TARGET variable to the right
cons_data <- cons_data[, c(which(colnames(cons_data) != "Classification of Lesion*"), which(colnames(cons_data) == "Classification of Lesion*"))]
# Make sure that Sample variable is in the left and the categorical variable is the second one
cons_data <- cons_data[, c(which(colnames(cons_data) == "Lameness score (5)"), which(colnames(cons_data) != "Lameness score (5)"))]
cons_data <- cons_data[, c(which(colnames(cons_data) == "Sample #"), which(colnames(cons_data) != "Sample #"))]

n1=nrow(cons_data)

# Read Prediction data
pred_data <- as.data.frame(read_excel("./DATA/PREDICTION_SAMPLE_1.xlsx"))
# Remove N:L (not used, correlation)
pred_data$`N:L` <- NULL
# Make sure that Sample variable is in the left and the categorical variable is the second one
pred_data <- pred_data[, c(which(colnames(pred_data) == "Lameness score (5)"), which(colnames(pred_data) != "Lameness score (5)"))]
pred_data <- pred_data[, c(which(colnames(pred_data) == "Sample #"), which(colnames(pred_data) != "Sample #"))]

###################################################################################
## Data Names/Types ###############################################################
###################################################################################

# Lesion column in Construction (not present in Prediction data)
colnames(cons_data)[colnames(cons_data)=="Classification of Lesion*"] <- "Lesion"

# Other columns
change_names <- function(ds) {
  
  colnames(ds)[colnames(ds)=="Sample #"] <- "Sample"
  colnames(ds)[colnames(ds)=="Lameness score (5)"] <- "Lameness"
  colnames(ds)[colnames(ds)=="Rectal Temp (F)"] <- "Rectaltemp"
  colnames(ds)[colnames(ds)=="SAA"] <- "SAA" #Same name
  colnames(ds)[colnames(ds)=="Sub P"] <- "SubP"
  colnames(ds)[colnames(ds)=="Hapto"] <- "Hapto" #Same name
  colnames(ds)[colnames(ds)=="Hair"] <- "Hair" #Same name
  colnames(ds)[colnames(ds)=="RBC (10^12/L 5 - 10)"] <- "RBC"
  colnames(ds)[colnames(ds)=="MCV (fl 40 - 60)"] <- "MCV"
  colnames(ds)[colnames(ds)=="HCT (L/L 024 - 046)"] <- "HCT"
  colnames(ds)[colnames(ds)=="PLT (10^9/L 100 - 800)"] <- "PLT"
  colnames(ds)[colnames(ds)=="MPV (fl 0 - 80)"] <- "MPV"
  colnames(ds)[colnames(ds)=="HGB (g/L 8 - 15)"] <- "HGB"
  colnames(ds)[colnames(ds)=="WBC (10^9/L 4 - 12)"] <- "WBC"
  colnames(ds)[colnames(ds)=="LYM (10^9/L 25 - 75)"] <- "LYM"
  colnames(ds)[colnames(ds)=="MONO (10^9/L 0 - 08)"] <- "MONO"
  colnames(ds)[colnames(ds)=="GRAN (10^9/L 06 - 4)"] <- "GRAN"
  colnames(ds)[colnames(ds)=="Cortisol (mmol/L)"] <- "Cortisol"
  
  return(ds)
}

cons_data <- change_names(cons_data)
pred_data <- change_names(pred_data)
rm(change_names)

# Lameness to factor, other variables must be numeric
# Construction data:
cons_data[["Lameness"]]=as.factor(cons_data[["Lameness"]])
for (icol in 1:ncol(cons_data)) {
  if (colnames(cons_data)[icol] != "Lameness"){cons_data[,icol] <- sapply(cons_data[,icol],as.numeric)}
}
# Prediction data: we must assign the same factor levels to Lameness 
pred_data[["Lameness"]]=factor(as.factor(pred_data[["Lameness"]]),levels=levels(cons_data[["Lameness"]]))
for (icol in 1:ncol(pred_data)) {
  if (colnames(pred_data)[icol] != "Lameness"){pred_data[,icol] <- sapply(pred_data[,icol],as.numeric)}
}
  
rm(icol)

###################################################################################
## Data Cleaning ##################################################################
###################################################################################

#####################
# MPV too high or 0 #
#####################

# Clear Construction data
cons_data=cons_data[cons_data$MPV<=50 | is.na(cons_data$MPV),]
cons_data[(cons_data$MPV == 0 & !is.na(cons_data$MPV)),"MPV"] <- NA

n2=nrow(cons_data)

# For Prediction, construct a flag
pred_data$Warning_MPV <- ifelse(pred_data$MPV>50, "MPV > 50", NA)
pred_data[(pred_data$MPV==0 & !is.na(pred_data$MPV)),"Warning_MPV"] <- "MPV = 0"
pred_data[(pred_data$MPV==0 & !is.na(pred_data$MPV)),"MPV"] <- NA

####################################
# Drop rows with >5 missing values #
####################################

# Clear Construction data
cons_data <- cons_data[rowSums(is.na(cons_data)) <= 5, ]
n3=nrow(cons_data)

# For Prediction, construct a flag
pred_data$Warning_Miss <- ifelse(rowSums(is.na(pred_data[,2:18])) > 5, "Miss > 5", NA)

###################################################################################
## Data Transformation ############################################################
###################################################################################

Transf <- function(ds){
  
  ds[["SAA"]]=sqrt(sqrt(ds[["SAA"]])) # ^(1/4)
  ds[["SubP"]]=sqrt(ds[["SubP"]]) # ^(1/2)
  ds[["Hapto"]]=log(ds[["Hapto"]]) # Log
  ds[["Hair"]]=log(ds[["Hair"]]) # Log
  ds[["MPV"]]=log(ds[["MPV"]]) # Log
  ds[["Cortisol"]]=log(ds[["Cortisol"]]) # Log
  ds$Ilr1 <- (1/sqrt(2))*log(ds$GRAN/ds$LYM)
  ds$Ilr2 <- (1/sqrt(6))*log((ds$GRAN*ds$LYM)/(ds$MONO^2))
  
  ds$GRAN <- NULL
  ds$MONO <- NULL
  ds$LYM <- NULL
  
  ds <- ds[,c(which(colnames(ds) != "Lesion"), which(colnames(ds) == "Lesion"))]
  
  return(ds)
}

cons_data <- Transf(cons_data)
pred_data <- Transf(pred_data)

rm(Transf)

# Move Warnings to the end for pred
pred_data <- pred_data[,c(which(colnames(pred_data) != "Warning_MPV"), which(colnames(pred_data) == "Warning_MPV"))]
pred_data <- pred_data[,c(which(colnames(pred_data) != "Warning_Miss"), which(colnames(pred_data) == "Warning_Miss"))]

###################################################################################
## Missing Imputation #############################################################
###################################################################################

library("mice")

# Since Sample & Lameness are on the left of the dataframe 
# and Lesion (cons) & Warnings (pred) are on the right 
# it suffices to save the positions only for construction (pred shuould be the same)
predictors_num <- numeric()
colClass <- sapply(cons_data,class)
for (icol in 1:ncol(cons_data))
{
  if ((colClass[icol] == "numeric" & colnames(cons_data)[icol] != "Sample" & colnames(cons_data)[icol] != "Lesion")) {predictors_num=c(predictors_num,icol)}
}

rm(icol,colClass)

# Concatenate construction and prediction data (only numeric predictors)
n_cons=nrow(cons_data)
n_pred=nrow(pred_data)
rownames(cons_data) <- 1:n_cons
rownames(pred_data) <- (n_cons+1):(n_cons+n_pred)

full_data_imp <- rbind(cons_data[,predictors_num],pred_data[,predictors_num])
# Create Vector (TRUE/FALSE) for ignoring prediction data to impute
ig=c(rep(FALSE, n_cons), rep(TRUE, n_pred))

# Create imputer
imp <- mice(data=full_data_imp, ignore = ig, method = "norm.nob", m = 1, seed=1234)

# Apply imputer
full_data_imp <- complete(imp)

# Separate imputed dataframes
cons_data_imp=full_data_imp[1:n_cons,]
pred_data_imp=full_data_imp[(n_cons+1):(n_cons+n_pred),]

# Reconstruct full data
cons_data <-cbind(cons_data[,c("Sample","Lameness")],cons_data_imp,cons_data["Lesion"])
pred_data <-cbind(pred_data[,c("Sample","Lameness")],pred_data_imp,pred_data[,c("Warning_MPV","Warning_Miss")])

rm(full_data_imp,cons_data_imp,pred_data_imp,ig,imp,predictors_num,n_cons,n_pred)

###################################################################################
## Final Data Recap ###############################################################
###################################################################################

print(sprintf("Total CONSTRUCTION data rows: %s  /  After MPV cleaning: %s  /  After missing cleaning %s",n1,n2,n3))
print(sprintf("Total PREDICTION data rows: %s",nrow(pred_data)))
rm(n1,n2,n3)


###################################################################################
###################################################################################
################################### 2- LASSO ######################################
###################################################################################
###################################################################################

library(dplyr)
library(InformationValue)
library(glmnet)
library(Metrics)
library(rpart)
library(rlang)

# Define GINI function
GINI <- function(pred,TARGET){return(100*(2*auc(TARGET,pred)-1))}

# Initialize results table
Results=data.frame(N=1:19,Features = c("(Intercept)",colnames(cons_data[,3:17]),"W_Lameness","ALERTS","GINI"))

# Lesion dictionary
Lesions <- vector(mode="list", length=8)
names(Lesions) <- c("Sound", "Foot Rot", "Dermatitis", "P3 Necrosis", "Joint Infection", "Injury", "PLI", "Other")
Lesions[[1]] <- 0 
Lesions[[2]] <- 7
Lesions[[3]] <- c(8,9)
Lesions[[4]] <- 6
Lesions[[5]] <- 1
Lesions[[6]] <- 2
Lesions[[7]] <- 3
Lesions[[8]] <- 4


for (nLesion in 1:8)
{
  # Lesion
  print(sprintf("Lesion %s of 8: - %s",nLesion,names(Lesions[nLesion])))
  
  # Define datsets
  model_data <- cons_data[,1:18]
  apli_data <- pred_data[,1:17]
  
  # Construct TARGET
  model_data$TARGET <- ifelse((cons_data$Lesion %in% Lesions[[nLesion]]), 1, 0)
  model_data$Lesion <- NULL
  
  ###################################################################################
  ## WoE for Lameness ###############################################################
  ###################################################################################
  
  # Construct WOE table for Lameness
  WOE_Lameness=WOETable(model_data$Lameness,model_data$TARGET,valueOfGood=1)
  WOE_Lameness <- WOE_Lameness[,c("CAT","WOE")]
  names(WOE_Lameness)=c("Lameness","W_Lameness")
  
  ######################
  # Aply to dataframes #
  ######################
  
  # Training
  model_data=merge(model_data,WOE_Lameness,by = 'Lameness', all.x = TRUE)
  model_data <- model_data[order(model_data$Sample),]
  rownames(model_data) <- NULL
  model_data$Lameness <- NULL
  # Move TARGET to the end for cons
  model_data <- model_data[, c(which(colnames(model_data) != "TARGET"), which(colnames(model_data) == "TARGET"))]
  
  # Prediction
  apli_data=merge(apli_data,WOE_Lameness,by = 'Lameness', all.x = TRUE)
  apli_data <- apli_data[order(apli_data$Sample),]
  rownames(apli_data) <- NULL
  apli_data$Lameness <- NULL

  rm(WOE_Lameness)
  
  ###################################################################################
  ## Alerts #########################################################################
  ###################################################################################
  
  create_alert <- function (ds_train,ds_test,col,Lesion_Rate,alerts)
  {
    tree_data=ds_train[c("Sample",col,"TARGET")]
    tree_data$TARGET_fac=as.factor(tree_data$TARGET)
    t=rpart(as.formula(paste0("TARGET_fac ~",col)),data=tree_data, method = "class", control = list(maxdepth = 1,cp=-1))
    spl=t$splits[,'index']
    
    alert_left=filter(tree_data, !!sym(col)<=spl)
    alert_right=filter(tree_data, !!sym(col)>spl)

    TOT_left=nrow(alert_left)
    TOT_right=nrow(alert_right)

    Lesioned_left=sum(alert_left$TARGET)
    Lesioned_right=sum(alert_right$TARGET)

    U="NO"
    
    if (nrow(filter(tree_data, !!sym(col)<spl)) <= nrow(filter(tree_data, !!sym(col)>=spl)))
    {
      Acti="<"
      Relative_Lesion_Rate=100*Lesioned_left/TOT_left/Lesion_Rate
      N_Acti=TOT_left
      
      if (Relative_Lesion_Rate > 1.5 | Relative_Lesion_Rate < 2/3 & Relative_Lesion_Rate > 0)
      {
        U="YES"
        new_col=paste0("A_",col)
        ds_train[,new_col] <- 0
        ds_train[(!(is.na(ds_train[,col])) & (ds_train[,col]<spl)),new_col] <- log(Relative_Lesion_Rate)
        ds_test[,new_col] <- 0
        ds_test[(!(is.na(ds_test[,col])) & (ds_test[,col]<spl)),new_col] <- log(Relative_Lesion_Rate)
      }
    } else
    {
      Acti=">="
      Relative_Lesion_Rate=100*Lesioned_right/TOT_right/Lesion_Rate
      N_Acti=TOT_right
      
      if (Relative_Lesion_Rate > 1.5 | Relative_Lesion_Rate < 2/3 & Relative_Lesion_Rate > 0)
      {
        U="YES"
        new_col=paste0("A_",col)
        ds_train[,new_col] <- 0
        ds_train[(!(is.na(ds_train[,col])) & (ds_train[,col]>=spl)),new_col] <- log(Relative_Lesion_Rate)
        ds_test[,new_col] <- 0
        ds_test[(!(is.na(ds_test[,col])) & (ds_test[,col]>=spl)),new_col] <- log(Relative_Lesion_Rate)
      }
    }
    
    P_Act=100*N_Acti/(TOT_right+TOT_left)
    alerta=data.frame("Alerta"=col,"Split"=spl,"Side"=Acti,"N_Act"=N_Acti,"P_Act"=P_Act,"RLR"=Relative_Lesion_Rate,"Use"=U)
    alerts <- rbind(alerts,alerta)
    
    ds_train[,col]<-NULL
    ds_test[,col]<-NULL
    
    return(list(ds_train,ds_test,alerts))
  }
  
  # Define Summarize dataset
  alerts <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(alerts) <- c("Alerta", "Split", "Side", "N_Act","P_Act","LRL","Use")
  
  # Initialize with susceptible features
  ds_train=model_data
  ds_train$W_Lameness <- NULL
  ds_test=apli_data[,1:16]
  
  # Overall Lesion Rate
  Lesion_Rate=100*sum(ds_train$TARGET)/nrow(ds_train)  
  
  # Loop through every possible candidate
  for (col in colnames(model_data[,2:16]))
  {
    if (class(model_data[[col]])=='numeric')
    {
      output=create_alert(ds_train,ds_test,col,Lesion_Rate,alerts)
      ds_train=output[[1]]
      ds_test=output[[2]]
      alerts=output[[3]]
      rm(output)
    }
  }
  
  # Rename alerts summary to keep 
  assign(paste0("alerts_",names(Lesions[nLesion])),alerts)
  
  # Calculate New "alerts" Feature
  model_data$ALERTS <- rowSums (ds_train[,3:ncol(ds_train)])
  apli_data$ALERTS <- rowSums (ds_test[,2:ncol(ds_test)])
  
  # Reorder
  model_data <- model_data[, c(which(colnames(model_data) != "W_Lameness"), which(colnames(model_data) == "W_Lameness"))]
  model_data <- model_data[, c(which(colnames(model_data) != "TARGET"), which(colnames(model_data) == "TARGET"))]
  apli_data <- apli_data[, c(which(colnames(apli_data) != "W_Lameness"), which(colnames(apli_data) == "W_Lameness"))]
  
  rm(col,create_alert,alerts,ds_train,ds_test,Lesion_Rate) 
  
  ###################################################################################
  ## Construct Model ################################################################
  ###################################################################################
  
  # Predictors (Do not use Lameness for "Sound")
  x_linear=model_data[,2:18]
  if (names(Lesions[nLesion])=="Sound") {x_linear=model_data[,2:17]}
  y_linear=model_data[,c("TARGET")]
  
  # Cross-Validate LASSO model
  set.seed(1234)
  cvlasso <- cv.glmnet(x=as.matrix(x_linear), y=y_linear, alpha=1, family="binomial", type.measure="auc", nfolds=5) 
  plot(cvlasso, ylab = paste("AUC - ",names(Lesions[nLesion])))
  
  # Save Predictors
  Coef_aux <- as.matrix(coef(cvlasso, s = "lambda.min"))
  Coef <- data.frame(a=row.names(Coef_aux),b=Coef_aux[,"1"])
  
  ###################################################################################
  ## Predictions ####################################################################
  ###################################################################################
  
  model_data[[names(Lesions[nLesion])]]=predict(cvlasso, s="lambda.min", newx=as.matrix(x_linear), type="response")
  cons_data[[names(Lesions[nLesion])]]=100*predict(cvlasso, s="lambda.min", newx=as.matrix(x_linear), type="response")
  if (names(Lesions[nLesion])=="Sound")
  {
    pred_data[[names(Lesions[nLesion])]]=100*predict(cvlasso, s="lambda.min", newx=as.matrix(apli_data[,2:17]), type="response")
  } else
  {
    pred_data[[names(Lesions[nLesion])]]=100*predict(cvlasso, s="lambda.min", newx=as.matrix(apli_data[,2:18]), type="response")
  }
  
  # Add GINI to Results
  G <- GINI(model_data[[names(Lesions[nLesion])]],model_data$TARGET)
  print(sprintf("GINI: %s",G))
  Coef2 <- data.frame(a="GINI",b=round(G,2))
  Coef <- rbind(Coef,Coef2)
  names(Coef)<-c("Features",names(Lesions[nLesion]))
  Results <- merge(Results, Coef, by="Features", all=TRUE)
  
  rm(x_linear,y_linear,model_data,apli_data,cvlasso,Coef,Coef_aux,Coef2,G)
  
}

rm(nLesion,GINI)

# Final LASSO results
Results <- Results[order(Results$N),]
Results$N <- NULL
rownames(Results) <- NULL
Results[Results == 0] <- NA

# Add Real Lesion Description (Construction)
cons_data$Lesion_desc <- sapply(cons_data$Lesion, function(x) names(Lesions[match(x,Lesions)]))
cons_data$Lesion_desc <- ifelse(cons_data$Lesion == 8 | cons_data$Lesion == 9,"Dermatitis",cons_data$Lesion_desc)


# Top & 2nd-Top Predicted Lesions
#Top
cons_data$model_pred <- colnames(cons_data[,19:26])[apply(cons_data[,19:26],1,which.max)] 
pred_data$model_pred <- colnames(pred_data[,20:27])[apply(pred_data[,20:27],1,which.max)]
# 2nd
v <- apply(cons_data[,19:26], 1,function(x) which(x == sort(x, decreasing = TRUE)[2]))
cons_data$model_pred2 <- colnames(cons_data[,19:26])[v]
rm(v)
v <- apply(pred_data[,20:27], 1,function(x) which(x == sort(x, decreasing = TRUE)[2]))
pred_data$model_pred2 <- colnames(pred_data[,20:27])[v]
rm(v)

# Accuracy
cons_data$Acc <- ifelse(cons_data$Lesion_desc == cons_data$model_pred,1,0)
cons_data$Acc2 <- ifelse(cons_data$Lesion_desc == cons_data$model_pred | cons_data$Lesion_desc == cons_data$model_pred2,1,0)

for (nLesion in 1:8)
{
  
  # Define datset
  mini <- cons_data[cons_data$Lesion %in% Lesions[[nLesion]],c("Sample","Lesion_desc","model_pred","model_pred2","Acc","Acc2")]
  

  cat(sprintf("##################################### %s #####################################\n",names(Lesions[nLesion])))
  cat(sprintf("Total cases: %s / Correct predictions %s / Acc: %s / Extended Acc: %s\n",nrow(mini),sum(mini$Acc),round(100*mean(mini$Acc),2),round(100*mean(mini$Acc2),2)))
  cat("\n")
  cat("First Prediction")
  print(table(mini$Lesion_desc,mini$model_pred))
  cat("\n")
  cat("Second Prediction")
  print(table(mini$Lesion_desc,mini$model_pred2))
  cat("\n")

}

rm(mini,nLesion,Lesions)

###################################################################################
###################################################################################
#################################### 3- LDA #######################################
###################################################################################
###################################################################################

library(MASS)

# Define datasets
model_data=cons_data[,c(3:18)]
apli_data=pred_data[,c(2:17)]

# Define TARGET 
model_data$TARGET <- ifelse(model_data$Lesion == 0, "1",ifelse(model_data$Lesion == 7, "2",ifelse(model_data$Lesion == 8 | model_data$Lesion == 9, "3", "0")))
model_data$Lesion <- NULL

# Train LDA model
lda_fit <- lda(TARGET ~ ., data=model_data, CV=F)

# Make predictions
# Construction data
model_data$pred <- predict(lda_fit,model_data)$class
cons_data$LDA_pred <- ifelse(model_data$pred == "1","Sound",ifelse(model_data$pred == "2","Foot Rot",ifelse(model_data$pred == "3","Dermatitis","All Other")))
# Prediction data
apli_data$pred <- predict(lda_fit,apli_data)$class
pred_data$LDA_pred <- ifelse(apli_data$pred == "1","Sound",ifelse(apli_data$pred == "2","Foot Rot",ifelse(apli_data$pred == "3","Dermatitis","All Other")))

rm(model_data,apli_data,lda_fit)

# Accuracy
cons_data$Lesion_desc_agrup=ifelse(cons_data$Lesion_desc %in% c("Sound","Foot Rot","Dermatitis"),cons_data$Lesion_desc,"All Other")
cons_data$AccLDA <- ifelse(cons_data$LDA_pred == cons_data$Lesion_desc_agrup,1,0)
Lesions_agrup=c("Sound","Foot Rot","Dermatitis","All Other")

for (Lesion in Lesions_agrup)
{
  
  # Define datset
  mini <- cons_data[cons_data$Lesion_desc_agrup == Lesion,c("Sample","Lesion_desc_agrup","LDA_pred","AccLDA")]
  
  
  cat(sprintf("##################################### %s #####################################\n",Lesion))
  cat(sprintf("Total cases: %s / Correct predictions %s / Acc: %s\n",nrow(mini),sum(mini$AccLDA),round(100*mean(mini$AccLDA),2)))
  cat("LDA Prediction")
  print(table(mini$Lesion_desc,mini$LDA_pred))
  cat("\n")
  
}

rm(mini,Lesion,Lesions_agrup)
cons_data$Acc <- NULL
cons_data$Acc2 <- NULL
cons_data$AccLDA <- NULL
cons_data$Lesion_desc_agrup <- NULL





cons_data$Acc <- ifelse(cons_data$Lesion_desc == cons_data$model_pred,1,0)
cons_data$Acc2 <- ifelse(cons_data$Lesion_desc == cons_data$model_pred | cons_data$Lesion_desc == cons_data$model_pred2,1,0)

mean(cons_data$Acc)
mean(cons_data$Acc2)

