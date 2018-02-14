#comap.R
library(reshape2)
library(MASS)

cosine.similarity = function(X, Y)
{
  return(as.matrix(X) %*% t(as.matrix(Y))/
           (norm(as.matrix(X), type = "2")*norm(as.matrix(Y), type = "2")))
}


#read in the data and name the columns
data_colnames = c("MSN","StateCode","Year","Data")
dictionary_colnames = c("MSN","Description","Unit")
data = read.csv2("~/Desktop/COMAP 2018/ProblemCData_seseds.csv", header = TRUE, sep = ",", dec = ".", numerals = "no.loss")
dictionary = read.csv2("~/Desktop/COMAP 2018/ProblemCData_msncodes.csv", header = TRUE, sep = ",")

#type the data
data$MSN = as.character(data$MSN)
data$StateCode = as.factor(data$StateCode)
data$Year = as.numeric(data$Year)
data$Data = as.numeric(data$Data)

dictionary$MSN = as.character(dictionary$MSN)
dictionary$Description = as.character(dictionary$Description)
dictionary$Unit = as.character(dictionary$Unit)

#Subset the long data by state
table(data$StateCode)
AZ.long = subset(data, StateCode == "AZ")
CA.long = subset(data, StateCode == "CA")
NM.long = subset(data, StateCode == "NM")
TX.long = subset(data, StateCode == "TX")

#Remove the state code column since it's all one value for these tables
AZ.long$StateCode = NULL
CA.long$StateCode = NULL
NM.long$StateCode = NULL
TX.long$StateCode = NULL

#Cast the data into wide format - http://seananderson.ca/2013/10/19/reshape.html
AZ.wide = dcast(AZ.long, Year ~ MSN)
CA.wide = dcast(CA.long, Year ~ MSN)
NM.wide = dcast(NM.long, Year ~ MSN)
TX.wide = dcast(TX.long, Year ~ MSN)

#Remove NA values
AZ.wide[is.na(AZ.wide)] = 0
CA.wide[is.na(CA.wide)] = 0
NM.wide[is.na(NM.wide)] = 0
TX.wide[is.na(TX.wide)] = 0

#Principle components analysis on the wide data
#Only years since 1989
AZ.pr = prcomp(AZ.wide[30:50,-1], center = FALSE, scale. = FALSE)
CA.pr = prcomp(CA.wide[30:50,-1], center = FALSE, scale. = FALSE)
NM.pr = prcomp(NM.wide[30:50,-1], center = FALSE, scale. = FALSE)
TX.pr = prcomp(TX.wide[30:50,-1], center = FALSE, scale. = FALSE)


######Arizona#########
summary(AZ.pr) 
#proportion of variance explained by PC1: 0.9927
#plot proportion of variance explained
par(mfrow = c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
par(mai=c(1.02,0.82,0.82,0.42))    
plot(AZ.pr, main = "Proportion of total variance explained", xlab = "Principle components", col = "blueviolet") # print this

#Extract value of the rotated variables
AZ.pc = data.frame(AZ.pr$x)
#add year to the rotated variables data frame
AZ.pc$Year = AZ.wide$Year[30:50]
#Plot Year vs PC1
plot(AZ.pc$Year ~ AZ.pc$PC1, xlab = "Principle component 1", ylab = "Year", col = "springgreen4", pch = 19)
#biplot
plot(AZ.pc$PC2 ~ AZ.pc$PC1, xlab = "Principle component 1", ylab = "Principle component 2", col = "springgreen4", pch = 19)
with(AZ.pc, text(PC2 ~ PC1, labels = Year, pos = 4))

#Extract first principle component loadings
rot = as.matrix(AZ.pr$rotation[,1])
#Check that AZ.pc[,1] == as.matrix(AZ.wide[,-1]) %*% rot
prcomp1 = as.matrix(AZ.wide[30:50,-1]) %*% rot
#Make a data frame out of prcomp1
AZ.prcomp1.df = data.frame(prcomp1)
#Add year column
AZ.prcomp1.df$Year = AZ.wide$Year[30:50]
#Make a linear model to predict prcomp1 from year
AZ.lm.pr = lm(prcomp1 ~ Year, data = AZ.prcomp1.df)
summary(AZ.lm.pr)
#Multiple R-squared:  0.9675,	Adjusted R-squared:  0.9658 
#p-value: 1.33e-15
#Plot Year vs PC1
plot(prcomp1 ~ Year, data = AZ.prcomp1.df, ylab = "Principle component 1", xlab = "Year", 
     col = "springgreen4", pch = 19, xlim = c(1989, 2009), main = "Arizona")
abline(AZ.lm.pr)

#Create newdata df for 2025, 2050
newdata.2025 = data.frame(numeric(1))
colnames(newdata.2025) = "Year"
newdata.2025$Year = 2025
newdata.2050 = data.frame(numeric(1))
colnames(newdata.2050) = "Year"
newdata.2050$Year = 2050

#predict prcomp1 for 2025,2050
newdata.2025$prcomp1.AZ = predict(AZ.lm.pr, newdata = newdata.2025)
newdata.2050$prcomp1.AZ = predict(AZ.lm.pr, newdata = newdata.2050)

#prediction for 2025 - use multiply prcomp1 by pseudoinverse of loadings matrix to get x vector precitions
pred = c( 2025, as.matrix(newdata.2025$prcomp1.AZ[1]) %*% ginv(as.matrix(AZ.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2025 predictions to AZ.wide
colnames(pred) = colnames(AZ.wide)
AZ.wide = rbind(AZ.wide, pred)
#prediction for 2050
pred = c(2050, as.matrix(newdata.2050$prcomp1.AZ[1]) %*% ginv(as.matrix(AZ.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2050 predictions to AZ.wide
colnames(pred) = colnames(AZ.wide)
AZ.wide = rbind(AZ.wide, pred)

#Choose important variables for Energy profile
imp.var = subset(dictionary, grepl("..TCB", MSN))$MSN
dict.imp = subset(dictionary, grepl("..TCB", MSN))
#Energy sources used to compute total energy consumption - p. 119 https://www.eia.gov/state/seds/sep_use/notes/use_technotes.pdf
# Non-Renewable Sources
# Fossil fuels:
# coal (CL)
# net imports of coal coke (U.S. only) - excluded from analysis since data is missing from state level
# natural gas excluding supplemental gaseous fuels (NN)
# petroleum products excluding fuel ethanol blended into motor gasoline (PM)
# Nuclear electric power (NU)
# Renewable Sources
# fuel ethanol minus denaturant (EM)
# geothermal direct use energy and geothermal heat pumps (GE)
# conventional hydroelectric power (HY)
# solar thermal direct use energy and photovoltaic electricity net generation (SO)
# electricity produced by wind (WY)
# wood and wood-derived fuels (WD)
# biomass waste (WS) 
imp.var = c("CLTCB", "NNTCB", "PMTCB", "NUTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "WDTCB", "WSTCB")
#Set color vector to match red colors to non-renewable and green colors to renewable
col.vect = c("tomato", "orangered", "red3", "violetred4",
             "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3", "springgreen4", "darkgreen")
imp.var = data.frame(imp.var, col.vect)
AZ.imp.var = grepl("CLTCB|NNTCB|PMTCB|NUTCB|EMTCB|GETCB|HYTCB|SOTCB|WYTCB|WDTCB|WSTCB", colnames(AZ.wide))
AZ.imp = AZ.wide[,AZ.imp.var]
AZ.imp$sum = AZ.imp$CLTCB + AZ.imp$EMTCB + AZ.imp$GETCB + AZ.imp$HYTCB + AZ.imp$NNTCB + AZ.imp$PMTCB + AZ.imp$SOTCB + AZ.imp$WYTCB
#Plot Arizona Energy Profile over time
par(mfrow = c(1,3))
barplot(as.numeric(AZ.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
    col = c("tomato", "orangered", "red3",
            "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
    main = "2009", ylim = c(0, 1100000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(AZ.imp[51,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
    col = c("tomato", "orangered", "red3",
            "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
    main = "2025", ylim = c(0, 1100000))
barplot(as.numeric(AZ.imp[52,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
    col = c("tomato", "orangered", "red3",
            "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
    main = "2050", ylim = c(0, 1100000))
#Create a Legend
par(mfrow = c(1,1))
plot(1, type="n", axes=FALSE, xlab="", ylab="")
legend("center", inset=.00001, title="Energy Source",
       c("Coal", "Natural Gas", "Petroleum Products", "Fuel Ethanol ", "Geothermal", "Hydroelectric", "Solar", "Wind"), fill=c("tomato", "orangered", "red3",
       "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"), horiz=FALSE, cex=0.8)

#IB. Develop a model to characterize how the energy profile of each of the four states has evolved from 1960 - 2009
AZ.wide$sum.clean = (AZ.wide$EMTCB + AZ.wide$GETCB + AZ.wide$HYTCB + AZ.wide$SOTCB + AZ.wide$WYTCB)
lm.AZ.clean = lm(sum.clean ~ Year + TPOPP, data = AZ.wide[1:50,])
#tried GDP and total energy consumption per real dollar of GDP (TETGR), not significant
summary(lm.AZ.clean)
#Multiple R-squared:  0.4971,	Adjusted R-squared:  0.4757 
#p-value: 9.65e-08
par(mfrow= c(1,2))
plot(AZ.wide$sum.clean[1:50] ~ AZ.wide$Year[1:50], xlab = "Year", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)
plot(AZ.wide$sum.clean[1:50] ~ AZ.wide$TPOPP[1:50], xlab = "Total Population", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)
#Adding climate data to the model such as yearly precipitation totals and sunny days might improve the Rsq since AZ is highly dependent on hydroelectric


#IB try again - stepwise regression of %clean energy on TPOPP, GDPRX, TETGR, Year, and any variable that ends in D
dict.IB = subset(dictionary, grepl("....D|TPOPP|GDPRX|TETGR", MSN))
AZ.wide$percent.clean = (AZ.wide$EMTCB + AZ.wide$GETCB + AZ.wide$HYTCB + AZ.wide$SOTCB + AZ.wide$WYTCB)/
  (AZ.imp$CLTCB + AZ.imp$EMTCB + AZ.imp$GETCB + AZ.imp$HYTCB + AZ.imp$NNTCB + AZ.imp$PMTCB + AZ.imp$SOTCB + AZ.imp$WYTCB)
AZ.imp.var.IB = grepl("....D|TPOPP|GDPRX|TETGR|Year|percent.clean", colnames(AZ.wide))
AZ.IB = AZ.wide[1:50,AZ.imp.var.IB]
AZ.IB.lm.full = lm(percent.clean ~ ., data = AZ.IB)
summary(AZ.IB.lm.full)
alias(AZ.IB.lm.full)
AZ.IB.step = stepAIC(AZ.IB.lm.full, trace = 0)
summary(AZ.IB.step)
par(mfrow = c(1,1))
plot(predict.lm(AZ.IB.step) ~ AZ.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Arizona", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)






######California#########
summary(CA.pr) 
#proportion of variance explained by PC1: 0.9978
#plot proportion of variance explained
par(mfrow = c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
par(mai=c(1.02,0.82,0.82,0.42))    
plot(CA.pr, main = "Proportion of total variance explained", xlab = "Principle components", col = "blueviolet") # print this

#Extract value of the rotated variables
CA.pc = data.frame(CA.pr$x)
#add year to the rotated variables data frame
CA.pc$Year = CA.wide$Year[30:50]
#Plot Year vs PC1
plot(CA.pc$Year ~ CA.pc$PC1, xlab = "Principle component 1", ylab = "Year", col = "springgreen4", pch = 19)
#biplot
plot(CA.pc$PC2 ~ CA.pc$PC1, xlab = "Principle component 1", ylab = "Principle component 2", col = "springgreen4", pch = 19)
with(CA.pc, text(PC2 ~ PC1, labels = Year, pos = 4))

#Extract first principle component loadings
rot = as.matrix(CA.pr$rotation[,1])
#Check that CA.pc[,1] == as.matrix(CA.wide[,-1]) %*% rot
prcomp1 = as.matrix(CA.wide[30:50,-1]) %*% rot
#Make a data frame out of prcomp1
CA.prcomp1.df = data.frame(prcomp1)
#Add year column
CA.prcomp1.df$Year = CA.wide$Year[30:50]
#Make a linear model to predict prcomp1 from year
CA.lm.pr = lm(prcomp1 ~ Year, data = CA.prcomp1.df)
summary(CA.lm.pr)
#Multiple R-squared:  0.7843,	Adjusted R-squared:  0.773 
#p-value: 9.436e-08
#Plot Year vs PC1
plot(prcomp1 ~ Year, data = CA.prcomp1.df, ylab = "Principle component 1", xlab = "Year", 
     col = "springgreen4", pch = 19, xlim = c(1989, 2009), main = "California")
abline(CA.lm.pr)

#predict prcomp1 for 2025,2050
newdata.2025$prcomp1.CA = predict(CA.lm.pr, newdata = newdata.2025)
newdata.2050$prcomp1.CA = predict(CA.lm.pr, newdata = newdata.2050)

#prediction for 2025 - use multiply prcomp1 by pseudoinverse of loadings matrix to get x vector precitions
pred = c( 2025, as.matrix(newdata.2025$prcomp1.CA[1]) %*% ginv(as.matrix(CA.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2025 predictions to AZ.wide
colnames(pred) = colnames(CA.wide)
CA.wide = rbind(CA.wide, pred)
#prediction for 2050
pred = c(2050, as.matrix(newdata.2050$prcomp1.CA[1]) %*% ginv(as.matrix(CA.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2050 predictions to AZ.wide
colnames(pred) = colnames(CA.wide)
CA.wide = rbind(CA.wide, pred)


imp.var = c("CLTCB", "NNTCB", "PMTCB", "NUTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "WDTCB", "WSTCB")
#Set color vector to match red colors to non-renewable and green colors to renewable
col.vect = c("tomato", "orangered", "red3", "violetred4",
             "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3", "springgreen4", "darkgreen")
imp.var = data.frame(imp.var, col.vect)
CA.imp.var = grepl("CLTCB|NNTCB|PMTCB|NUTCB|EMTCB|GETCB|HYTCB|SOTCB|WYTCB|WDTCB|WSTCB", colnames(CA.wide))
CA.imp = CA.wide[,CA.imp.var]
CA.imp$sum = CA.imp$CLTCB + CA.imp$EMTCB + CA.imp$GETCB + CA.imp$HYTCB + CA.imp$NNTCB + CA.imp$PMTCB + CA.imp$SOTCB + CA.imp$WYTCB
#Plot California Energy Profile over time
par(mfrow = c(1,3))
barplot(as.numeric(CA.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2009", ylim = c(0, 5000000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(CA.imp[51,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2025", ylim = c(0, 5000000))
barplot(as.numeric(CA.imp[52,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2050", ylim = c(0, 5000000))

#IB. Develop a model to characterize how the energy profile of each of the four states has evolved from 1960 - 2009
CA.wide$sum.clean = (CA.wide$EMTCB + CA.wide$GETCB + CA.wide$HYTCB + CA.wide$SOTCB + CA.wide$WYTCB)
lm.CA.clean = lm(sum.clean ~ Year + TPOPP, data = CA.wide[1:50,])
#tried GDP and total energy consumption per real dollar of GDP (TETGR), not significant
summary(lm.CA.clean)
#Multiple R-squared:  0.6313,	Adjusted R-squared:  0.6156 
#p-value: p-value: 6.543e-11
par(mfrow= c(1,2))
plot(CA.wide$sum.clean[1:50] ~ CA.wide$Year[1:50], xlab = "Year", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)
plot(CA.wide$sum.clean[1:50] ~ CA.wide$TPOPP[1:50], xlab = "Total Population", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)

#IB try again - stepwise regression of %clean energy on TPOPP, GDPRX, TETGR, Year, and any variable that ends in D
dict.IB = subset(dictionary, grepl("....D|TPOPP|GDPRX|TETGR", MSN))
CA.wide$percent.clean = (CA.wide$EMTCB + CA.wide$GETCB + CA.wide$HYTCB + CA.wide$SOTCB + CA.wide$WYTCB)/
  (CA.imp$CLTCB + CA.imp$EMTCB + CA.imp$GETCB + CA.imp$HYTCB + CA.imp$NNTCB + CA.imp$PMTCB + CA.imp$SOTCB + CA.imp$WYTCB)
CA.imp.var.IB = grepl("....D|TPOPP|GDPRX|TETGR|Year|percent.clean", colnames(CA.wide))
CA.IB = CA.wide[1:50,CA.imp.var.IB]
CA.IB.lm.full = lm(percent.clean ~ ., data = CA.IB)
summary(CA.IB.lm.full)
CA.IB.step = stepAIC(CA.IB.lm.full, trace = 0)
summary(CA.IB.step)
par(mfrow = c(1,1))
plot(predict.lm(CA.IB.step) ~ CA.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "California", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)


######New Mexico#########
summary(NM.pr) 
#proportion of variance explained by PC1: 0.9970
#plot proportion of variance explained
par(mfrow = c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
par(mai=c(1.02,0.82,0.82,0.42))    
plot(NM.pr, main = "Proportion of total variance explained", xlab = "Principle components", col = "blueviolet") # print this

#Extract value of the rotated variables
NM.pc = data.frame(NM.pr$x)
#add year to the rotated variables data frame
NM.pc$Year = NM.wide$Year[30:50]
#Plot Year vs PC1
plot(NM.pc$Year ~ NM.pc$PC1, xlab = "Principle component 1", ylab = "Year", col = "springgreen4", pch = 19)
#biplot
plot(NM.pc$PC2 ~ NM.pc$PC1, xlab = "Principle component 1", ylab = "Principle component 2", col = "springgreen4", pch = 19)
with(NM.pc, text(PC2 ~ PC1, labels = Year, pos = 4))

#Extract first principle component loadings
rot = as.matrix(NM.pr$rotation[,1])
#Check that NM.pc[,1] == as.matrix(NM.wide[,-1]) %*% rot
prcomp1 = as.matrix(NM.wide[30:50,-1]) %*% rot
#Make a data frame out of prcomp1
NM.prcomp1.df = data.frame(prcomp1)
#Add year column
NM.prcomp1.df$Year = NM.wide$Year[30:50]
#Make a linear model to predict prcomp1 from year
NM.lm.pr = lm(prcomp1 ~ Year, data = NM.prcomp1.df)
summary(NM.lm.pr)
#Multiple R-squared:  0.4299,	Adjusted R-squared:  0.3999 
#p-value: 0.001251
#Plot Year vs PC1
plot(prcomp1 ~ Year, data = NM.prcomp1.df, ylab = "Principle component 1", xlab = "Year", 
     col = "springgreen4", pch = 19, xlim = c(1989, 2009), main = "New Mexico")
abline(NM.lm.pr)

#predict prcomp1 for 2025,2050
newdata.2025$prcomp1.NM = predict(NM.lm.pr, newdata = newdata.2025)
newdata.2050$prcomp1.NM = predict(NM.lm.pr, newdata = newdata.2050)

#prediction for 2025 - use multiply prcomp1 by pseudoinverse of loadings matrix to get x vector precitions
pred = c( 2025, as.matrix(newdata.2025$prcomp1.NM[1]) %*% ginv(as.matrix(NM.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2025 predictions to AZ.wide
colnames(pred) = colnames(NM.wide)
NM.wide = rbind(NM.wide, pred)
#prediction for 2050
pred = c(2050, as.matrix(newdata.2050$prcomp1.NM[1]) %*% ginv(as.matrix(NM.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2050 predictions to AZ.wide
colnames(pred) = colnames(NM.wide)
NM.wide = rbind(NM.wide, pred)


imp.var = c("CLTCB", "NNTCB", "PMTCB", "NUTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "WDTCB", "WSTCB")
#Set color vector to match red colors to non-renewable and green colors to renewable
col.vect = c("tomato", "orangered", "red3", "violetred4",
             "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3", "springgreen4", "darkgreen")
imp.var = data.frame(imp.var, col.vect)
NM.imp.var = grepl("CLTCB|NNTCB|PMTCB|NUTCB|EMTCB|GETCB|HYTCB|SOTCB|WYTCB|WDTCB|WSTCB", colnames(NM.wide))
NM.imp = NM.wide[,NM.imp.var]
NM.imp$sum = NM.imp$CLTCB + NM.imp$EMTCB + NM.imp$GETCB + NM.imp$HYTCB + NM.imp$NNTCB + NM.imp$PMTCB + NM.imp$SOTCB + NM.imp$WYTCB
#Plot New Mexico Energy Profile over time
par(mfrow = c(1,3))
barplot(as.numeric(NM.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2009", ylim = c(0, 450000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(NM.imp[51,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2025", ylim = c(0, 450000))
barplot(as.numeric(NM.imp[52,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2050", ylim = c(0, 450000))

#IB. Develop a model to characterize how the energy profile of each of the four states has evolved from 1960 - 2009
NM.wide$sum.clean = (NM.wide$EMTCB + NM.wide$GETCB + NM.wide$HYTCB + NM.wide$SOTCB + NM.wide$WYTCB)
lm.NM.clean = lm(sum.clean ~ Year + TPOPP + I(TPOPP^2), data = NM.wide[1:50,])
#tried GDP and total energy consumption per real dollar of GDP (TETGR), not significant
summary(lm.NM.clean)
#Multiple R-squared:  0.8047,	Adjusted R-squared:  0.7919 
#p-value: p-value: 2.416e-16
par(mfrow= c(1,2))
plot(NM.wide$sum.clean[1:50] ~ NM.wide$Year[1:50], xlab = "Year", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)
plot(NM.wide$sum.clean[1:50] ~ NM.wide$TPOPP[1:50], xlab = "Total Population", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)
#Rapid change in the past few years so we added a second order term to the model for growth

#IB try again - stepwise regression of %clean energy on TPOPP, GDPRX, TETGR, Year, and any variable that ends in D
dict.IB = subset(dictionary, grepl("....D|TPOPP|GDPRX|TETGR", MSN))
NM.wide$percent.clean = (NM.wide$EMTCB + NM.wide$GETCB + NM.wide$HYTCB + NM.wide$SOTCB + NM.wide$WYTCB)/
  (NM.imp$CLTCB + NM.imp$EMTCB + NM.imp$GETCB + NM.imp$HYTCB + NM.imp$NNTCB + NM.imp$PMTCB + NM.imp$SOTCB + NM.imp$WYTCB)
NM.imp.var.IB = grepl("....D|TPOPP|GDPRX|TETGR|Year|percent.clean", colnames(NM.wide))
NM.IB = NM.wide[1:50,NM.imp.var.IB]
NM.IB.lm.full = lm(percent.clean ~ ., data = NM.IB)
summary(NM.IB.lm.full)
NM.IB.step = stepAIC(NM.IB.lm.full, trace = 0)
summary(NM.IB.step)
par(mfrow = c(1,1))
plot(predict.lm(NM.IB.step) ~ NM.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "New Mexico", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)



######Texax#########
summary(TX.pr) 
#proportion of variance explained by PC1: 0.9963
#plot proportion of variance explained
par(mfrow = c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
par(mai=c(1.02,0.82,0.82,0.42))    
plot(TX.pr, main = "Proportion of total variance explained", xlab = "Principle components", col = "blueviolet") # print this

#Extract value of the rotated variables
TX.pc = data.frame(TX.pr$x)
#add year to the rotated variables data frame
TX.pc$Year = TX.wide$Year[30:50]
#Plot Year vs PC1
plot(TX.pc$Year ~ TX.pc$PC1, xlab = "Principle component 1", ylab = "Year", col = "springgreen4", pch = 19)
#biplot
plot(TX.pc$PC2 ~ TX.pc$PC1, xlab = "Principle component 1", ylab = "Principle component 2", col = "springgreen4", pch = 19)
with(TX.pc, text(PC2 ~ PC1, labels = Year, pos = 4))

#Extract first principle component loadings
rot = as.matrix(TX.pr$rotation[,1])
#Check that TX.pc[,1] == as.matrix(TX.wide[,-1]) %*% rot
prcomp1 = as.matrix(TX.wide[30:50,-1]) %*% rot
#Make a data frame out of prcomp1
TX.prcomp1.df = data.frame(prcomp1)
#Add year column
TX.prcomp1.df$Year = TX.wide$Year[30:50]
#Make a linear model to predict prcomp1 from year
TX.lm.pr = lm(prcomp1 ~ Year, data = TX.prcomp1.df)
summary(TX.lm.pr)
#Multiple R-squared:  0.5563,	Adjusted R-squared:  0.5563 
#p-value: 0.0001038
#Plot Year vs PC1
plot(prcomp1 ~ Year, data = TX.prcomp1.df, ylab = "Principle component 1", xlab = "Year", 
     col = "springgreen4", pch = 19, xlim = c(1989, 2009), main = "Texas")
abline(TX.lm.pr)

#predict prcomp1 for 2025,2050
newdata.2025$prcomp1.TX = predict(TX.lm.pr, newdata = newdata.2025)
newdata.2050$prcomp1.TX = predict(TX.lm.pr, newdata = newdata.2050)

#prediction for 2025 - use multiply prcomp1 by pseudoinverse of loadings matrix to get x vector precitions
pred = c( 2025, as.matrix(newdata.2025$prcomp1.TX[1]) %*% ginv(as.matrix(TX.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2025 predictions to AZ.wide
colnames(pred) = colnames(TX.wide)
TX.wide = rbind(TX.wide, pred)
#prediction for 2050
pred = c(2050, as.matrix(newdata.2050$prcomp1.TX[1]) %*% ginv(as.matrix(TX.pr$rotation)[,1]))
pred = data.frame(t(pred))
#rbind 2050 predictions to AZ.wide
colnames(pred) = colnames(TX.wide)
TX.wide = rbind(TX.wide, pred)


imp.var = c("CLTCB", "NNTCB", "PMTCB", "NUTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "WDTCB", "WSTCB")
#Set color vector to match red colors to non-renewable and green colors to renewable
col.vect = c("tomato", "orangered", "red3", "violetred4",
             "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3", "springgreen4", "darkgreen")
imp.var = data.frame(imp.var, col.vect)
TX.imp.var = grepl("CLTCB|NNTCB|PMTCB|NUTCB|EMTCB|GETCB|HYTCB|SOTCB|WYTCB|WDTCB|WSTCB", colnames(TX.wide))
TX.imp = TX.wide[,TX.imp.var]
TX.imp$sum = TX.imp$CLTCB + TX.imp$EMTCB + TX.imp$GETCB + TX.imp$HYTCB + TX.imp$NNTCB + TX.imp$PMTCB + TX.imp$SOTCB + TX.imp$WYTCB
#Plot Texas Energy Profile over time
par(mfrow = c(1,3))
barplot(as.numeric(TX.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2009", ylim = c(0, 6500000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(TX.imp[51,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2025", ylim = c(0, 6500000))
barplot(as.numeric(TX.imp[52,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2050", ylim = c(0, 6500000))

#IB. Develop a model to characterize how the energy profile of each of the four states has evolved from 1960 - 2009
TX.wide$sum.clean = (TX.wide$EMTCB + TX.wide$GETCB + TX.wide$HYTCB + TX.wide$SOTCB + TX.wide$WYTCB)
lm.TX.clean = lm(sum.clean ~ Year + TPOPP + I(TPOPP^2), data = TX.wide[1:50,])
#tried GDP and total energy consumption per real dollar of GDP (TETGR), not significant
summary(lm.TX.clean)
#Multiple R-squared:  0.7012,	Adjusted R-squared:  0.6817 
#p-value: 4.001e-12
par(mfrow= c(1,2))
plot(TX.wide$sum.clean[1:50] ~ TX.wide$Year[1:50], xlab = "Year", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)
plot(TX.wide$sum.clean[1:50] ~ TX.wide$TPOPP[1:50], xlab = "Total Population", ylab = "Consumed Renewable Energy", 
     col = "springgreen4", pch = 19)
#Rapid change in the past few years so we added a second order term to the model for growth

#IB try again - stepwise regression of %clean energy on TPOPP, GDPRX, TETGR, Year, and any variable that ends in D
dict.IB = subset(dictionary, grepl("....D|TPOPP|GDPRX|TETGR", MSN))
TX.wide$percent.clean = (TX.wide$EMTCB + TX.wide$GETCB + TX.wide$HYTCB + TX.wide$SOTCB + TX.wide$WYTCB)/
  (TX.imp$CLTCB + TX.imp$EMTCB + TX.imp$GETCB + TX.imp$HYTCB + TX.imp$NNTCB + TX.imp$PMTCB + TX.imp$SOTCB + TX.imp$WYTCB)
TX.imp.var.IB = grepl("....D|TPOPP|GDPRX|TETGR|Year|percent.clean", colnames(TX.wide))
TX.IB = TX.wide[1:50,TX.imp.var.IB]
TX.IB.lm.full = lm(percent.clean ~ ., data = TX.IB)
summary(TX.IB.lm.full)
TX.IB.step = stepAIC(TX.IB.lm.full, trace = 0)
summary(TX.IB.step)
par(mfrow = c(1,1))
plot(predict.lm(TX.IB.step) ~ TX.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Texas", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)






############PLOTS#############
######Plot 2009 profiles against each other#######
par(mfrow = c(1,4))
barplot(as.numeric(AZ.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "Arizona", ylim = c(0, 6500000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(CA.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "California", ylim = c(0, 6500000))
barplot(as.numeric(NM.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "New Mexico", ylim = c(0, 6500000))
barplot(as.numeric(TX.imp[50,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "Texas", ylim = c(0, 6500000))

percent.renew.AZ = AZ.wide$sum.clean/AZ.imp$sum
percent.renew.CA = CA.wide$sum.clean/CA.imp$sum
percent.renew.NM = NM.wide$sum.clean/NM.imp$sum
percent.renew.TX = AZ.wide$sum.clean/TX.imp$sum


percent.renew.AZ[50]*100
percent.renew.CA[50]*100
percent.renew.NM[50]*100
percent.renew.TX[50]*100



#plot scree plots for PCA for four states
par(mfrow = c(1,4))
barplot(summary(AZ.pr)$importance[2,1:10], main = "Arizona", col = "blueviolet", ylim = c(0,1), ylab = "Proportion of variance explained", names.arg = NULL)
barplot(summary(CA.pr)$importance[2,1:10], main = "California", col = "blueviolet", ylim = c(0,1), names.arg = NULL)
barplot(summary(NM.pr)$importance[2,1:10], main = "New Mexico", col = "blueviolet", ylim = c(0,1), names.arg = NULL)
barplot(summary(TX.pr)$importance[2,1:10], main = "Texas", col = "blueviolet", ylim = c(0,1), names.arg = NULL)


AZ.pc1.propvar = summary(AZ.pr)$importance[2,1]
CA.pc1.propvar = summary(CA.pr)$importance[2,1]
NM.pc1.propvar = summary(NM.pr)$importance[2,1]
TX.pc1.propvar = summary(TX.pr)$importance[2,1]
propvar = c(AZ.pc1.propvar, CA.pc1.propvar, NM.pc1.propvar, TX.pc1.propvar)
names(propvar) = c("Arizona", "California", "New Mexico", "Texas")



#plot predicted IB vs actual for four states
par(mfrow = c(2,2))
plot(predict.lm(AZ.IB.step) ~ AZ.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Arizona", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)
plot(predict.lm(CA.IB.step) ~ CA.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "California", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)
plot(predict.lm(NM.IB.step) ~ NM.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "New Mexico", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)
plot(predict.lm(TX.IB.step) ~ TX.IB$percent.clean, xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Texas", xlim = c(0,0.3), ylim = c(0,0.3), 
     col = "springgreen4", pch = 19)
abline(0, 1)

#make table of RSq and adjRsr and Pvalue
IB = data.frame(Rsq = numeric(4), adjRsq = numeric(4), pvalue = numeric(4))
rownames(IB) = c("Arizona", "California", "New Mexico", "Texas")
IB[1,] = c(summary(AZ.IB.step)$r.squared, summary(AZ.IB.step)$adj.r.squared, 0.0007217)
IB[2,] = c(summary(CA.IB.step)$r.squared, summary(CA.IB.step)$adj.r.squared, 0.003428)
IB[3,] = c(summary(NM.IB.step)$r.squared, summary(NM.IB.step)$adj.r.squared, 0.003428)
IB[4,] = c(summary(TX.IB.step)$r.squared, summary(TX.IB.step)$adj.r.squared, 3.112e-10)
IB





###### IC ########
opt.percent = data.frame(matrix(NA, nrow = 5, ncol = 8))
colnames(opt.percent) = colnames(AZ.imp)[1:8]
rownames(opt.percent) = c("Default", "Arizona", "California", "New Mexico", "Texas")
opt.percent["Default",] = c(0,.2,.2,.2,0,0,.2,.2)
opt.percent["Arizona",] = c(0,.2,.2,.2,0,0,.2,.2)
opt.percent["California",] = c(0,.2,.2,.2,0,0,.2,.2)
opt.percent["New Mexico",] = c(0,.2,.2,.2,0,0,.2,.2)
opt.percent["Texas",] = c(0,.2,.2,.2,0,0,.2,.2)

opt.BTU = data.frame(matrix(NA, nrow = 4, ncol = 8))
colnames(opt.BTU) = colnames(AZ.imp)[1:8]
rownames(opt.BTU) = c("Arizona", "California", "New Mexico", "Texas")
opt.BTU["Arizona",] = AZ.imp$sum[50]*opt.percent["Arizona",]
opt.BTU["California",] = CA.imp$sum[50]*opt.percent["California",]
opt.BTU["New Mexico",] = NM.imp$sum[50]*opt.percent["New Mexico",]
opt.BTU["Texas",] = TX.imp$sum[50]*opt.percent["Texas",]

library(tcR)
AZ.perform = cosine.similarity(opt.BTU["Arizona",], AZ.imp[50,1:8])
CA.perform = cosine.similarity(opt.BTU["California",], CA.imp[50,1:8])
NM.perform = cosine.similarity(opt.BTU["New Mexico",], NM.imp[50,1:8])
TX.perform = cosine.similarity(opt.BTU["Texas",], TX.imp[50,1:8])

perform = c(AZ.perform, CA.perform, NM.perform, TX.perform)
names(perform) = c("Arizona", "California", "New Mexico", "Texas")


AZ.opt.BTU = data.frame(matrix(NA, nrow = 50, ncol = 8))
colnames(AZ.opt.BTU) = colnames(AZ.imp)[1:8]
rownames(AZ.opt.BTU) = 1960:2009
for(i in 1:nrow(AZ.opt.BTU))
{
  AZ.opt.BTU[i,] = AZ.imp$sum[i]*opt.percent["Arizona",]
}
AZ.opt.BTU$sim.score = 0
for(i in 1:nrow(AZ.opt.BTU))
{
  AZ.opt.BTU$sim.score[i] = cosine.similarity(AZ.opt.BTU[i,1:8], AZ.imp[i,1:8])
}

CA.opt.BTU = data.frame(matrix(NA, nrow = 50, ncol = 8))
colnames(CA.opt.BTU) = colnames(CA.imp)[1:8]
rownames(CA.opt.BTU) = 1960:2009
for(i in 1:nrow(CA.opt.BTU))
{
  CA.opt.BTU[i,] = CA.imp$sum[i]*opt.percent["California",]
}
CA.opt.BTU$sim.score = 0
for(i in 1:nrow(CA.opt.BTU))
{
  CA.opt.BTU$sim.score[i] = cosine.similarity(CA.opt.BTU[i,1:8], CA.imp[i,1:8])
}


NM.opt.BTU = data.frame(matrix(NA, nrow = 50, ncol = 8))
colnames(NM.opt.BTU) = colnames(NM.imp)[1:8]
rownames(NM.opt.BTU) = 1960:2009
for(i in 1:nrow(NM.opt.BTU))
{
  NM.opt.BTU[i,] = NM.imp$sum[i]*opt.percent["New Mexico",]
}
NM.opt.BTU$sim.score = 0
for(i in 1:nrow(NM.opt.BTU))
{
  NM.opt.BTU$sim.score[i] = cosine.similarity(NM.opt.BTU[i,1:8], NM.imp[i,1:8])
}


TX.opt.BTU = data.frame(matrix(NA, nrow = 50, ncol = 8))
colnames(TX.opt.BTU) = colnames(TX.imp)[1:8]
rownames(TX.opt.BTU) = 1960:2009
for(i in 1:nrow(TX.opt.BTU))
{
  TX.opt.BTU[i,] = TX.imp$sum[i]*opt.percent["Texas",]
}
TX.opt.BTU$sim.score = 0
for(i in 1:nrow(TX.opt.BTU))
{
  TX.opt.BTU$sim.score[i] = cosine.similarity(TX.opt.BTU[i,1:8], TX.imp[i,1:8])
}




target.AZ.2025 = 0.2
target.AZ.2050 = 0.4
target.CA.2025 = 0.2
target.CA.2050 = 0.4
target.NM.2025 = 0.2
target.NM.2050 = 0.4
target.TX.2025 = 0.2
target.TX.2050 = 0.4

get.target = function(target_cosine, profile_percent, year_consump)
{
  best_profile = year_consump*profile_percent
  for(i in 1:8)
  {
    if(best_profile[1,i] == 0)
      best_profile[1,i] = 1
  }
  target.unit = ginv(as.matrix(best_profile))* target_cosine * norm(best_profile, type="2")
  sum.target.unit = sum(target.unit)
  target.1 = target.unit/sum.target.unit
  target_profile = target.1 * year_consump
  target_profile = data.frame(t(target_profile))
  colnames(target_profile) = colnames(profile_percent)
  for(i in 1:8)
  {
    if(target_profile[1,i] <= 1)
      target_profile[1,i] = 0
  }
  return(target_profile)
}


AZ.target.2025 = get.target(target.AZ.2025, opt.percent["Arizona",], AZ.imp$sum[51])



cosine.similarity(AZ.target.2025, opt.percent["Arizona",]*AZ.imp$sum[51])





TX.opt.BTU = data.frame(matrix(NA, nrow = 50, ncol = 8))
colnames(TX.opt.BTU) = colnames(TX.imp)[1:8]
rownames(TX.opt.BTU) = 1960:2009
for(i in 1:nrow(TX.opt.BTU))
{
  TX.opt.BTU[i,] = TX.imp$sum[i]*opt.percent["Texas",]
}
TX.opt.BTU$sim.score = 0
for(i in 1:nrow(TX.opt.BTU))
{
  TX.opt.BTU$sim.score[i] = cosine.similarity(TX.opt.BTU[i,1:8], TX.imp[i,1:8], .do.norm = TRUE)
}




##TEST REVERSE PCR REGRESSION
AZ.2009.pred = as.matrix(AZ.pc["50","PC1"]) %*% ginv(as.matrix(AZ.pr$rotation)[,1])
colnames(AZ.2009.pred) = colnames(AZ.wide)[2:584]
cosine.similarity(AZ.2009.pred[1,], AZ.wide[50,2:584], .do.norm = TRUE)
AZ.2009.pred[1,"CLTCB"]
AZ.wide[50,"CLTCB"]
abs(AZ.2009.pred[1,"CLTCB"] - AZ.wide[50,"CLTCB"])/AZ.wide[50,"CLTCB"]


lm.confirm = lm(CLTCB ~ Year, data = AZ.wide[1:50,])
newdata.2025$CLTCB = predict(lm.confirm, newdata = newdata.2025)
newdata.2050$CLTCB = predict(lm.confirm, newdata = newdata.2050)

newdata.2025$CLTCB
AZ.wide[51,"CLTCB"]
abs(newdata.2025$CLTCB - AZ.wide[51,"CLTCB"])/newdata.2025$CLTCB

newdata.2050$CLTCB
AZ.wide[52,"CLTCB"]
abs(newdata.2050$CLTCB - AZ.wide[52,"CLTCB"])/newdata.2050$CLTCB


lm.confirm = lm(SOTCB ~ Year, data = AZ.wide[1:50,])
newdata.2025$SOTCB = predict(lm.confirm, newdata = newdata.2025)
newdata.2050$SOTCB = predict(lm.confirm, newdata = newdata.2050)

newdata.2025$SOTCB
AZ.wide[51,"SOTCB"]
abs(newdata.2025$SOTCB - AZ.wide[51,"SOTCB"])/newdata.2025$SOTCB

newdata.2050$SOTCB
AZ.wide[52,"SOTCB"]
abs(newdata.2050$SOTCB - AZ.wide[52,"SOTCB"])/newdata.2050$SOTCB








#IB try again - interpretable
pairs(~percent.clean + Year + TPOPP + GDPRX + TETGR + ESTCD, data = AZ.IB)
AZ.IB.intp = lm(percent.clean ~ Year + TPOPP + GDPRX + TETGR + ESTCD , data = AZ.IB)
summary(AZ.IB.intp)
#Multiple R-squared:  0.5395,	Adjusted R-squared:  0.4871 
#F-statistic: 10.31 on 5 and 44 DF,  p-value: 1.38e-06
library(glmnet)
grid = 10^seq(10,-2,length=100)
AZ.X = as.matrix(AZ.wide[1:50,c("Year", "TPOPP", "GDPRX", "TETGR", "ESTCD")])
AZ.Y = as.matrix(AZ.wide[1:50,"percent.clean"])
AZ.IB.ridge = cv.glmnet(AZ.X, AZ.Y, alpha=0, lambda=grid, nfolds = 3)
plot(AZ.IB.ridge)
plot(glmnet(AZ.X, AZ.Y, alpha=0, lambda=grid))
abline(v = AZ.IB.ridge$lambda.min)
AZ.out = glmnet(AZ.X, AZ.Y, alpha=0, lambda = AZ.IB.ridge$lambda.min)
AZ.out$dev.ratio
predict(AZ.out, type="coef", s = AZ.IB.ridge$lambda.min)
plot(predict(AZ.out, type="response", newx = AZ.X, s = AZ.IB.ridge$lambda.min) ~ AZ.Y, 
     xlim = c(0,0.3), ylim = c(0,0.3), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Arizona", col = "springgreen4", pch = 19)
abline(0,1)


pairs(~percent.clean + Year + TPOPP + GDPRX + TETGR + ESTCD, data = CA.IB)
CA.IB.intp = lm(percent.clean ~ Year + TPOPP + GDPRX + TETGR + ESTCD, data = CA.IB)
CA.IB.intp = stepAIC(CA.IB.intp, trace = 0)
summary(CA.IB.intp)
#Multiple R-squared:  0.5785,	Adjusted R-squared:  0.4704 
#F-statistic: 5.353 on 10 and 39 DF,  p-value: 6.076e-05
CA.X = as.matrix(CA.wide[1:50,c("Year", "TPOPP", "GDPRX", "TETGR", "ESTCD")])
CA.Y = as.matrix(CA.wide[1:50,"percent.clean"])
CA.IB.ridge = cv.glmnet(CA.X, CA.Y, alpha=0, lambda=grid, nfolds = 3)
plot(CA.IB.ridge)
plot(glmnet(CA.X, CA.Y, alpha=0, lambda=grid))
abline(v = CA.IB.ridge$lambda.min)
CA.out = glmnet(CA.X, CA.Y, alpha=0, lambda = CA.IB.ridge$lambda.min)
preplot(predict(CA.out, type="response", newx = CA.X, s = CA.IB.ridge$lambda.min) ~ CA.Y, 
     xlim = c(0,0.3), ylim = c(0,0.3), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "California", col = "springgreen4", pch = 19)
abline(0,1)


pairs(~percent.clean + Year + TPOPP + GDPRX + TETGR + ESTCD, data = NM.IB)
NM.IB.intp = lm(percent.clean ~ Year + TPOPP + GDPRX + TETGR + ESTCD, data = NM.IB)
NM.IB.intp = stepAIC(NM.IB.intp, trace = 0)
summary(NM.IB.intp)
#Multiple R-squared:  0.6747,	Adjusted R-squared:  0.6377 
#F-statistic: 18.25 on 5 and 44 DF,  p-value: 9.002e-10
NM.X = as.matrix(NM.wide[1:50,c("Year", "TPOPP", "GDPRX", "TETGR", "ESTCD")])
NM.Y = as.matrix(NM.wide[1:50,"percent.clean"])
NM.IB.ridge = cv.glmnet(NM.X, NM.Y, alpha=0, lambda=grid, nfolds = 3)
plot(NM.IB.ridge)
plot(glmnet(NM.X, NM.Y, alpha=0, lambda=grid))
abline(v = NM.IB.ridge$lambda.min)
NM.out = glmnet(NM.X, NM.Y, alpha=0, lambda = NM.IB.ridge$lambda.min)
predict(NM.out, type="coef", s = NM.IB.ridge$lambda.min)
plot(predict(NM.out, type="response", newx = NM.X, s = NM.IB.ridge$lambda.min) ~ NM.Y, 
     xlim = c(0,0.3), ylim = c(0,0.3), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "New Mexico", col = "springgreen4", pch = 19)
abline(0,1)


pairs(~percent.clean + Year + TPOPP + GDPRX + TETGR + ESTCD, data = TX.IB)
TX.IB.intp = lm(percent.clean ~ Year + TPOPP + GDPRX + TETGR + ESTCD, data = TX.IB)
TX.IB.intp = stepAIC(TX.IB.intp, trace = 0)
summary(TX.IB.intp)
#Multiple R-squared:  0.6366,	Adjusted R-squared:  0.5953 
#F-statistic: 15.42 on 5 and 44 DF,  p-value: 9.454e-09
TX.X = as.matrix(TX.wide[1:50,c("Year", "TPOPP", "GDPRX", "TETGR", "ESTCD")])
TX.Y = as.matrix(TX.wide[1:50,"percent.clean"])
TX.IB.ridge = cv.glmnet(TX.X, TX.Y, alpha=0, lambda=grid, nfolds = 3)
plot(TX.IB.ridge)
plot(glmnet(TX.X, TX.Y, alpha=0, lambda=grid))
abline(v = TX.IB.ridge$lambda.min)
TX.out = glmnet(TX.X, TX.Y, alpha=0, lambda = TX.IB.ridge$lambda.min)
predict(TX.out, type="coef", s = TX.IB.ridge$lambda.min)
plot(predict(TX.out, type="response", newx = TX.X, s = TX.IB.ridge$lambda.min) ~ TX.Y, 
     xlim = c(0,0.3), ylim = c(0,0.3), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Texas", col = "springgreen4", pch = 19)
abline(0,1)


#plot glmnet
par(mfrow= c(2,2))
plot(glmnet(AZ.X, AZ.Y, alpha=0, lambda=grid), main = "Arizona")
plot(glmnet(CA.X, CA.Y, alpha=0, lambda=grid), main = "California")
plot(glmnet(NM.X, NM.Y, alpha=0, lambda=grid), main = "New Mexico")
plot(glmnet(TX.X, TX.Y, alpha=0, lambda=grid), main = "Texas")

#plot glmnet fits
plot(predict(AZ.out, type="response", newx = AZ.X, s = AZ.IB.ridge$lambda.min) ~ AZ.Y, 
     xlim = c(0,0.2), ylim = c(0,0.2), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Arizona", col = "springgreen4", pch = 19)
abline(0,1)
plot(predict(CA.out, type="response", newx = CA.X, s = CA.IB.ridge$lambda.min) ~ CA.Y, 
     xlim = c(0,0.2), ylim = c(0,0.2), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "California", col = "springgreen4", pch = 19)
abline(0,1)
plot(predict(NM.out, type="response", newx = NM.X, s = NM.IB.ridge$lambda.min) ~ NM.Y, 
     xlim = c(0,0.2), ylim = c(0,0.2), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "New Mexico", col = "springgreen4", pch = 19)
abline(0,1)
plot(predict(TX.out, type="response", newx = TX.X, s = TX.IB.ridge$lambda.min) ~ TX.Y, 
     xlim = c(0,0.2), ylim = c(0,0.2), xlab = "Percent clean energy", 
     ylab = "Predicted percent clean energy", main = "Texas", col = "springgreen4", pch = 19)
abline(0,1)

AZ.out$dev.ratio
CA.out$dev.ratio
NM.out$dev.ratio
TX.out$dev.ratio









############GOALS
maxdiff = function(df.imp)
{
  maxdiffs = numeric(8)
  names(maxdiffs) = colnames(df.imp)[1:8]
  for(i in 1:8)
  {
    diffs = numeric(49)
    for(j in 1:49)
    {
      diffs[j] = df.imp[j+1,i] - df.imp[j,i]
    }
    maxdiffs[i] = max(diffs)
  }
  return(maxdiffs)
}

AZ.maxdiffs = maxdiff(AZ.imp)
CA.maxdiffs = maxdiff(CA.imp)
NM.maxdiffs = maxdiff(NM.imp)
TX.maxdiffs = maxdiff(TX.imp)

AZ.2025.goals = data.frame(16*AZ.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])
CA.2025.goals = data.frame(16*CA.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])
NM.2025.goals = data.frame(16*NM.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])
TX.2025.goals = data.frame(16*TX.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])


AZ.2050.goals = data.frame((25*1.1)*AZ.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")] + 16*AZ.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])
CA.2050.goals = data.frame((25*1.1)*CA.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")] + 16*CA.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])
NM.2050.goals = data.frame((25*1.1)*NM.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")] + 16*NM.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])
TX.2050.goals = data.frame((25*1.1)*TX.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")] + 16*TX.maxdiffs[c("EMTCB", "GETCB", "SOTCB", "WYTCB")])

AZ.goals = AZ.imp[50:52,]
CA.goals = CA.imp[50:52,]
NM.goals = NM.imp[50:52,]
TX.goals = TX.imp[50:52,]

empty = data.frame(matrix(NA, nrow = 2, ncol = 9))
colnames(empty) = colnames(AZ.goals)

AZ.goals = rbind(AZ.goals, empty)
CA.goals = rbind(CA.goals, empty)
NM.goals = rbind(NM.goals, empty)
TX.goals = rbind(TX.goals, empty)

AZ.goals[4,c(2,3,7,8)] = AZ.goals[1,c(2,3,7,8)] + AZ.2025.goals[1,]
AZ.goals[5,c(2,3,7,8)] = AZ.goals[1,c(2,3,7,8)] + AZ.2050.goals[1,]
AZ.goals[4,4] = max(AZ.imp$HYTCB[1:50])
AZ.goals[5,4] = max(AZ.imp$HYTCB[1:50])
AZ.goals[4,9] = AZ.goals[2,9]
AZ.goals[5,9] = AZ.goals[3,9]
AZ.goals$sum.clean = AZ.goals$EMTCB + AZ.goals$GETCB + AZ.goals$HYTCB + AZ.goals$SOTCB + AZ.goals$WYTCB
AZ.goals[4,1] = max(c(AZ.imp[51,1] - AZ.goals$sum.clean[4]/3, 0))
AZ.goals[4,5] = max(c(AZ.imp[51,5] - AZ.goals$sum.clean[4]/3, 0))
AZ.goals[4,6] = max(c(AZ.imp[51,6] - AZ.goals$sum.clean[4]/3, 0))
AZ.goals[5,1] = max(c(AZ.imp[52,1] - AZ.goals$sum.clean[5]/3, 0))
AZ.goals[5,5] = max(c(AZ.imp[52,5] - AZ.goals$sum.clean[5]/3, 0))
AZ.goals[5,6] = max(c(AZ.imp[52,6] - AZ.goals$sum.clean[5]/3, 0))


CA.goals[4,c(2,3,7,8)] = CA.goals[1,c(2,3,7,8)] + CA.2025.goals[1,]
CA.goals[5,c(2,3,7,8)] = CA.goals[1,c(2,3,7,8)] + CA.2050.goals[1,]
CA.goals[4,4] = max(CA.imp$HYTCB[1:50])
CA.goals[5,4] = max(CA.imp$HYTCB[1:50])
CA.goals[4,9] = CA.goals[2,9]
CA.goals[5,9] = CA.goals[3,9]
CA.goals$sum.clean = CA.goals$EMTCB + CA.goals$GETCB + CA.goals$HYTCB + CA.goals$SOTCB + CA.goals$WYTCB
CA.goals[4,1] = max(c(CA.imp[51,1] - CA.goals$sum.clean[4]/3, 0))
CA.goals[4,5] = max(c(CA.imp[51,5] - CA.goals$sum.clean[4]/3, 0))
CA.goals[4,6] = max(c(CA.imp[51,6] - CA.goals$sum.clean[4]/3, 0))
CA.goals[5,1] = max(c(CA.imp[52,1] - CA.goals$sum.clean[5]/3, 0))
CA.goals[5,5] = max(c(CA.imp[52,5] - CA.goals$sum.clean[5]/3, 0))
CA.goals[5,6] = max(c(CA.imp[52,6] - CA.goals$sum.clean[5]/3, 0))



NM.goals[4,c(2,3,7,8)] = NM.goals[1,c(2,3,7,8)] + NM.2025.goals[1,]
NM.goals[5,c(2,3,7,8)] = NM.goals[1,c(2,3,7,8)] + NM.2050.goals[1,]
NM.goals[4,4] = max(NM.imp$HYTCB[1:50])
NM.goals[5,4] = max(NM.imp$HYTCB[1:50])
NM.goals[4,9] = NM.goals[2,9]
NM.goals[5,9] = NM.goals[3,9]
NM.goals$sum.clean = NM.goals$EMTCB + NM.goals$GETCB + NM.goals$HYTCB + NM.goals$SOTCB + NM.goals$WYTCB
NM.goals[4,1] = max(c(NM.imp[51,1] - NM.goals$sum.clean[4]/3, 0))
NM.goals[4,5] = max(c(NM.imp[51,5] - NM.goals$sum.clean[4]/3, 0))
NM.goals[4,6] = max(c(NM.imp[51,6] - NM.goals$sum.clean[4]/3, 0))
NM.goals[5,1] = max(c(NM.imp[52,1] - NM.goals$sum.clean[5]/3, 0))
NM.goals[5,5] = max(c(NM.imp[52,5] - NM.goals$sum.clean[5]/3, 0))
NM.goals[5,6] = max(c(NM.imp[52,6] - NM.goals$sum.clean[5]/3, 0))




TX.goals[4,c(2,3,7,8)] = TX.goals[1,c(2,3,7,8)] + TX.2025.goals[1,]
TX.goals[5,c(2,3,7,8)] = TX.goals[1,c(2,3,7,8)] + TX.2050.goals[1,]
TX.goals[4,4] = max(TX.imp$HYTCB[1:50])
TX.goals[5,4] = max(TX.imp$HYTCB[1:50])
TX.goals[4,9] = TX.goals[2,9]
TX.goals[5,9] = TX.goals[3,9]
TX.goals$sum.clean = TX.goals$EMTCB + TX.goals$GETCB + TX.goals$HYTCB + TX.goals$SOTCB + TX.goals$WYTCB
TX.goals[4,1] = max(c(TX.imp[51,1] - TX.goals$sum.clean[4]/3, 0))
TX.goals[4,5] = max(c(TX.imp[51,5] - TX.goals$sum.clean[4]/3, 0))
TX.goals[4,6] = max(c(TX.imp[51,6] - TX.goals$sum.clean[4]/3, 0))
TX.goals[5,1] = max(c(TX.imp[52,1] - TX.goals$sum.clean[5]/3, 0))
TX.goals[5,5] = max(c(TX.imp[52,5] - TX.goals$sum.clean[5]/3, 0))
TX.goals[5,6] = max(c(TX.imp[52,6] - TX.goals$sum.clean[5]/3, 0))




par(mfrow = c(1,3))
barplot(as.numeric(AZ.goals[1,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2009", ylim = c(0, 600000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(AZ.goals[4,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2025", ylim = c(0, 600000))
barplot(as.numeric(AZ.goals[5,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2050", ylim = c(0, 600000))




barplot(as.numeric(CA.goals[1,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2009", ylim = c(0, 4000000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(CA.goals[4,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2025", ylim = c(0, 4000000))
barplot(as.numeric(CA.goals[5,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2050", ylim = c(0, 4000000))


barplot(as.numeric(NM.goals[1,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2009", ylim = c(0, 360000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(NM.goals[4,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2025", ylim = c(0, 360000))
barplot(as.numeric(NM.goals[5,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2050", ylim = c(0, 360000))



barplot(as.numeric(TX.goals[1,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2009", ylim = c(0, 5300000), ylab = "Total Energy Consumption in BTU")
barplot(as.numeric(TX.goals[4,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2025", ylim = c(0, 5300000))
barplot(as.numeric(TX.goals[5,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB")]), 
        col = c("tomato", "orangered", "red3",
                "olivedrab4", "olivedrab3", "olivedrab2", "springgreen2", "springgreen3"),
        main = "2050", ylim = c(0, 5300000))


goal.sim = function(goals.vect, opt.percent.vect)
{
  opt.vect = opt.percent.vect*goals.vect[,"sum"]
  sim = cosine.similarity(opt.vect, goals.vect[,1:8])
  return(sim)
}





goals.table = rbind(AZ.goals[4:5,1:9], CA.goals[4:5,1:9])
goals.table = rbind(goals.table, NM.goals[4:5,1:9])
goals.table = rbind(goals.table, TX.goals[4:5,1:9])
rownames(goals.table) = c("Arizona 2025", "Arizona 2050",
                          "California 2025", "California 2050",
                          "New Mexico 2025", "New Mexico 2050",
                          "Texas 2025", "Texas 2050")
goals.table$sim = 0
for(i in 1:nrow(goals.table))
{
  goals.table$sim[i] = goal.sim(goals.table[i,], opt.percent["Default",])
}
write.table(goals.table[,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB","sim")], file = "/Users/katiebalcewicz/desktop/COMAP 2018/goals_table.csv", row.names = TRUE, col.names = TRUE, sep = ",")


rownames(AZ.goals) = c("2009", "2025 Predicted", "2050 Predicted", "2025 Goal", "2050 Goal")
rownames(CA.goals) = c("2009", "2025 Predicted", "2050 Predicted", "2025 Goal", "2050 Goal")
rownames(NM.goals) = c("2009", "2025 Predicted", "2050 Predicted", "2025 Goal", "2050 Goal")
rownames(TX.goals) = c("2009", "2025 Predicted", "2050 Predicted", "2025 Goal", "2050 Goal")

AZ.goals$similarity = 0
for(i in 1:nrow(AZ.goals))
{
  AZ.goals$similarity[i] = goal.sim(AZ.goals[i,], opt.percent["Default",])
}

CA.goals$similarity = 0
for(i in 1:nrow(CA.goals))
{
  CA.goals$similarity[i] = goal.sim(CA.goals[i,], opt.percent["Default",])
}

NM.goals$similarity = 0
for(i in 1:nrow(NM.goals))
{
  NM.goals$similarity[i] = goal.sim(NM.goals[i,], opt.percent["Default",])
}

TX.goals$similarity = 0
for(i in 1:nrow(TX.goals))
{
  TX.goals$similarity[i] = goal.sim(TX.goals[i,], opt.percent["Default",])
}










AZ.goals = AZ.goals[,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "similarity")]
CA.goals = CA.goals[,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "similarity")]
NM.goals = NM.goals[,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "similarity")]
TX.goals = TX.goals[,c("CLTCB", "NNTCB", "PMTCB", "EMTCB", "GETCB", "HYTCB", "SOTCB", "WYTCB", "similarity")]

AZ.goals.out = round(AZ.goals[,1:8]/rowSums(AZ.goals[,1:8]), 4)
AZ.goals.out$similarity = AZ.goals$similarity
colnames(AZ.goals.out) = c("Coal", "Natural Gas", "Petroleum Products", "Fuel Ethanol ", "Geothermal", "Hydroelectric", "Solar", "Wind", "Profile Score")
write.table(AZ.goals.out, file = "/Users/katiebalcewicz/desktop/COMAP 2018/arizona.csv", row.names = TRUE, col.names = TRUE, sep = ",")



CA.goals.out = round(CA.goals[,1:8]/rowSums(CA.goals[,1:8]), 4)
CA.goals.out$similarity = CA.goals$similarity
colnames(CA.goals.out) = c("Coal", "Natural Gas", "Petroleum Products", "Fuel Ethanol ", "Geothermal", "Hydroelectric", "Solar", "Wind", "Profile Score")
write.table(CA.goals.out, file = "/Users/katiebalcewicz/desktop/COMAP 2018/california.csv", row.names = TRUE, col.names = TRUE, sep = ",")



NM.goals.out = round(NM.goals[,1:8]/rowSums(NM.goals[,1:8]), 4)
NM.goals.out$similarity = NM.goals$similarity
colnames(NM.goals.out) = c("Coal", "Natural Gas", "Petroleum Products", "Fuel Ethanol ", "Geothermal", "Hydroelectric", "Solar", "Wind", "Profile Score")
write.table(NM.goals.out, file = "/Users/katiebalcewicz/desktop/COMAP 2018/new mexico.csv", row.names = TRUE, col.names = TRUE, sep = ",")


TX.goals.out = round(TX.goals[,1:8]/rowSums(TX.goals[,1:8]), 4)
TX.goals.out$similarity = TX.goals$similarity
colnames(TX.goals.out) = c("Coal", "Natural Gas", "Petroleum Products", "Fuel Ethanol ", "Geothermal", "Hydroelectric", "Solar", "Wind", "Profile Score")
write.table(TX.goals.out, file = "/Users/katiebalcewicz/desktop/COMAP 2018/texas.csv", row.names = TRUE, col.names = TRUE, sep = ",")


