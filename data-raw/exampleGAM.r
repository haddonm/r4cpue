

library(r4cpue)
library(r4maps)
library(mgcv)
library(nlme)
library(gamm4)

data(sps)
pick <- which((sps$Depth < 300) & (sps$Depth > 50))
sps1 <- droplevels(sps[pick,])
sps1$DepCat <- NA
sps1$DepCat <- trunc(sps1$Depth/25) * 25
sps1$LnCE <- NA
pick <- which((sps1$catch_kg > 0) & (sps1$Effort > 0))
sps1$LnCE[pick] <- log(sps1$catch_kg[pick]/sps1$Effort[pick])


labelM <‐ c("Year","Vessel","DepCat","Zone","Month","DayNight","Month:Zone")
sps2 <‐ makecategorical(labelM[1:6],sps1) # convert variables to factors
mods <‐ makemodels(labelM)
out1 <‐ standLM(mods,sps2,"FlatheadLM")
cat("Optimum Model: ",as.character(mods[[out1$Optimum]]),"\n\n")
cat("Summary of Optimum Model \n")
summary(out1)
cat("\n\n Anova of Optimum Model \n\n")
anova(out1$optModel)
model1 <‐ out1$optModel
plotstand(out1,bars=TRUE)



labelM <‐ c("Year","Vessel","DepCat","Zone","Month","DayNight","Month:Zone")
sps2 <‐ makecategorical(labelM[1:6],sps1) # convert variables to factors


model7 <‐ gam(LnCE ~ s(Long,Lat) + Vessel + DepCat + Month + Year, data = sps2)

answer7 <‐ getfact(model7,"Year")

answer7
round(answer7,3)
plotprep(width=5,height=8)
plot(model7)
addpoints(sps2)
abline(mod8,col=4)


labelM <‐ c("Year","Vessel","DepCat","Zone","Month","DayNight")
sps3 <‐ makecategorical(labelM[1:6],sps1) # convert variables to factors
ypts <‐ c(-42,-45.5); xpts <- c(144.0, 146.5)
mod8 <‐ lm(ypts ~ xpts)
coef(mod8)
sps3$LR <‐ NA
sps3$LR <‐ coef(mod8)[2] * sps3$Long + coef(mod8)[1]
pick <‐ which(sps3$LR > sps3$Lat)
if (length(pick) > 0) sps3 <‐ sps3[‐pick,]
model8 <‐ gam(LnCE ~ s(Long,Lat) + Vessel + DepCat + Zone + Month + Year, data = sps3)

answer8 <‐ getfact(model8,"Year")


plotprep(width=5,height=8)
plot(model8)
addpoints(sps3)

plotprep(width=7,height=4)
ymax <- max(answer8[,"Scaled"],na.rm=TRUE)*1.05
plot(2003:2014,answer8[,"Scaled"],type="p",pch=16,col=2,ylim=c(0,ymax),
     panel.first=grid())
lines(2003:2014,answer8[,"Scaled"],lwd=1,col=2)
lines(2003:2014,answer7[,"Scaled"],lwd=1,col=3)



















