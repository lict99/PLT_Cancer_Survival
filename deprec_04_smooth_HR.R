
# env settings ------------------------------------------------------------

library(smoothHR)
library(survival)
library(Hmisc)
library(magrittr)
load("01/tidy_data_diagnosis_after_attending.RData")
dir.create("04", FALSE)

# calculating -------------------------------------------------------------
tidy_data_dia_after_att <- transform(
  tidy_data_dia_after_att,
  futime = Date_end - Date_of_diagnosis
)

fit <- coxph(
  Surv(futime, OS) ~ rcspline.eval(Platelet, nk = 3, inclx = TRUE),
  data = tidy_data_dia_after_att,
  x = TRUE
)

hr <- smoothHR(data = tidy_data_dia_after_att, coxfit = fit)

p <- summary(hr$coxfit)

refvalue <- 400
smoothlogHR.point <- predict(
  hr,
  predictor = "Platelet",
  pred.value = refvalue,
  prediction.values = seq(min(tidy_data_dia_after_att$Platelet, na.rm = T), max(tidy_data_dia_after_att$Platelet, na.rm = T), length.out = 10000),
  conf.level = 0.95
) %>% 
  as.data.frame()


# plotting ----------------------------------------------------------------

violet <- "#89439B"

pdf("04/pan_cancer_HR_with_RCS_nk3_ref400.pdf",width = 5,height = 5)
par(mar = c(5, 4, 4, 4) + 0.3)
par(xpd=NA)

ylim.bot <- min(exp(smoothlogHR.point$LnHR),exp(smoothlogHR.point$`lower .95`),exp(smoothlogHR.point$`upper .95`))
ylim.top <- max(exp(smoothlogHR.point$LnHR),exp(smoothlogHR.point$`lower .95`),exp(smoothlogHR.point$`upper .95`))

dens <- density(tidy_data_dia_after_att$Platelet, na.rm = TRUE)

plot(dens$x,dens$y, col=ggplot2::alpha(violet,0.5), type="l", xlab = "", ylab = "",xaxt="n",yaxt="n")
polygon(dens$x,dens$y,col = ggplot2::alpha(violet,0.5),border = ggplot2::alpha(violet,0.5))
axis(side=4, at = pretty(range(dens$y))[-length(pretty(range(dens$y)))])
mtext("Fraction of population (Density)", side=4, line=2)

par(new=TRUE)
plot(smoothlogHR.point[,1],exp(smoothlogHR.point$LnHR), 
     xlab = "platelet",ylab = paste0("HR"),
     type = "l",ylim = c(ylim.bot,ylim.top),
     col="red",lwd=2) 
lines(smoothlogHR.point[,1],exp(smoothlogHR.point$`lower .95`),lty=2,lwd=1.5)
lines(smoothlogHR.point[,1],exp(smoothlogHR.point$`upper .95`),lty=2,lwd=1.5)
lines(x=range(smoothlogHR.point[,1]),y=c(1,1),lty=3,col="grey40",lwd=1.3)
#points(refvalue,1,pch=16,cex=1.2)
#text(refvalue + 5, 1.5, paste0("refvalue = ",refvalue))

legend("topright",
       paste0("P-overall ",ifelse(round(p$logtest[3],3) < 0.001,"< 0.001",round(p$logtest[3],3)),
              "\nP-non-linear = ",ifelse(round(p$coefficients[2,5],3) < 0.001,"< 0.001",round(p$coefficients[2,5],3))),
       bty="n",cex=0.8)

#legend("bottomleft",lty = c(1,2),col = c("red","black"),
#       c("Estimation","95% CI"),
#       bty="n",cex=0.8)
dev.off()
