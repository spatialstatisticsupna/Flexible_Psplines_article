rm(list=ls())
library(INLA)
library(RColorBrewer)
library(sf)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#################################################
## Load the data and its asociated cartography ##
#################################################
load("../data/Carto_ESP.Rdata")

head(Carto_ESP)
plot(Carto_ESP$geometry, axes=TRUE)

Data <- read.table("../data/BreastCancer_ESP.txt", header=TRUE)
Data <- Data[order(Data$Prov,Data$Year,Data$Age.group), ] 
str(Data)

T <- length(unique(Data$Year))
S <- length(unique(Data$Province))
A <- length(unique(Data$Age.group))

t.from <- min(Data$Year)
t.to <- max(Data$Year)


############################################################################################
## Load the final model fitted with INLA: Model9-SI (with posterior patterns),            ##
## which can be download from: https://emi-sstcdapp.unavarra.es/Flexible_Psplines_article ##
############################################################################################
load(url("https://emi-sstcdapp.unavarra.es/Flexible_Psplines_article/Model9_SI.Rdata"))

Model <- M9.SI
summary(Model)


#############################################################################################################################
## Figure 3: Posterior mean of female breast cancer spatial effects common to all years and age groups (left),             ##
##           posterior mean of temporal effects (with 95% credible bands) common to all provinces and age groups (middle), ##
##           and posterior mean of the age pattern (with 95% credible bands) common to all years and provinces (right).    ##
#############################################################################################################################


## Spatial pattern ##
Carto_ESP$zeta <- c(unlist(lapply(Model$marginals.lincomb.derived[2:(S+1)], function(x) inla.emarginal(exp,x))), NA)

paleta <- brewer.pal(6,"YlOrRd")
values <- c(0.7,0.8,0.9,1,1.1,1.2,1.3)

Map.fig1 <- tm_shape(Carto_ESP) +
  tm_polygons(col="zeta", title="", palette=paleta,
              colorNA=NULL, showNA=FALSE,
              legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title="Spatial pattern", main.title.position=0.2,
            legend.outside=T, legend.outside.position="right", legend.frame=F)
print(Map.fig1)


## Temporal pattern ##
temporal <- unlist(lapply(Model$marginals.lincomb.derived[(S+2):(S+T+1)], function(x) inla.emarginal(exp,x)))
aux <- lapply(Model$marginals.lincomb.derived[(S+2):(S+T+1)], function(x) inla.tmarginal(exp,x))
q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025,x)))
q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975,x)))

x <- 1:T
plot(range(x),c(0.7,1.3),type="n",xlab="",ylab="", xaxt="n", cex.axis=1.5)
title(main="Temporal pattern", cex.main=2.5, line=2)
axis(1, at=seq(1,T,4), labels=seq(t.from,t.to,4), las=0, cex.axis=1.5)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(temporal, lwd=2)
abline(h=1, lty=2)


## Age-group pattern ##
edad <- unlist(lapply(Model$marginals.lincomb.derived[(S+T+2):(S+T+A+1)], function(x) inla.emarginal(exp,x)))
aux <- lapply(Model$marginals.lincomb.derived[(S+T+2):(S+T+A+1)], function(x) inla.tmarginal(exp,x))
q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025,x)))
q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975,x)))

x <- 1:A
plot(range(x),c(0,4),type="n",xlab="",ylab="", xaxt="n", cex.axis=1.5)
title(main="Age pattern", cex.main=2.5, line=2)
axis(1, at=seq(1,A), labels=unique(Data$Age.label), las=2, cex.axis=1.5)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
polygon(X.Vec, Y.Vec, col = "red4", border = NA)
#lines(edad, lwd=1)
abline(h=1,lty=2)


###################################################################################################
## Figure 4: Posterior mean of the area-age interaction effect in female breast cancer mortality ##
###################################################################################################
delta <- unlist(lapply(Model$marginals.lincomb.derived[(S+T+A+S*T+2):(S+T+A+S*T+S*A+1)], function(x) inla.emarginal(exp,x)))
risk.matrix <- rbind(matrix(delta,S,A,byrow=T), rep(NA,A))

risk.data.frame <- data.frame(Carto_ESP$NAME,risk.matrix)
colnames(risk.data.frame) <- c("NAME",paste("Age",seq(1,A),sep=""))

carto <- merge(Carto_ESP,risk.data.frame)

paleta <- brewer.pal(8,"YlOrRd")
values <- c(0.79,0.84,0.89,0.94,1.00,1.06,1.12,1.19,1.26)

Map.fig4 <- tm_shape(carto) +
  tm_polygons(col=paste("Age",seq(1,A),sep=""), title="",
              colorNA=NULL, showNA=FALSE,
              palette=paleta, legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title="Area-age interaction", main.title.position=0.2, panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            panel.labels=unique(Data$Age.label),
            legend.outside.size=0.3, outer.margins=c(0.0,0.03,0.01,0)) + 
  tm_facets(nrow=4, ncol=3)

print(Map.fig4)


########################################################################################################
## Figure 5: Posterior mean estimate of time-age interaction effect in female breast cancer mortality ##
########################################################################################################
delta <- unlist(lapply(Model$marginals.lincomb.derived[(S+T+A+S*T+S*A+2):(S+T+A+S*T+S*A+T*A+1)], function(x) inla.emarginal(exp,x)))
risk.matrix <- matrix(delta,A,T)

plot(c(1,T),c(0.7,1.4),type="n", xlab="",ylab="", xaxt="n", main="Age-year interaction")
axis(1, at=seq(1,T,4), labels=seq(t.from,t.to,4), las=0)

for(k in 1:A){
  lines(risk.matrix[k,], lwd=2, col=k, lty=c(rep(1,8),rep(2,3))[k])
}

legend(x=4, y=1.4, legend=unique(Data$Age.label)[1:5], col=1:5, lwd=2, lty=1, bty="n")
legend(x=12, y=1.4, legend=unique(Data$Age.label)[6:11], col=6:11, lwd=2, lty=rep(1:2,each=3), bty="n")


#########################################################################################
## Figure 6: Temporal evolution of mortality rates and 95% credible bands by age group ##
##           and province for three selected provinces: Madrid, Navarra, and Huelva.   ##
#########################################################################################
Data$fitted.values <- Model$summary.fitted.values$mean[1:(S*T*A)]
Data$fitted.values.q1 <- Model$summary.fitted.values$'0.025quant'[1:(S*T*A)]
Data$fitted.values.q2 <- Model$summary.fitted.values$'0.975quant'[1:(S*T*A)]
Data <- Data[order(Data$Year,Data$Province,Data$Age.group), ]

risk.M9 <- array(Data$fitted.values*100000,dim=c(A,S,T))
q1.M9 <- array(Data$fitted.values.q1*100000,dim=c(A,S,T))
q2.M9 <- array(Data$fitted.values.q2*100000,dim=c(A,S,T))

paleta <- c("#a6cee3","#b2df8a","#fb9a99","#fdbf6f","#ff7f00","#1f78b4","#33a02c","#e31a1c","#6a3d9a","#a6761d","#999999")

par(mfrow=c(3,2), pty="m")
for (j in c(28,31,21)) {
  x <- 1:T
  plot(c(1,T), c(0,70), type="n", main=Carto_ESP$NAME[j], xlab="", ylab="", xaxt="n", cex.main=2, cex.axis=1.5)
  axis(1, at=seq(1,T,5), labels=as.character(seq(t.from,t.to,5)), las=0, cex.axis=1.5)
  # legend("topright", legend=unique(Data$Age.label)[1:5], col=paleta[1:5], lwd=2, cex=1.2)
  
  for (i in 1:5){
    q1 <- q1.M9[i,j,]
    q2 <- q2.M9[i,j,]
    risk <- risk.M9[i,j,]
    
    X.Vec <- c(x, tail(x, 1), rev(x), x[1])
    Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
    
    color=rgb(t(col2rgb(paleta[i])), maxColorValue=255, alpha=150)
    polygon(X.Vec, Y.Vec, col=color, border = NA)
    
    lines(risk, col=paleta[i], lwd=2)
  }
  
  x <- 1:T
  plot(c(1,T), c(30,220), type="n", main=Carto_ESP$NAME[j], xlab="", ylab="", xaxt="n", cex.main=2, cex.axis=1.5)
  axis(1, at=seq(1,T,5), labels=as.character(seq(t.from,t.to,5)), las=0, cex.axis=1.5)
  # legend("topleft", legend=unique(Data$Age.label)[6:11], col=paleta[6:11], lwd=2, cex=1.2)
  
  for (i in 6:11){
    q1 <- q1.M9[i,j,]
    q2 <- q2.M9[i,j,]
    risk <- risk.M9[i,j,]
    
    X.Vec <- c(x, tail(x, 1), rev(x), x[1])
    Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
    
    color=rgb(t(col2rgb(paleta[i])), maxColorValue=255, alpha=150)
    polygon(X.Vec, Y.Vec, col=color, border = NA)
    
    lines(risk, col=paleta[i], lwd=2)
  }
}


#############################################################################################################
## Figure 7: Age-specific breast cancer mortality patterns in Spain (1985-2010) for the age groups [45-50) ##
#############################################################################################################
i <- 3 ## Age-group [45-50)

risk.matrix <- rbind(risk.M9[i,,],rep(NA,T))
risk.data.frame <- data.frame(Carto_ESP$NAME,risk.matrix)
colnames(risk.data.frame)<- c("NAME",paste("Year",seq(1,T),sep=""))

carto <- merge(Carto_ESP,risk.data.frame)


Obs.est <- (Data$fitted.values/Data$Pop)[Data$Age.group==i]
rate.age2 <- apply(matrix(Data$Obs[Data$Age.group==i],nrow=S,ncol=T),2,mean)/apply(matrix(Data$Pop[Data$Age.group==i],nrow=S,ncol=T),2,mean)*100000

paleta <- brewer.pal(8,"RdBu")[8:1]
values <- c(round(c(mean(rate.age2)/2,mean(rate.age2)/1.75,mean(rate.age2)/1.5,mean(rate.age2)/1.25, mean(rate.age2),mean(rate.age2)*1.25,mean(rate.age2)*1.5,mean(rate.age2)*1.75,mean(rate.age2)*2)))

Map.fig7 <- tm_shape(carto) +
  tm_polygons(col=paste("Year",seq(1,T),sep=""), title="",
              colorNA=NULL, showNA=FALSE,
              palette=paleta, legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") + 
  tm_layout(main.title=paste("Age-group",unique(Data$Age.label)[i]),
            main.title.position="center", panel.label.size=1.5,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            panel.labels=as.character(seq(t.from,t.to)),
            legend.outside.size=0.3, outer.margins=c(0.0,0.03,0.01,0)) + 
  tm_facets(nrow=4, ncol=7)

print(Map.fig7)
