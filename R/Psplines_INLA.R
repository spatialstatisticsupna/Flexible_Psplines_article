#########################################################################
####	INLA code to fit the models described in Goicoa et al. (2019)  ####
####	https://doi.org/10.1177/0962280217726802                       ####
#########################################################################
rm(list=ls())

## Download and install the R-INLA package (stable or testing version)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/testing"), dep=TRUE)
library(INLA)
library(splines)


###################################################
##  1. Read the data to be analyzed  		         ##
##  NOTE: Data must be ordered according to the  ##
##        Kronecker product defined in Section 3 ##
###################################################
Data <- read.table("../data/BreastCancer_ESP.txt", header=TRUE)
Data <- Data[order(Data$Year,Data$Province,Data$Age.group), ] 

T <- length(unique(Data$Year))
S <- length(unique(Data$Province))
A <- length(unique(Data$Age.group))


########################################################
##  2. Construct the B-spline basis for time and age  ##
########################################################
p <- 3 		    ## Cubic B-splines
q.time <- 7		## Number of internal intervals for time
q.age <- 3		## Number of internal intervals for age


## Marginal basis for time-periods ##
xt <- unique(Data$Year)
xt <- (xt-min(xt))/(max(xt)-min(xt))	## Scaled into the [0,1] interval

dis.t <- (max(xt)-min(xt))/q.time
xl.t <- min(xt)-dis.t*0.05
xr.t <- max(xt)+dis.t*0.05
dx.t <- (xr.t-xl.t)/q.time
knots.t <- seq(xl.t-p*dx.t, xr.t+p*dx.t, by=dx.t)

Bt <- spline.des(knots.t,xt,p+1)$design
kt <- ncol(Bt)


## Marginal basis for age-groups ##
xa <- unique(Data$Age.group)
xa <- (xa-min(xa))/(max(xa)-min(xa))	## Scaled into the [0,1] interval

dis.a <- (max(xa)-min(xa))/q.age
xl.a <- min(xa)-dis.a*0.05
xr.a <- max(xa)+dis.a*0.05
dx.a <- (xr.a-xl.a)/q.age
knots.a <- seq(xl.a-p*dx.a, xr.a+p*dx.a, by=dx.a)

Ba <- spline.des(knots.a,xa,p+1)$design
ka <- ncol(Ba)


###################################################################
##  3. Define the structure matrices for the main random effects ##
###################################################################

## Load the spatial neighbourhood matrix ##
g <- inla.read.graph("../data/Esp_prov_nb.inla")

Qs <- matrix(0, g$n, g$n)
for (i in 1:g$n){
	Qs[i,i]=g$nnbs[[i]]
	Qs[i,g$nbs[[i]]]=-1
}

Q.Leroux <- diag(S)-Qs


## Structure matrix for time effect (random walk of first order) ##
Dt <- diff(diag(kt),differences=1)
Pt <- t(Dt)%*%Dt


## Structure matrix for age effect (random walk of first order) ##
Da <- diff(diag(ka),differences=1)
Pa <- t(Da)%*%Da


#########################################################################################
##  4. Define appropriate hyperprior distributions using the "expression()" function   ##
##  - Unif(0,Inf) for standard deviations                                              ##
##	- Unif(0,1) for the spatial smoothing parameter						                         ##
#########################################################################################
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

lunif = "expression:
    a = 1;
    b = 1;
    beta = exp(theta)/(1+exp(theta));
    logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
    log_jacobian = log(beta*(1-beta));
    return(logdens+log_jacobian)"


###########################################################################
##  5. Define the 'formula' object for the models described in Section 3 ##
###########################################################################
ones.T <- matrix(1,T,1)
ones.S <- matrix(1,S,1)
ones.A <- matrix(1,A,1)

ones.kt <- matrix(1,kt,1)
ones.ka <- matrix(1,ka,1)


## M1: log(r_it) = intercept + phi_i + f(xt) + f(xa)
#####################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xt <- kronecker(Bt,kronecker(ones.S,ones.A))
Xa <- kronecker(ones.T,kronecker(ones.S,Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt+ka)),
			ID.area=c(NA,1:S,rep(NA,kt+ka)),
			ID.year=c(rep(NA,1+S),1:kt,rep(NA,ka)),
			ID.age=c(rep(NA,1+S+kt),1:ka))

f.M1 <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
		f(ID.age, model="rw1", hyper=list(prec=list(prior=sdunif)))

M1 <- inla(f.M1, family="poisson", data=Data.inla, E=E,
	     control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xt,Xa), link=1, cdf=c(log(1))),
	     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
	     control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M2: log(r_it) = intercept + phi_i + f(xt) + fs(xa)
######################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xt <- kronecker(Bt,kronecker(ones.S,ones.A))
Xsa <- kronecker(ones.T,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt+S*ka)),
			ID.area=c(NA,1:S,rep(NA,kt+S*ka)),
			ID.year=c(rep(NA,1+S),1:kt,rep(NA,S*ka)),
			ID.area.age=c(rep(NA,1+S+kt),1:(S*ka)))

R <- kronecker(diag(S),Pa)
r.def <- S
A.constr <- kronecker(diag(S),t(ones.ka))

f.M2 <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
		f(ID.area.age, model="generic0", Cmatrix=R, rankdef=r.def,
		  constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  extraconstr=list(A=A.constr, e=rep(0,S)))

M2 <- inla(f.M2, family="poisson", data=Data.inla, E=E,
	     control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xt,Xsa), link=1, cdf=c(log(1))),
	     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
	     control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M3: log(r_it) = intercept + phi_i + fs(xt) + f(xa)
######################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xa <- kronecker(ones.T,kronecker(ones.S,Ba))
Xst <- kronecker(Bt,kronecker(diag(S),ones.A))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+ka+S*kt)),
			ID.area=c(NA,1:S,rep(NA,ka+S*kt)),
			ID.age=c(rep(NA,1+S),1:ka,rep(NA,S*kt)),
			ID.area.year=c(rep(NA,1+S+ka),1:(S*kt)))

R <- kronecker(Pt,diag(S))
r.def <- S
A.constr <- kronecker(t(ones.kt),diag(S))

f.M3 <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.age, model="rw1", hyper=list(prec=list(prior=sdunif))) +
		f(ID.area.year, model="generic0", Cmatrix=R, rankdef=r.def,
		  constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  extraconstr=list(A=A.constr, e=rep(0,S)))

M3 <- inla(f.M3, family="poisson", data=Data.inla, E=E,
	     control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xa,Xst), link=1, cdf=c(log(1))),
	     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
	     control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M4: log(r_it) = intercept + phi_i + fs(xt) + fs(xa)
#######################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xst <- kronecker(Bt,kronecker(diag(S),ones.A))
Xsa <- kronecker(ones.T,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+S*kt+S*ka)),
			ID.area=c(NA,1:S,rep(NA,S*kt+S*ka)),
			ID.area.year=c(rep(NA,1+S),1:(S*kt),rep(NA,S*ka)),
			ID.area.age=c(rep(NA,1+S+S*kt),1:(S*ka)))

R1 <- kronecker(Pt,diag(S))
r1.def <- S
A1.constr <- kronecker(t(ones.kt),diag(S))

R2 <- kronecker(diag(S),Pa)
r2.def <- S
A2.constr <- kronecker(diag(S),t(ones.ka))

f.M4 <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.area.year, model="generic0", Cmatrix=R1, rankdef=r1.def,
		  constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  extraconstr=list(A=A1.constr, e=rep(0,S))) +
		f(ID.area.age, model="generic0", Cmatrix=R2, rankdef=r2.def,
		  constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  extraconstr=list(A=A2.constr, e=rep(0,S)))

M4 <- inla(f.M4, family="poisson", data=Data.inla, E=E,
	     control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xst,Xsa), link=1, cdf=c(log(1))),
	     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
	     control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M5-FRS: log(r_it) = intercept + phi_i + f(xt,xa)
####################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xta <- kronecker(Bt,kronecker(ones.S,Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt*ka)),
			ID.area=c(NA,1:S,rep(NA,kt*ka)),
			ID.year.age=c(rep(NA,1+S),1:(kt*ka)))

R <- kronecker(Pt,diag(ka)) + kronecker(diag(kt),Pa)

f.M5.FRS <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.year.age, model="generic0", Cmatrix=R, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif)))

M5.FRS <- inla(f.M5.FRS, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xta), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M5-SI: log(r_it) = intercept + phi_i + f(xt,xa)
###################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xta <- kronecker(Bt,kronecker(ones.S,Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt*ka)),
			ID.area=c(NA,1:S,rep(NA,kt*ka)),
			ID.year.age=c(rep(NA,1+S),1:(kt*ka)))

R <- list(inla.as.sparse(kronecker(Pt,diag(ka))),
	    inla.as.sparse(kronecker(diag(kt),Pa)))

f.M5.SI <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.year.age, model="generic3", Cmatrix=R, constr=TRUE,
		  hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))

M5.SI <- inla(f.M5.SI, family="poisson", data=Data.inla, E=E,
		  control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xta), link=1, cdf=c(log(1))),
		  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		  control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M6-FRS: log(r_it) = intercept + phi_i + fs(xt,xa)
#####################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S),1:(kt*S*ka)))

R <- kronecker(Pt,kronecker(diag(S),diag(ka))) + kronecker(diag(kt),kronecker(diag(S),Pa))
r.def <- S
A.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M6.FRS <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.year.area.age, model="generic0", Cmatrix=R, rankdef=r.def,
		  constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  extraconstr=list(A=A.constr, e=rep(0,S)))

M6.FRS <- inla(f.M6.FRS, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M6-SI: log(r_it) = intercept + phi_i + fs(xt,xa)
####################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S),1:(kt*S*ka)))

R <- list(inla.as.sparse(kronecker(Pt,kronecker(diag(S),diag(ka)))),
	    inla.as.sparse(kronecker(diag(kt),kronecker(diag(S),Pa))))
r.def <- S
A.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M6.SI <- O ~ -1 + Intercept +
		f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		f(ID.year.area.age, model="generic3", Cmatrix=R, constr=TRUE,
		  hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)),
		  extraconstr=list(A=A.constr, e=rep(0,S)))

M6.SI <- inla(f.M6.SI, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M7-FRS: log(r_it) = intercept + phi_i + f(xt) + f(xa) + fs(xt,xa)
#####################################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xt <- kronecker(Bt,kronecker(ones.S,ones.A))
Xa <- kronecker(ones.T,kronecker(ones.S,Ba))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt+ka+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,kt+ka+kt*S*ka)),
			ID.year=c(rep(NA,1+S),1:kt,rep(NA,ka+kt*S*ka)),
			ID.age=c(rep(NA,1+S+kt),1:ka,rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S+kt+ka),1:(kt*S*ka)))

R <- kronecker(Pt,kronecker(diag(S),diag(ka))) + kronecker(diag(kt),kronecker(diag(S),Pa))
r.def <- S
A.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M7.FRS <- O ~ -1 + Intercept +
		    f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  	hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		    f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
		    f(ID.age, model="rw1", hyper=list(prec=list(prior=sdunif))) +
		    f(ID.year.area.age, model="generic0", Cmatrix=R, rankdef=r.def,
		  	constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  	extraconstr=list(A=A.constr, e=rep(0,S)))

M7.FRS <- inla(f.M7.FRS, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xt,Xa,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M7-SI: log(r_it) = intercept + phi_i + f(xt) + f(xa) + fs(xt,xa)
####################################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xt <- kronecker(Bt,kronecker(ones.S,ones.A))
Xa <- kronecker(ones.T,kronecker(ones.S,Ba))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt+ka+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,kt+ka+kt*S*ka)),
			ID.year=c(rep(NA,1+S),1:kt,rep(NA,ka+kt*S*ka)),
			ID.age=c(rep(NA,1+S+kt),1:ka,rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S+kt+ka),1:(kt*S*ka)))

R <- list(inla.as.sparse(kronecker(Pt,kronecker(diag(S),diag(ka)))),
	    inla.as.sparse(kronecker(diag(kt),kronecker(diag(S),Pa))))
r.def <- S
A.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M7.SI <- O ~ -1 + Intercept +
		    f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  	hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		    f(ID.year, model="rw1", hyper=list(prec=list(prior=sdunif))) +
		    f(ID.age, model="rw1", hyper=list(prec=list(prior=sdunif))) +
		    f(ID.year.area.age, model="generic3", Cmatrix=R, constr=TRUE,
		  	hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)),
		  	extraconstr=list(A=A.constr, e=rep(0,S)))

M7.SI <- inla(f.M7.SI, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xt,Xa,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M8-FRS: log(r_it) = intercept + phi_i + fs(xt) + fs(xa) + fs(xt,xa)
#######################################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xst <- kronecker(Bt,kronecker(diag(S),ones.A))
Xsa <- kronecker(ones.T,kronecker(diag(S),Ba))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+S*kt+S*ka+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,S*kt+S*ka+kt*S*ka)),
			ID.area.year=c(rep(NA,1+S),1:(S*kt),rep(NA,S*ka+kt*S*ka)),
			ID.area.age=c(rep(NA,1+S+S*kt),1:(S*ka),rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S+S*kt+S*ka),1:(kt*S*ka)))

R1 <- kronecker(Pt,diag(S))
r1.def <- S
A1.constr <- kronecker(t(ones.kt),diag(S))

R2 <- kronecker(diag(S),Pa)
r2.def <- S
A2.constr <- kronecker(diag(S),t(ones.ka))

R3 <- kronecker(Pt,kronecker(diag(S),diag(ka))) + kronecker(diag(kt),kronecker(diag(S),Pa))
r3.def <- S
A3.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M8.FRS <- O ~ -1 + Intercept +
		    f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  	hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		    f(ID.area.year, model="generic0", Cmatrix=R1, rankdef=r1.def,
		  	constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  	extraconstr=list(A=A1.constr, e=rep(0,S))) +
		    f(ID.area.age, model="generic0", Cmatrix=R2, rankdef=r2.def,
		  	constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  	extraconstr=list(A=A2.constr, e=rep(0,S))) +
		    f(ID.year.area.age, model="generic0", Cmatrix=R3, rankdef=r3.def,
		  	constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  	extraconstr=list(A=A3.constr, e=rep(0,S)))

M8.FRS <- inla(f.M8.FRS, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xst,Xsa,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M8-SI: log(r_it) = intercept + phi_i + fs(xt) + fs(xa) + fs(xt,xa)
#######################################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xst <- kronecker(Bt,kronecker(diag(S),ones.A))
Xsa <- kronecker(ones.T,kronecker(diag(S),Ba))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+S*kt+S*ka+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,S*kt+S*ka+kt*S*ka)),
			ID.area.year=c(rep(NA,1+S),1:(S*kt),rep(NA,S*ka+kt*S*ka)),
			ID.area.age=c(rep(NA,1+S+S*kt),1:(S*ka),rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S+S*kt+S*ka),1:(kt*S*ka)))

R1 <- kronecker(Pt,diag(S))
r1.def <- S
A1.constr <- kronecker(t(ones.kt),diag(S))

R2 <- kronecker(diag(S),Pa)
r2.def <- S
A2.constr <- kronecker(diag(S),t(ones.ka))

R3 <- list(inla.as.sparse(kronecker(Pt,kronecker(diag(S),diag(ka)))),
	    inla.as.sparse(kronecker(diag(kt),kronecker(diag(S),Pa))))
r3.def <- S
A3.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M8.SI <- O ~ -1 + Intercept +
		    f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  	hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		    f(ID.area.year, model="generic0", Cmatrix=R1, rankdef=r1.def,
		  	constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  	extraconstr=list(A=A1.constr, e=rep(0,S))) +
		    f(ID.area.age, model="generic0", Cmatrix=R2, rankdef=r2.def,
		  	constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  	extraconstr=list(A=A2.constr, e=rep(0,S))) +
		    f(ID.year.area.age, model="generic3", Cmatrix=R3, constr=TRUE,
		  	hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)),
		  	extraconstr=list(A=A3.constr, e=rep(0,S)))

M8.SI <- inla(f.M8.SI, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xst,Xsa,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M9-FRS: log(r_it) = intercept + phi_i + f(xt,xa) + fs(xt,xa)
################################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xta <- kronecker(Bt,kronecker(ones.S,Ba))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt*ka+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,kt*ka+kt*S*ka)),
			ID.year.age=c(rep(NA,1+S),1:(kt*ka),rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S+kt*ka),1:(kt*S*ka)))

R1 <- kronecker(Pt,diag(ka)) + kronecker(diag(kt),Pa)

R2 <- kronecker(Pt,kronecker(diag(S),diag(ka))) + kronecker(diag(kt),kronecker(diag(S),Pa))
r.def <- S
A.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M9.FRS <- O ~ -1 + Intercept +
		    f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  	hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		    f(ID.year.age, model="generic0", Cmatrix=R1, constr=TRUE,
		  	hyper=list(prec=list(prior=sdunif))) +
		    f(ID.year.area.age, model="generic0", Cmatrix=R2, rankdef=r.def,
		  	constr=TRUE, hyper=list(prec=list(prior=sdunif)),
		  	extraconstr=list(A=A.constr, e=rep(0,S)))

M9.FRS <- inla(f.M9.FRS, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xta,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## M9-SI: log(r_it) = intercept + phi_i + f(xt,xa) + fs(xt,xa)
###############################################################
Xs <- kronecker(ones.T,kronecker(diag(S),ones.A))
Xta <- kronecker(Bt,kronecker(ones.S,Ba))
Xtsa <- kronecker(Bt,kronecker(diag(S),Ba))

Data.inla <- list(O=Data$Obs, E=Data$Pop,
			Intercept=c(1,rep(NA,S+kt*ka+kt*S*ka)),
			ID.area=c(NA,1:S,rep(NA,kt*ka+kt*S*ka)),
			ID.year.age=c(rep(NA,1+S),1:(kt*ka),rep(NA,kt*S*ka)),
			ID.year.area.age=c(rep(NA,1+S+kt*ka),1:(kt*S*ka)))

R1 <- list(inla.as.sparse(kronecker(Pt,diag(ka))),
	     inla.as.sparse(kronecker(diag(kt),Pa)))

R2 <- list(inla.as.sparse(kronecker(Pt,kronecker(diag(S),diag(ka)))),
	    inla.as.sparse(kronecker(diag(kt),kronecker(diag(S),Pa))))

A.constr <- kronecker(t(ones.kt),kronecker(diag(S),t(ones.ka)))

f.M9.SI <- O ~ -1 + Intercept +
		    f(ID.area, model="generic1", Cmatrix=Q.Leroux, constr=TRUE,
		  	hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif))) +
		    f(ID.year.age, model="generic3", Cmatrix=R1, constr=TRUE,
		 	hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif))) + 
		    f(ID.year.area.age, model="generic3", Cmatrix=R2, constr=TRUE,
		  	hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)),
		  	extraconstr=list(A=A.constr, e=rep(0,S))) 

M9.SI <- inla(f.M9.SI, family="poisson", data=Data.inla, E=E,
		   control.predictor=list(compute=TRUE, A=cbind(rep(1,T*S*A),Xs,Xta,Xtsa), link=1, cdf=c(log(1))),
		   control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
		   control.inla=list(strategy="simplified.laplace", verbose=TRUE))


## Save INLA models ##
MODELS <- list(M1=M1, M2=M2, M3=M3, M4=M4,
               M5.FRS=M5.FRS, M5.SI=M5.SI,
               M6.FRS=M5.FRS, M6.SI=M5.SI,
               M7.FRS=M5.FRS, M7.SI=M5.SI,
               M8.FRS=M5.FRS, M8.SI=M5.SI,
               M9.FRS=M5.FRS, M9.SI=M5.SI)
               
save(MODELS, file="Psplines_INLA.Rdata")
