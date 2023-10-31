#By gender
####set up####
library(rethinking)

D<-read.csv(".../data.csv")
D$SIstd<-standardize(D$SI)
D$log_pop_fem<-log(D$pop_fem)
D$log_pop_male<-log(D$pop_male)
D$LM<-standardize(D$Lr100_male)
D$LF<-standardize(D$Lr100_fem)
D$rate_fem<-(D$d_fem*100000)/D$pop_fem
D$rate_male<-(D$d_male*100000)/D$pop_male

Dmat=read.csv(".../Dmatcol.csv")
rownames(Dmat)<-c("AL", "AK" ,  "AZ" ,  "AR" ,  "CA" ,  "CO" ,  "CT" ,  "DE" ,  "DC" ,  "FL" ,  "GA" ,  "HI" ,  "ID"  , "IL" ,  "IN" ,  "IA" ,  "KS" ,  "KY" ,  "LA" ,  "ME" ,  "MD",   "MA" ,  "MI"  , "MN" ,  "MS" ,  "MO" ,  "MT"  , "NE"  , "NV"  , "NH" ,  "NJ"  , "NM"  , "NY" ,  "NC" ,  "ND"  , "OH" ,"OK"  , "OR" ,  "PA"  , "RI" ,  "SC"  , "SD" ,  "TN"  , "TX" ,  "UT" ,  "VT" ,  "VA" ,  "WA" ,  "WV" ,  "WI"  , "WY")
Dmat=as.matrix(Dmat)

d<-list(
	SI=D$SIstd,
	Sf=D$d_fem,
	Sm=D$d_male,
	Lf=D$LF,
	Lm=D$LM,
	Pm=D$log_pop_male,
	Pf=D$log_pop_fem,
	state=D$id,
	Dmat=Dmat
	)

####models###########
#MALE
#spatial
mMs<-ulam(
	alist(
		Sm~dpois(lambda),
		log(lambda)<-Pm + a + k[state] + b1*SI +b2*Lm,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		b1~dnorm(0,1),
		b2~dnorm(0,1),
		a~dnorm(-8,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mMs)
postM<-extract.samples(mMs)

#Gamma-Poisson
mMg<-ulam(
	alist(
		Sm~dgampois(lambda, phi),
		log(lambda)<-Pm + a + b1*SI +b2*Lm,
		b1~dnorm(0,1),
		b2~dnorm(0,1),
		phi~dexp(1),
		a~dnorm(-8,1)
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mMg)
postMg<-extract.samples(mMg)



#FEMALE##################
#spatial
mFs<-ulam(
	alist(
		Sf~dpois(lambda),
		log(lambda)<-Pf + a + k[state] + b1*SI +b2*Lf,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(2),
		rhosq~dexp(0.5),
		b1~dnorm(0,1),
		b2~dnorm(0,1),
		a~dnorm(-9,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95),  iter=4000, warmup=1000, log_lik=TRUE)

precis(mFs)
postF<-extract.samples(mFs)

#Gamma-Poisson
mFg<-ulam(
	alist(
		Sf~dgampois(lambda, phi),
		log(lambda)<-Pf + a + b1*SI +b2*Lf,
		b1~dnorm(0,1),
		b2~dnorm(0,1),
		phi~dexp(1),
		a~dnorm(-9,1)
	),data=d, chains=4 ,cores=4, log_lik=TRUE)

precis(mFg)
postFg<-extract.samples(mFg)

