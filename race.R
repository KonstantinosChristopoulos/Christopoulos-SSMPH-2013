#By race
####set up####
library(rethinking)
D<-read.csv(".../data.csv")

D$SIstd<-standardize(D$SI)

D$log_pop_hisp<-log(D$pop_hispanic)
D$log_pop_white<-log(D$pop_white)
D$log_pop_AIAN<-log(D$pop_AIAN)
D$log_pop_NHPI<-log(D$pop_NHPI)
D$log_pop_Asian<-log(D$pop_Asian)
D$log_pop_black<-log(D$pop_black)

D$LH<-standardize(D$Lr100_hisp)
D$LW<-standardize(D$Lr100_white)
D$LNHPI<-standardize(D$Lr100_NHPI)
D$LAIAN<-standardize(D$Lr100_AIAN)
D$LA<-standardize(D$Lr100_Asian)
D$LB<-standardize(D$Lr100_black)

D$rate_white<-(D$d_white*100000)/D$pop_white
D$rate_black<-(D$d_black*100000)/D$pop_black
D$rate_Asian<-(D$d_Asian*100000)/D$pop_Asian
D$rate_AIAN<-(D$d_AIAN*100000)/D$pop_AIAN
D$rate_NHPI<-(D$d_NHPI*100000)/D$pop_NHPI
D$rate_hisp<-(D$d_hisp*100000)/D$pop_hisp

Dmat=read.csv(".../Dmatcol.csv")
rownames(Dmat)<-c("AL", "AK" ,  "AZ" ,  "AR" ,  "CA" ,  "CO" ,  "CT" ,  "DE" ,  "DC" ,  "FL" ,  "GA" ,  "HI" ,  "ID"  , "IL" ,  "IN" ,  "IA" ,  "KS" ,  "KY" ,  "LA" ,  "ME" ,  "MD",   "MA" ,  "MI"  , "MN" ,  "MS" ,  "MO" ,  "MT"  , "NE"  , "NV"  , "NH" ,  "NJ"  , "NM"  , "NY" ,  "NC" ,  "ND"  , "OH" ,"OK"  , "OR" ,  "PA"  , "RI" ,  "SC"  , "SD" ,  "TN"  , "TX" ,  "UT" ,  "VT" ,  "VA" ,  "WA" ,  "WV" ,  "WI"  , "WY")
Dmat=as.matrix(Dmat)

d<-list(
	SI=D$SIstd,
	Sw=D$d_white,
	Sb=D$d_black,
	Sh=D$d_hisp,
	Saian=D$d_AIAN,
	Snhpi=D$d_NHPI,
	Sa=D$d_Asian,
	Pw=D$log_pop_white,
	Pb=D$log_pop_black,
	Ph=D$log_pop_hisp,
	Paian=D$log_pop_AIAN,
	Pnhpi=D$log_pop_NHPI,
	Pa=D$log_pop_Asian,
	Lw=D$LW,
	Lb=D$LB,
	Lh=D$LH,
	La=D$LA,
	Laian=D$LAIAN,
	Lnhpi=D$LNHPI,
	LaianR=D$Lr100_AIAN,
	LnhpiR=D$Lr100_NHPI,
	lPaian=D$log_pop_AIAN,
    lPnhpi=D$log_pop_NHPI,
	state=D$id,
	Dmat=Dmat
	)

#####models#######
############################WHITE######################################
#spatial
mWs<-ulam(
	alist(
		Sw~dpois(lambda),
		log(lambda)<-Pw + a + k[state] + b1*SI +b2*Lw,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mWs)
postW<-extract.samples(mWs)

#Gamma-Poisson
mWg<-ulam(
	alist(
		Sw~dgampois(lambda, phi),
		log(lambda)<-Pw + a + b1*SI +b2*Lw,
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8,1),
		phi~dexp(1)		
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mWg)
postWg<-extract.samples(mWg)

#################################BLACK#################################################
mBs<-ulam(
	alist(
		Sb~dpois(lambda),
		log(lambda)<-Pb + a + k[state] + b1*SI +b2*Lb,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8.5,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mBs)
postB<-extract.samples(mBs)

#Gamma-Poisson
mBg<-ulam(
	alist(
		Sb~dgampois(lambda, phi),
		log(lambda)<-Pb + a + b1*SI +b2*Lb,
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8.5,1),
		phi~dexp(1)		
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mBg)
postBg<-extract.samples(mBg)

########################################HISPANIC########################################
#spatial
mHs<-ulam(
	alist(
		Sb~dpois(lambda),
		log(lambda)<-Ph + a + k[state] + b1*SI +b2*Lh,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8.5,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mHs)
postHs<-extract.samples(mHs)


#Gamma-Poisson
mHg<-ulam(
	alist(
		Sb~dgampois(lambda, phi),
		log(lambda)<-Ph + a + b1*SI +b2*Lh,
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8.5,1),
		phi~dexp(1)		
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mHg)
postHg<-extract.samples(mHg)

######################################ASIAN#############################################
mAs<-ulam(
	alist(
		Sa~dpois(lambda),
		log(lambda)<-Pa + a + k[state] + b1*SI +b2*La,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8.5,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mAs)
postA<-extract.samples(mAs)

#Gamma-Poisson
mAg<-ulam(
	alist(
		Sa~dgampois(lambda, phi),
		log(lambda)<-Pa + a + b1*SI +b2*La,
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8.5,1),
		phi~dexp(1)		
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mAg)
postAg<-extract.samples(mAg)

############################################AIAN#######################################
#spatial
mAIANs<-ulam(
	alist(
		Sa~dpois(lambda),
		log(lambda)<-Paian + a + k[state] + b1*SI +b2*Laian,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8,1)
	),data=d, chains=4 ,cores=2, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mAIANs)
postAIANs<-extract.samples(mAIANs)

#Gamma-Poisson
mAIANg<-ulam(
	alist(
		Sa~dgampois(lambda, phi),
		log(lambda)<-Paian + a + b1*SI +b2*Laian,
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8,1),
		phi~dexp(1)		
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mAIANg)
postAIANg<-extract.samples(mAIANg)
	
#zero inflated Poisson
mAIANzi<-ulam(
	alist(
		Sa~dzipois(p, lambda),
		logit(p)<-c + c1*Laian + c2*lPaian, 
		log(lambda)<-Paian + a + b1*SI +b2*Laian,
		c(b1,b2,c1)~dnorm(0,1),
		c2~dnorm(-1,0.5),
		a~dnorm(-8,1),
		c~dnorm(-2,1)
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mAIANzi)
postAIANzi<-extract.samples(mAIANzi)

##############################################NHPI#####################################
#spatial
mNHPIs<-ulam(
	alist(
		Sa~dpois(lambda),
		log(lambda)<-Pnhpi + a + k[state] + b1*SI +b2*Lnhpi,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-9,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mNHPIs)
postNHPIs<-extract.samples(mNHPIs)

#Gamma-Poisson 
mNHPIg<-ulam(
	alist(
		Sa~dgampois(lambda, phi),
		log(lambda)<-Pnhpi + a + b1*SI +b2*Lnhpi,
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-9,1),
		phi~dexp(1)		
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mNHPIg)
postNHPIg<-extract.samples(mNHPIg)
	
#zero-inflated Poisson
mNHPIzi<-ulam(
	alist(
		Sa~dzipois(p, lambda),
		logit(p)<-c + c1*Lnhpi + c2*lPnhpi, 
		log(lambda)<-Pnhpi + a + b1*SI +b2*Lnhpi,
		c(b1,b2,c1)~dnorm(0,1),
		c2~dnorm(-1,0.5),
		a~dnorm(-8,1),
		c~dnorm(0,1)
	),data=d, chains=4 ,cores=4, log_lik=TRUE)		

precis(mNHPIzi)
postNHPIzi<-extract.samples(mNHPIzi)	
