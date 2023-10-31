#Read README file
#Total population 
library(rethinking)
D<-read.csv(".../data.csv") #download from "Data" folder add the path in your computer
str(D)

#Some data transformation for the analysis and the graphs
D$SIstd<-standardize(D$SI) #standardization of the stringency index (SI) variable
D$LT<-standardize(D$Lr100_tot) #standardization of the lagged rates (LT) variable
D$log_pop_tot<-log(D$pop_tot) #create the population offset
D$rate_tot<-(D$d_tot*100000)/D$pop_tot #create rates per 100k for the plot

#Spatial distance matrix.
Dmat=read.csv(".../Dmatcol.csv") # download from "Data" folder add the path in your computer
rownames(Dmat)<-c("AL", "AK" ,  "AZ" ,  "AR" ,  "CA" ,  "CO" ,  "CT" ,  "DE" ,  "DC" ,  "FL" ,  "GA" ,  "HI" ,  "ID"  , "IL" ,  "IN" ,  "IA" ,  "KS" ,  "KY" ,  "LA" ,  "ME" ,  "MD" ,  "MA" ,  "MI" , "MN" ,  "MS" ,  "MO" ,  "MT"  , "NE"  , "NV"  , "NH" ,  "NJ"  , "NM" , "NY" ,  "NC" ,  "ND"  , "OH" ,"OK"  , "OR" ,  "PA"  , "RI" ,  "SC"  , "SD" ,  "TN" , "TX" ,  "UT" ,  "VT" ,  "VA" ,  "WA" ,  "WV" ,  "WI"  , "WY")
Dmat=as.matrix(Dmat)

#Make list for the analysis
d<-list(
	S=D$d_tot,
	SI=D$SIstd,
	P=D$log_pop_tot,
	LS=D$LT,
	state=D$id,
	Dmat=Dmat
	)

#Model 1. Gaussian process for the spatial model
mTs<-ulam(
	alist(
		S~dpois(lambda),
		log(lambda)<-P + a + k[state] + b1*SI +b2*LS,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		b1~dnorm(0,1),
		b2~dnorm(0,1),
		a~dnorm(-8,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)

precis(mTs,2)
posts<-extract.samples(mTs)

#Model 2. Gamma-Poisson model 
mTg<-ulam(
	alist(
		S~dgampois(lambda, phi),
		log(lambda)<-P + a + b1*SI +b2*LS,
		b1~dnorm(0,1),
		b2~dnorm(0,1),
		phi~dexp(1),
		a~dnorm(-8,1)
	),data=d, chains=4 ,cores=4, log_lik=TRUE)
	
precis(mTg)
postg<-extract.samples(mTg)


