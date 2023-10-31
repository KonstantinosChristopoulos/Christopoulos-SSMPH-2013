#Code for the figures 2-5. For Figure 1 see the maps.r file.
#Run the corresponding R files (total,gender,race) first
library(rethinking)
D<-read.csv(".../data.csv")

D$SIstd<-standardize(D$SI)

D$LT<-standardize(D$Lr100_tot)
D$log_pop_tot<-log(D$pop_tot)

D$log_pop_fem<-log(D$pop_fem)
D$log_pop_male<-log(D$pop_male)
D$LM<-standardize(D$Lr100_male)
D$LF<-standardize(D$Lr100_fem)

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

#Figure 2. ##########################
#load distance matrix
Dmat=read.csv(".../Dmatcol.csv")
rownames(Dmat)<-c("AL", "AK" ,  "AZ" ,  "AR" ,  "CA" ,  "CO" ,  "CT" ,  "DE" ,  "DC" ,  "FL" ,  "GA" ,  "HI" ,  "ID"  , "IL" ,  "IN" ,  "IA" ,  "KS" ,  "KY" ,  "LA" ,  "ME" ,  "MD" ,  "MA" ,  "MI" , "MN" ,  "MS" ,  "MO" ,  "MT"  , "NE"  , "NV"  , "NH" ,  "NJ"  , "NM" , "NY" ,  "NC" ,  "ND"  , "OH" ,"OK"  , "OR" ,  "PA"  , "RI" ,  "SC"  , "SD" ,  "TN" , "TX" ,  "UT" ,  "VT" ,  "VA" ,  "WA" ,  "WV" ,  "WI"  , "WY")
Dmat=as.matrix(Dmat)

d<-list(
	S=D$d_tot,
	SI=D$SIstd,
	P=D$log_pop_tot,
	LS=D$LT,
	Sf=D$d_fem,
	Sm=D$d_male,
	Lf=D$LF,
	Lm=D$LM,
	Pm=D$log_pop_male,
	Pf=D$log_pop_fem,
	Sw=D$d_white,
	Sb=D$d_black,
	Sh=D$d_hisp,
	Saian=D$d_AIAN,
	Sa=D$d_Asian,
	Pw=D$log_pop_white,
	Pb=D$log_pop_black,
	Ph=D$log_pop_hisp,
	Paian=D$log_pop_AIAN,
	Pa=D$log_pop_Asian,
	Lw=D$LW,
	Lb=D$LB,
	Lh=D$LH,
	La=D$LA,
	Laian=D$LAIAN,
	state=D$id,
	Dmat=Dmat
	)

#### Gaussian process models####
#TOTAL
Total<-ulam(
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

#MALE
Male<-ulam(
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

#FEMALE
Female<-ulam(
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

#ASIAN
Asian<-ulam(
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



#AIAN
AIAN<-ulam(
	alist(
		Sa~dpois(lambda),
		log(lambda)<-Paian + a + k[state] + b1*SI +b2*Laian,
		vector[51]:k~multi_normal(0,SIGMA),
		matrix[51,51]:SIGMA<-cov_GPL2(Dmat, etasq, rhosq, 0.01 ),
		etasq~dexp(3),
		rhosq~dexp(1),
		c(b1,b2)~dnorm(0,1),
		a~dnorm(-8,1)
	),data=d, chains=4 ,cores=4, control=list(adapt_delta=0.95), iter=4000, warmup=1000, log_lik=TRUE)



#BLACK
Black<-ulam(
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


#HISPANIC
Hispanic<-ulam(
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


#WHITE
White<-ulam(
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


#extract samples from posterior
postT<-extract.samples(Total)
postM<-extract.samples(Male)
postF<-extract.samples(Female)
postA<-extract.samples(Asian)
postAIAN<-extract.samples(AIAN)
postB<-extract.samples(Black)
postH<-extract.samples(Hispanic)
postW<-extract.samples(White)

plot<-list(
	Total=postT$b1,
	Male=postM$b1,
	Female=postF$b1,
		Asian=postA$b1,
		AIAN=postAIAN$b1,
		Black=postB$b1,
		Hispanic=postH$b1,
		White=postW$b1
		)

plot(precis(plot))

#Figure 3. ##########################Total
#Run total.r first
#Figure 3a posterior densities
dens(posts$b1)
dens(postg$b1, col=rangi2, add=TRUE)
abline(v=0, lty=2)

#Figure 3b (Posterior prediction from Model 1)
lambda<-function(SI) exp(posts$a + posts$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
mu<-sapply(SI.seq,lambda)
lmu<-apply(mu,2,mean)
lci<-apply(mu,2,PI)
lmu1k<-lmu*100000
lci1k<-lci*100000
k<-PSIS(mTs, pointwise=TRUE) $k

plot(d$SI,D$rate_tot, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=col.alpha(rangi2), pch=1 ,lwd=2, ylim=c(0,60) , cex=1+normalize(k))
lines(SI.seq, lmu1k ,lty=1 ,lwd=1.5)
shade(lci1k, SI.seq ,xpd=TRUE)

#Figure 3c (Posterior prediction from Model 2)
lambda<-function(SI) exp(postg$a + postg$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
mu<-sapply(SI.seq,lambda)
lmu<-apply(mu,2,mean)
lci<-apply(mu,2,PI)
lmu1k<-lmu*100000
lci1k<-lci*100000
k<-PSIS(mTg, pointwise=TRUE) $k

plot(d$SI,D$rate_tot, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=col.alpha(rangi2), pch=1 ,lwd=2, ylim=c(0,60) , cex=1+normalize(k))
lines(SI.seq, lmu1k ,lty=1 ,lwd=1.5)
shade(lci1k, SI.seq ,xpd=TRUE)

#Figure 4. ##########################Gender
#Run gender.r first
#Figure 4a posterior densities
dens(postM$b1)
dens(postF$b1, col=rangi2, add=TRUE)
abline(v=0, lty=2)

#Figure 4b Posterior prediction from Model 1 for males
postM<-extract.samples(mMs)
lambda<-function(SI) exp(postM$a + postM$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
mu<-sapply(SI.seq,lambda)
lmu<-apply(mu,2,mean)
lci<-apply(mu,2,PI)
lmu1k<-lmu*100000
lci1k<-lci*100000
k<-PSIS(mMs, pointwise=TRUE) $k

plot(d$SI,D$rate_male, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=rangi2, pch=1 ,lwd=2, ylim=c(0,100) , cex=1+normalize(k))

lines(SI.seq, lmu1k ,lty=1 ,lwd=1.5)
shade(lci1k, SI.seq ,xpd=TRUE)

identify(x= d$SI ,y= D$rate_male, labels= D$Abv ,cex=0.8)
dev.print(pdf,'plot_male_s.pdf')

#Figure 4c Posterior prediction from Model 1 for females
lambda<-function(SI) exp(postF$a + postF$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
mu<-sapply(SI.seq,lambda)
lmu<-apply(mu,2,mean)
lci<-apply(mu,2,PI)
lmu1k<-lmu*100000
lci1k<-lci*100000
k<-PSIS(mFs, pointwise=TRUE) $k

plot(d$SI,D$rate_fem, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=rangi2, pch=1 ,lwd=2, ylim=c(0,100) , cex=1+normalize(k))

lines(SI.seq, lmu1k ,lty=1 ,lwd=1.5)
shade(lci1k, SI.seq ,xpd=TRUE)

identify(x= d$SI ,y= D$rate_fem , labels= D$Abv ,cex=0.8)

#Figure 5. ##########################Race
#Run race.r first
#Figure 5a Posterior prediction for Asian
lambda<-function(SI) exp(postA$a + postA$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
mu<-sapply(SI.seq,lambda)
lmu<-apply(mu,2,mean)
lci<-apply(mu,2,PI)
lmu1k<-lmu*100000
lci1k<-lci*100000
k<-PSIS(mAs, pointwise=TRUE) $k

plot(d$SI,D$rate_Asian, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=rangi2, pch=1 ,lwd=2, ylim=c(0,70) , cex=1+normalize(k))
lines(SI.seq, lmu1k ,lty=1 ,lwd=1.5)
shade(lci1k, SI.seq ,xpd=TRUE)
identify(x= d$SI ,y= D$rate_Asian, labels= D$Abv)

#Figure 5b Posterior prediction for American Indian or Alaskan Native
lambdaS<-function(SI) exp(postAIANs$a + postAIANs$b1*SI)
lambdaG<-function(SI) exp(postAIANg$a + postAIANg$b1*SI)
lambdaZI<-function(SI) exp(postAIANzi$a + postAIANzi$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)

#spatial
muS<-sapply(SI.seq,lambdaS)
lmuS<-apply(muS,2,mean)
lciS<-apply(muS,2,PI)
lmu1kS<-lmuS*100000
lci1kS<-lciS*100000

#Gamma-Poisson
muG<-sapply(SI.seq,lambdaG)
lmuG<-apply(muG,2,mean)
lciG<-apply(muG,2,PI)
lmu1kG<-lmuG*100000
lci1kG<-lciG*100000

#ZI Poisson
muZI<-sapply(SI.seq,lambdaZI)
lmuZI<-apply(muZI,2,mean)
lciZI<-apply(muZI,2,PI)
lmu1kZI<-lmuZI*100000
lci1kZI<-lciZI*100000

#plot
plot(d$SI,D$rate_AIAN, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=rangi2, pch=1 ,lwd=2, ylim=c(0,120))

#spatial
lines(SI.seq, lmu1kS ,lty=1 ,lwd=1.5)
shade(lci1kS, SI.seq ,xpd=TRUE)
#Gamma-Poisson
lines(SI.seq, lmu1kG ,lty=2 ,lwd=1.5)
shade(lci1kG, SI.seq ,xpd=TRUE)
#ZI Poisson
lines(SI.seq, lmu1kZI ,lty=3 ,lwd=1.5)
shade(lci1kZI, SI.seq ,xpd=TRUE)

identify(x= d$SI ,y= D$rate_AIAN, labels= D$Abv)


#Figure 5c Posterior prediction for Black or African American
lambda<-function(SI) exp(postB$a + postB$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
mu<-sapply(SI.seq,lambda)
lmu<-apply(mu,2,mean)
lci<-apply(mu,2,PI)
lmu1k<-lmu*100000
lci1k<-lci*100000
k<-PSIS(mBs, pointwise=TRUE) $k

plot(d$SI,D$rate_black, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=rangi2, pch=1 ,lwd=2, ylim=c(0,45) , cex=1+normalize(k))
lines(SI.seq, lmu1k ,lty=1 ,lwd=1.5)
shade(lci1k, SI.seq ,xpd=TRUE)
identify(x= d$SI ,y= D$rate_black, labels= D$Abv)

#Figure 5d Posterior prediction for Hispanic origin
lambdaG<-function(SI) exp(postHg$a + postHg$b1*SI)
lambdaS<-function(SI) exp(postHs$a + postHs$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
#spatial
muS<-sapply(SI.seq,lambdaS)
lmuS<-apply(muS,2,mean)
lciS<-apply(muS,2,PI)
lmu1kS<-lmuS*100000
lci1kS<-lciS*100000
#gamma-poisson
muG<-sapply(SI.seq,lambdaG)
lmuG<-apply(muG,2,mean)
lciG<-apply(muG,2,PI)
lmu1kG<-lmuG*100000
lci1kG<-lciG*100000

plot(d$SI,D$rate_hisp, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=rangi2, pch=1 ,lwd=2, ylim=c(0,40))
#spatial
lines(SI.seq, lmu1kS ,lty=1 ,lwd=1.5)
shade(lci1kS, SI.seq ,xpd=TRUE)
#gamma-poisson
lines(SI.seq, lmu1kG ,lty=2 ,lwd=1.5)
shade(lci1kG, SI.seq ,xpd=TRUE)

identify(x= d$SI ,y= D$rate_hisp, labels= D$Abv)

#Figure 5e Posterior prediction for White
lambda<-function(SI) exp(postW$a + postW$b1*SI)
ns<-100
SI.seq<-seq(from=-3, to=3 , length.out=ns)
mu<-sapply(SI.seq,lambda)
lmu<-apply(mu,2,mean)
lci<-apply(mu,2,PI)
lmu1k<-lmu*100000
lci1k<-lci*100000
k<-PSIS(mWs, pointwise=TRUE) $k

plot(d$SI,D$rate_white, xlab="Stringency Index (std)", ylab="Suicide deaths per 100k", 
	col=rangi2, pch=1 ,lwd=2, ylim=c(0,80) , cex=1+normalize(k))
lines(SI.seq, lmu1k ,lty=1 ,lwd=1.5)
shade(lci1k, SI.seq ,xpd=TRUE)
identify(x= d$SI ,y= D$rate_white , labels= D$Abv)




