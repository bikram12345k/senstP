## DATA lssinc07.csv downloaded from www.rerf.org

## Load data
lss07 = read.csv("lssinc07.csv", header=T)

## Only use the participants NIC or in city with exposure
##	between 200milligray and 4000milligray in Hiroshima
lss07 = lss07[(lss07$un4gy!=0) & lss07$distcat!=2,]
lss07 = rbind(lss07[(lss07$distcat==1) & lss07$dcat>=12,],
			lss07[(lss07$distcat==3),])
lss07hir = lss07[lss07$city==1,]


## Create Table 1
lss07hir$ageatexp = lss07hir$agxcat
lss07hir$ageatexp = 1*(lss07hir$ageatexp <=4 ) + 
	 1*(lss07hir$ageatexp <= 13) + 
	1*(lss07hir$ageatexp>0)

lss07hir$str = as.factor(lss07hir$ageatexp):as.factor(lss07hir$sex)

library(plyr)
tab = ddply( lss07hir, .(sex, ageatexp, distcat),
	summarise,
	subjects = sum(subjects), solid = sum(solid) )

for(i in 1:(nrow(tab)/2)){
	idx = (i-1)*2 + 1:2
	temp = tab[idx, 4:5]

	odd1 = temp[1,2]/(temp[1,1] - temp[1,2])
	odd2 = temp[2,2]/(temp[2,1] - temp[2,2])

	print( (temp[1,2]/temp[1,1])/(temp[2,2]/temp[2,1]))
}	

## Table 1
cbind(tab, apply(tab, 1, function(x) x[5]/x[4]*100))


##  Stratified data on age-at-exposure and sex
strData <- ddply(lss07hir, .(sex, agxcat, distcat),
				summarize,
			csolid = sum(solid), nocsolid = sum(subjects)-sum(solid))


## Data for the 'All participants' comparison
strAna <-  array(NA, c(2,2,nrow(strData)/2))

for(i in 1:(dim(strAna)[3])){
	idx = 0+(i-1)*2+(1:2)

	strAna[1,1,i] = strData[idx[1], 4]	
	strAna[1,2,i] = strData[idx[1], 5]	

	strAna[2,1,i] = strData[idx[2], 4]	
	strAna[2,2,i] = strData[idx[2], 5]	
}

#cbind(1+seq(0,.45,.01), round(sapply(1+seq(0,.45,.01), function(gam) mh(strAna, Gamma=gam)$pval),4))


### Data in the six subgroups.
y = 8		## <= 20 years male
yp = 26		## <= 65 years male
z = 38		## <= 20 years female
zp = 56		## <= 65 years female
## Male
## Young
## rows 1--8
strAna1 <- array(NA, c(2,2,y/2))

for(i in 1:(dim(strAna1)[3])){
	idx = 0+(i-1)*2+(1:2)

	strAna1[1,1,i] = strData[idx[1], 4]	
	strAna1[1,2,i] = strData[idx[1], 5]	

	strAna1[2,1,i] = strData[idx[2], 4]	
	strAna1[2,2,i] = strData[idx[2], 5]	
}

## Male
## Adult
## rows 9--26
strAna2 <- array(NA, c(2,2,(yp-y)/2))

for(i in 1:(dim(strAna2)[3])){
	idx = y+(i-1)*2+(1:2)

	strAna2[1,1,i] = strData[idx[1], 4]	
	strAna2[1,2,i] = strData[idx[1], 5]	

	strAna2[2,1,i] = strData[idx[2], 4]	
	strAna2[2,2,i] = strData[idx[2], 5]	
}


## Male
## Old
## rows 27--30
strAna3 <- array(NA, c(2,2,(30-yp)/2))

for(i in 1:(dim(strAna3)[3])){
	idx = yp+(i-1)*2+(1:2)

	strAna3[1,1,i] = strData[idx[1], 4]	
	strAna3[1,2,i] = strData[idx[1], 5]	

	strAna3[2,1,i] = strData[idx[2], 4]	
	strAna3[2,2,i] = strData[idx[2], 5]	
}


## Female
## Young
## rows 31--38
strAna4 <- array(NA, c(2,2,(z-31+1)/2))

for(i in 1:(dim(strAna4)[3])){
	idx = 30+(i-1)*2+(1:2)

	strAna4[1,1,i] = strData[idx[1], 4]	
	strAna4[1,2,i] = strData[idx[1], 5]	

	strAna4[2,1,i] = strData[idx[2], 4]	
	strAna4[2,2,i] = strData[idx[2], 5]	
}

## Female
## Adult
## rows 39--56
strAna5 <- array(NA, c(2,2,(zp-z)/2))

for(i in 1:(dim(strAna5)[3])){
	idx = z+(i-1)*2+(1:2)

	strAna5[1,1,i] = strData[idx[1], 4]	
	strAna5[1,2,i] = strData[idx[1], 5]	

	strAna5[2,1,i] = strData[idx[2], 4]	
	strAna5[2,2,i] = strData[idx[2], 5]	
}


## Female
## Old
## rows 57--60
strAna6 <- array(NA, c(2,2,(60-zp)/2))

for(i in 1:(dim(strAna6)[3])){
	idx = zp+(i-1)*2+(1:2)

	strAna6[1,1,i] = strData[idx[1], 4]	
	strAna6[1,2,i] = strData[idx[1], 5]	

	strAna6[2,1,i] = strData[idx[2], 4]	
	strAna6[2,2,i] = strData[idx[2], 5]	
}

## Final analysis
library(sensitivity2x2xk)
library(struncatedP)
library(partitions)

gamSeq = c(1,  1.2, 1.4, 1.45, 1.5, 1.55)

res = matrix(NA, 0, 6)

cbind(1+seq(0,.45,.01), round(sapply(1+seq(0,.45,.01), function(gam) mh(strAna1, Gamma=gam)$pval),4))

for(i in 1:length(gamSeq)){
gam = gamSeq[i]
p <- c(mh(strAna1, Gamma=gam)$pval,
mh(strAna2, Gamma=gam)$pval,
mh(strAna3, Gamma=gam)$pval,

mh(strAna4, Gamma=gam)$pval,
mh(strAna5, Gamma=gam)$pval,
mh(strAna6, Gamma=gam)$pval)

p = abs(pmin(1,p))

xiVIp <- function (x, s = 1/2) 
    s * exp(- 10 * x)


res = rbind(res, 
c(mh(strAna, Gamma=gam)$pval,
6*min(sort(p)/(1:6)),
truncatedP(p, trunc=1),
truncatedP(p, trunc=.1),
getstruncatedPval1(p, foo(xiVIp, trunc=.1), N=25, trunc=.1),
getstruncatedPval1(p, foo(xiVIpf, trunc=.1), N=25, trunc=.1)))

}

rownames(res) = gamSeq
colnames(res) = c('All', 'Simes', 'Fisher', 'truncated', 'struncated')

round(res, 4)

## Use the following 

gamSeq = c(1,  1.2, 1.4, 1.45, 1.5, 1.55)

res = matrix(NA, 0, 5)

cbind(1+seq(0,.45,.01), round(sapply(1+seq(0,.45,.01), function(gam) mh(strAna1, Gamma=gam)$pval),4))

for(i in 1:length(gamSeq)){
gam = gamSeq[i]
p <- c(mh(strAna1, Gamma=gam)$pval,
mh(strAna2, Gamma=gam)$pval,
mh(strAna3, Gamma=gam)$pval,

mh(strAna4, Gamma=gam)$pval,
mh(strAna5, Gamma=gam)$pval,
mh(strAna6, Gamma=gam)$pval)

p = abs(pmin(1,p))

xiVIp <- function (x, s = 1/2) 
    s * exp(- 10 * x)

foo <- function(f, trunc=.15)
			function(x) ifelse(x<.5,f((x-trunc)/(.5-trunc)),0)


res = rbind(res, 
c(mh(strAna, Gamma=gam)$pval,
6*min(sort(p)/(1:6)),
truncatedP(p, trunc=1),
truncatedP(p, .1),
getstruncatedPval1(p, foo(xiVp, trunc=.15), N=20, trunc=.15)#,
#getstruncatedPval1(p, foo(xiVIpf, trunc=.1), N=25, trunc=.1))
))

}

rownames(res) = gamSeq
colnames(res) = c('All', 'Simes', 'Fisher', 'truncated', 'struncated')

round(res, 4)
