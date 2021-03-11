args = commandArgs(trailingOnly=TRUE)
split=as.numeric(args)[1]

library(senstrat)
library(sensitivitymw)
library(sensitivitymv)
library(struncatedP)
library(partitions)
library(foreach)
library(doSNOW)
library(poolr)
#library(beepr)

set.seed(50*split)

cl<-makeCluster(6) 
registerDoSNOW(cl)


I = 400

## Choose ONE
theta = c(0, 0, 0, 0.15, 0.3, 0.6)	# Alternative
#theta = rep(0, 6)	# NULL
#theta = rep(0.15, 6)	# NULL


K = length(theta)
p = c(20, 20, 20, 20, 10)
gamSeq = c(seq(1, 2.2, .1))#[1:6]

Itr = 100


Simes = Fisher = tP = tP2 = SumP = stouf= minP = array(NA, c(Itr, K, length(gamSeq)))
stPI = stPII = stPIII = stPIV = stPV = stPVI = array(NA, c(Itr, K, length(gamSeq)))
stPIf = stPIIf = stPIIIf = stPIVf = stPVf = stPVIf = array(NA, c(Itr, K, length(gamSeq)))

temp = Sys.time()
itr=1
Itr1 = Itr
while(itr<= Itr1){
#for(itr in 1:Itr){
		

	#theta = rep(0, rbinom(1,4,.25))
	#theta = c(theta, rnorm(K-length(theta), .25, 1))
	
	u <- runif(I, 0, 1)

	gp = unlist(sapply(u, function(ui) which(ui<(1:K)/K)[1]))

	R1 = rnorm(I, 0, 1) #rpois(I, 3)/sqrt(3)
	R2 = theta[gp] + rnorm(I, 0, 1)

	combinedPvals_allgam <- NULL

	combinedPvals_allgam <- foreach(gammaID = 1:length(gamSeq),
		.combine = 'cbind', .packages = c("senstrat",
			'sensitivitymw','foreach',
			'sensitivitymv','partitions','struncatedP')) %dopar%  {
		

		gam = gamSeq[gammaID]

		pvals = c()

		for(g in 1:K)
			pvals = c(pvals, senmw(cbind(R2[gp==g], R1[gp==g]), gamma=gam, m1=7, m2=8, m=8)$pval)



		
		n = length(pvals)
		#for(u in 1:3){#(n-1)){
		#combinedPvals = 
		#foreach(u=1:1, .combine = 'cbind', 
		#		.packages = c('sensitivitymv','foreach','partitions','struncatedP')) %dopar% {
		u = 1		
				if(u>1){
					temp = rep(0, 16)
				
				} else {
			foo <- function(f, trunc=.15)
						function(x) ifelse(x<.5,f((x-trunc)/(.5-trunc)),0)

			#cat("\nu =",u)
			pvalsord = tail(sort(pvals), n-u+1)
			
			temp1 = c((n-u+1)*min(pvalsord/(1:(n-u+1))),
				truncatedP(pvalsord,trunc=1),
				truncatedP(pvalsord, trunc=.2),
				ifelse(sum(pvalsord) > (n-u+1)*((n-u+1)/(n-u+1+1))^(n-u+1), 1,
						min(sum(pvalsord)^(n-u+1)/factorial(n-u+1),1)))
				
				xi = c('xiIp', 'xiIIp', 'xiIIIp', 'xiIVp',
						'xiVp', 'xiVIp', 'xiIpf', 'xiIIpf', 'xiIIIpf', 'xiIVpf',
						'xiVpf', 'xiVIpf')
						
			#tempstp <- foreach(j=1:12, .combine = 'c',
			#	.packages = c('partitions','struncatedP')) %dopar% {
				for(j in 1:12){
					strcall = paste0('getstruncatedPval1(pvalsord, foo(',
							xi[j],', trunc=.1), N=20, trunc=.1)')
					tempstp = eval(parse(text=strcall))
					temp1  = c(temp1, tempstp)
				}		
				temp1 = c(temp1, truncatedP(pvalsord, trunc=.1),
							poolr::stouffer(pvalsord)$p,
							poolr::tippett(pvalsord)$p
							)
				
			}
			#temp1
			
		#}
		combinedPvals = temp1
		
		#combinedPvals_allgam = cbind(combinedPvals_allgam, combinedPvals)
		combinedPvals
	}
		
	Simes[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[1,], ncol=length(gamSeq))
	Fisher[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[2,], ncol=length(gamSeq))
	tP[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[3,], ncol=length(gamSeq))
		tP2[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[17,], ncol=length(gamSeq))
	SumP[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[4,], ncol=length(gamSeq))
	
	stPI[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[5,], ncol=length(gamSeq))
	stPII[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[6,], ncol=length(gamSeq))
	stPIII[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[7,], ncol=length(gamSeq))
	stPIV[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[8,], ncol=length(gamSeq))
	stPV[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[9,], ncol=length(gamSeq))
	stPVI[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[10,], ncol=length(gamSeq))
	
	stPIf[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[11,], ncol=length(gamSeq))
	stPIIf[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[12,], ncol=length(gamSeq))
	stPIIIf[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[13,], ncol=length(gamSeq))
	stPIVf[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[14,], ncol=length(gamSeq))
	stPVf[itr, 1:1, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[15,], ncol=length(gamSeq))
	stPVIf[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[16,], ncol=length(gamSeq))
		
	stouf[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[18,], ncol=length(gamSeq))
	minP[itr, 1:1, 1:length(gamSeq)] = matrix(combinedPvals_allgam[19,], ncol=length(gamSeq))
	
	if( ( (itr/Itr1)*100 )%%5 == 0){
			val = (Sys.time()-temp)*(Itr1-itr)/(Itr1*.05)
			cat(' -',( (itr/Itr1)*100 ),'%- ', as.numeric(val),attr(val,"units"))
			temp = Sys.time()
			gc()
	}
	
	itr = itr+1
}

stopCluster(cl)

#beep(3)

#apply(Simes<0.05, c(2,3), mean, na.rm=TRUE)
#apply(tP<0.05, c(2, 3), mean, na.rm=TRUE)
#apply(stP<0.05, c(2, 3), mean, na.rm=TRUE)
#apply(SumP<0.05, c(2, 3), mean, na.rm=TRUE)

nmethods = 19
res <- matrix(NA, length(gamSeq), nmethods*K)
rownames(res) <- gamSeq
 
res[,seq(1, nmethods*K, nmethods)] = t(round(apply(Simes<0.05, c(2,3), sum, na.rm=TRUE), 3))

res[,seq(2, nmethods*K, nmethods)] = t(round(apply(SumP<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(3, nmethods*K, nmethods)] = t(round(apply(tP<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(4, nmethods*K, nmethods)] = t(round(apply(stPI<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(5, nmethods*K, nmethods)] = t(round(apply(stPII<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(6, nmethods*K, nmethods)] = t(round(apply(stPIII<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(7, nmethods*K, nmethods)] = t(round(apply(stPIV<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(8, nmethods*K, nmethods)] = t(round(apply(stPV<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(9, nmethods*K, nmethods)] = t(round(apply(stPVI<0.05, c(2,3), sum, na.rm=TRUE), 3))

res[,seq(10, nmethods*K, nmethods)] = t(round(apply(Fisher<0.05, c(2,3), sum, na.rm=TRUE), 3))
 
res[,seq(11, nmethods*K, nmethods)] = t(round(apply(stPIf<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(12, nmethods*K, nmethods)] = t(round(apply(stPIIf<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(13, nmethods*K, nmethods)] = t(round(apply(stPIIIf<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(14, nmethods*K, nmethods)] = t(round(apply(stPIVf<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(15, nmethods*K, nmethods)] = t(round(apply(stPVf<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(16, nmethods*K, nmethods)] = t(round(apply(stPVIf<0.05, c(2,3), sum, na.rm=TRUE), 3))

res[,seq(17, nmethods*K, nmethods)] = t(round(apply(tP2<0.05, c(2,3), sum, na.rm=TRUE), 3))

res[,seq(18, nmethods*K, nmethods)] = t(round(apply(stouf<0.05, c(2,3), sum, na.rm=TRUE), 3))
res[,seq(19, nmethods*K, nmethods)] = t(round(apply(minP<0.05, c(2,3), sum, na.rm=TRUE), 3))
write.table(rbind(Itr1, res[,1:nmethods], NA), col.names =FALSE, file = 'simulation_res.txt', append=TRUE)


