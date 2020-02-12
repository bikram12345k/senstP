library(senstrat)
library(sensitivitymw)
library(sensitivitymv)
library(struncatedP)
library(partitions)
library(foreach)
library(doSNOW)


cl<-makeCluster(16) 
registerDoSNOW(cl)


I = 1000

## Choose ONE
theta = c(0, 0.05, .05, .15, .3)	# Alternative
theta = rep(0, 5)	# NULL


K = length(theta)
p = c(20, 20, 20, 20, 10)
gamSeq = c(seq(1, 4, .2))

Itr = 14000


Simes = Fisher = tP = SumP = array(NA, c(Itr, K, length(gamSeq)))
stPI = stPII = stPIII = stPIV = stPV = stPVI = array(NA, c(Itr, K, length(gamSeq)))
stPIf = stPIIf = stPIIIf = stPIVf = stPVf = stPVIf = array(NA, c(Itr, K, length(gamSeq)))

temp = Sys.time()
itr=1
Itr1 = Itr
while(itr<= Itr1){
#for(itr in 1:Itr){
		
	gp = sample(1:K, I, replace = TRUE, prob = p)

	
	R1 = rnorm(I, 0, 1/sqrt(2))
	R2 = theta[gp] + rnorm(I, 0, 1/sqrt(2))


	#temp = Sys.time()
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
		combinedPvals = foreach(u=1:3, .combine = 'cbind', 
				.packages = c('sensitivitymv','foreach','partitions','struncatedP')) %dopar% {
				
				if(u>1){
					temp = rep(0, 16)
				
				} else {
			foo <- function(f, trunc=.2)
						function(x) f((x-trunc)/(1-trunc))

			#cat("\nu =",u)
			pvalsord = tail(sort(pvals), n-u+1)
			
			temp = c((n-u+1)*min(pvalsord/(1:(n-u+1))),
				truncatedP(pvalsord,trunc=1),
				truncatedP(pvalsord, trunc=.1),
				ifelse(sum(pvalsord) > (n-u+1)*((n-u+1)/(n-u+1+1))^(n-u+1), 1,
						min(sum(pvalsord)^(n-u+1)/factorial(n-u+1),1)))
				
				xi = c('xiIp', 'xiIIp', 'xiIIIp', 'xiIVp',
						'xiVp', 'xiVIp', 'xiIpf', 'xiIIpf', 'xiIIIpf', 'xiIVpf',
						'xiVpf', 'xiVIpf')
						
			tempstp <- foreach(j=1:12, .combine = 'c',
				.packages = c('partitions','struncatedP')) %dopar% {

					strcall = paste0('getstruncatedPval1(pvalsord, foo(',
							xi[j],', trunc=.1), N=20, trunc=.1)')
					eval(parse(text=strcall))
				}
							
				temp  = c(temp, tempstp)
			}
			temp
			
		}
		
		
		combinedPvals
	}
		
	Simes[itr, 1:3, 1:length(gamSeq)] = matrix(combinedPvals_allgam[1,], ncol=length(gamSeq))
	Fisher[itr, 1:3, 1:length(gamSeq)] = matrix(combinedPvals_allgam[2,], ncol=length(gamSeq))
	tP[itr, 1:3, 1:length(gamSeq)] = matrix(combinedPvals_allgam[3,], ncol=length(gamSeq))
	SumP[itr, 1:3, 1:length(gamSeq)] = matrix(combinedPvals_allgam[4,], ncol=length(gamSeq))
	
	stPI[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[5,], ncol=length(gamSeq))
	stPII[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[6,], ncol=length(gamSeq))
	stPIII[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[7,], ncol=length(gamSeq))
	stPIV[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[8,], ncol=length(gamSeq))
	stPV[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[9,], ncol=length(gamSeq))
	stPVI[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[10,], ncol=length(gamSeq))
	
	stPIf[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[11,], ncol=length(gamSeq))
	stPIIf[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[12,], ncol=length(gamSeq))
	stPIIIf[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[13,], ncol=length(gamSeq))
	stPIVf[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[14,], ncol=length(gamSeq))
	stPVf[itr, 1:3, 1:length(gamSeq)] =  matrix(combinedPvals_allgam[15,], ncol=length(gamSeq))
	stPVIf[itr, 1:3, 1:length(gamSeq)] = matrix(combinedPvals_allgam[16,], ncol=length(gamSeq))
		
	
	if( ( (itr/Itr1)*100 )%%5 == 0){
			val = (Sys.time()-temp)*(Itr1-itr)/(Itr1*.05)
			cat(' -',( (itr/Itr1)*100 ),'%- ', as.numeric(val),attr(val,"units"))
			temp = Sys.time()
			gc()
	}
	
	itr = itr+1
}

#apply(Simes<0.05, c(2,3), mean, na.rm=TRUE)
#apply(tP<0.05, c(2, 3), mean, na.rm=TRUE)
#apply(stP<0.05, c(2, 3), mean, na.rm=TRUE)
#apply(SumP<0.05, c(2, 3), mean, na.rm=TRUE)

nmethods = 16
res <- matrix(NA, length(gamSeq), nmethods*K)
rownames(res) <- gamSeq
 
res[,seq(1, nmethods*K, nmethods)] = t(round(apply(Simes<0.05, c(2,3), mean, na.rm=TRUE), 3))

res[,seq(2, nmethods*K, nmethods)] = t(round(apply(SumP<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(3, nmethods*K, nmethods)] = t(round(apply(tP<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(4, nmethods*K, nmethods)] = t(round(apply(stPI<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(5, nmethods*K, nmethods)] = t(round(apply(stPII<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(6, nmethods*K, nmethods)] = t(round(apply(stPIII<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(7, nmethods*K, nmethods)] = t(round(apply(stPIV<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(8, nmethods*K, nmethods)] = t(round(apply(stPV<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(9, nmethods*K, nmethods)] = t(round(apply(stPVI<0.05, c(2,3), mean, na.rm=TRUE), 3))

res[,seq(10, nmethods*K, nmethods)] = t(round(apply(Fisher<0.05, c(2,3), mean, na.rm=TRUE), 3))
 
res[,seq(11, nmethods*K, nmethods)] = t(round(apply(stPIf<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(12, nmethods*K, nmethods)] = t(round(apply(stPIIf<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(13, nmethods*K, nmethods)] = t(round(apply(stPIIIf<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(14, nmethods*K, nmethods)] = t(round(apply(stPIVf<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(15, nmethods*K, nmethods)] = t(round(apply(stPVf<0.05, c(2,3), mean, na.rm=TRUE), 3))
res[,seq(16, nmethods*K, nmethods)] = t(round(apply(stPVIf<0.05, c(2,3), mean, na.rm=TRUE), 3))

i=1

round(res[,(nmethods*i):(nmethods*(i-1)+1)][, c(setdiff(1:13, 7), 7, 14)]*100) 
