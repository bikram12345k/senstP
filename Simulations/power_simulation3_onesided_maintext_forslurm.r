args = commandArgs(trailingOnly=TRUE)
split=as.numeric(args)[1]

### Vary sample size
library(senstrat)
library(sensitivitymw)
library(sensitivitymv)
library(struncatedP)
library(partitions)
library(foreach)
library(doSNOW)
library(BSDA)

#library(beepr)

set.seed(50*split)

cl<-makeCluster(6) 
registerDoSNOW(cl)

resall = NULL
resall_1 = NULL

## Choose ONE
#theta = c(-.5, -.25, -.1, 0.1, .15, .25)	# Alternative
theta = c(-1, -.5, 0.1, 0.15, .2, .3)	# Alternative
#theta = rep(0, 6)	# NULL


K = length(theta)

for(I in seq(10,100,5)*K){
	cat('\n',I,'\n')
	

	p = c(20, 20, 20, 20, 10)
	gamSeq = c(seq(1, 1.25, .25))	# For this simulation we do not use gamma, but keep it to be of length 2

	Itr = 65


	Simes = Fisher = tP = tP2 = SumP = array(NA, c(Itr, K, length(gamSeq)))
	stPI = stPII = stPIII = stPIV = stPV = stPVI = array(NA, c(Itr, K, length(gamSeq)))
	stPIf = stPIIf = stPIIIf = stPIVf = stPVf = stPVIf = array(NA, c(Itr, K, length(gamSeq)))

	temp = Sys.time()
	itr=1
	Itr1 = Itr
	while(itr<= Itr1){
	#for(itr in 1:Itr){
			
		# Create a stratified fisher's randomization structure
		pvals = NULL
		for(gp in 1:K){
			N = I/K

			x = rnorm(N, m=theta[gp], 1)
			
			#pvals = c(pvals, t.test(x, alternative='greater')$p.value)	
			pvals = c(pvals, z.test(x, sigma.x=1, alternative='greater')$p.value)
		}
		
		pvals = pmin(1, pvals)

		combinedPvals_allgam <- NULL

			
		n = length(pvals)
		
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
		
			for(j in 1:12){
				strcall = paste0('getstruncatedPval1(pvalsord, foo(',
						xi[j],', trunc=.1), N=20, trunc=.1)')
				tempstp = eval(parse(text=strcall))
				temp1  = c(temp1, tempstp)
			}		
			temp1 = c(temp1, truncatedP(pvalsord, trunc=.1))
		
		}
			
		combinedPvals = temp1
			
		combinedPvals_allgam = cbind(combinedPvals, combinedPvals)
			
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
			
		
		if( ( (itr/Itr1)*100 )%%5 == 0){
				val = (Sys.time()-temp)*(Itr1-itr)/(Itr1*.05)
				cat(' -',( (itr/Itr1)*100 ),'%- ', as.numeric(val),attr(val,"units"))
				temp = Sys.time()
				gc()
		}
		
		itr = itr+1
	}


	#beep(3)

	#apply(Simes<0.05, c(2,3), mean, na.rm=TRUE)
	#apply(tP<0.05, c(2, 3), mean, na.rm=TRUE)
	#apply(stP<0.05, c(2, 3), mean, na.rm=TRUE)
	#apply(SumP<0.05, c(2, 3), mean, na.rm=TRUE)

	nmethods = 17
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


	resall = rbind(resall, c(I, Itr1, res[1,1:nmethods]))
	
	nmethods = 17
	res <- matrix(NA, length(gamSeq), nmethods*K)
	rownames(res) <- gamSeq
	 
	res[,seq(1, nmethods*K, nmethods)] = t(round(apply(Simes<0.1, c(2,3), sum, na.rm=TRUE), 3))

	res[,seq(2, nmethods*K, nmethods)] = t(round(apply(SumP<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(3, nmethods*K, nmethods)] = t(round(apply(tP<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(4, nmethods*K, nmethods)] = t(round(apply(stPI<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(5, nmethods*K, nmethods)] = t(round(apply(stPII<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(6, nmethods*K, nmethods)] = t(round(apply(stPIII<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(7, nmethods*K, nmethods)] = t(round(apply(stPIV<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(8, nmethods*K, nmethods)] = t(round(apply(stPV<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(9, nmethods*K, nmethods)] = t(round(apply(stPVI<0.1, c(2,3), sum, na.rm=TRUE), 3))

	res[,seq(10, nmethods*K, nmethods)] = t(round(apply(Fisher<0.1, c(2,3), sum, na.rm=TRUE), 3))
	 
	res[,seq(11, nmethods*K, nmethods)] = t(round(apply(stPIf<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(12, nmethods*K, nmethods)] = t(round(apply(stPIIf<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(13, nmethods*K, nmethods)] = t(round(apply(stPIIIf<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(14, nmethods*K, nmethods)] = t(round(apply(stPIVf<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(15, nmethods*K, nmethods)] = t(round(apply(stPVf<0.1, c(2,3), sum, na.rm=TRUE), 3))
	res[,seq(16, nmethods*K, nmethods)] = t(round(apply(stPVIf<0.1, c(2,3), sum, na.rm=TRUE), 3))

	res[,seq(17, nmethods*K, nmethods)] = t(round(apply(tP2<0.1, c(2,3), sum, na.rm=TRUE), 3))


	resall_1 = rbind(resall_1, c(I, Itr1, res[1,1:nmethods]))
}

write.table(resall, col.names =FALSE, row.names=FALSE,
	file = 'simulation_res_onesided.txt', append=TRUE)

write.table(resall_1, col.names =FALSE, row.names=FALSE,
	file = 'simulation_res_onesided_1.txt', append=TRUE)
	
stopCluster(cl)

