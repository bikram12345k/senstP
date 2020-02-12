library(senstrat)
library(sensitivitymw)
library(sensitivitymv)
library(struncatedP)
library(partitions)
library(foreach)
library(doSNOW)

lead <- read.csv('lead_absorption1.csv')



gamSeq = seq(1, 5, .05)

res_original = matrix(NA, length(gamSeq), 5)

for(gamma in gamSeq){

	pvals <- c()

	## Exposed child vs control child
	pvals <- c(pvals, senmw(cbind(lead$lead[1:33], lead$lead[(33+1):(33+33)]), 
								gamma=gamma, m1=7, m2=8, m=8)$pval)

	## High and Medium vs Low exposure
	who <- lead$exposed_gp==1
	st = rep(1, sum(who))
	z = lead[who, 'job_exposure'] %in% c("h","m")
	sc = rscores(lead$lead[who], z, st, method="U788")
	pvals <- c(pvals, sen2sample(sc, z, alternative="greater", gamma=gamma )$pval)

	## High vs Medium exposure
	who <- lead$exposed_gp==1 & lead$job_exposure != "l"
	z = lead[who, 'job_exposure'] == "h"
	st = rep(1, sum(who))
	sc = rscores(lead$lead[who], z, st, method="U788")
	pvals <- c(pvals, sen2sample(sc, z, alternative="greater", gamma=gamma )$pval)

	## Good vs Moderate or poor hygine
	who <- lead$exposed_gp==1 & lead$job_exposure == "h"
	st = rep(1, sum(who))
	z = lead[who, 'hygiene'] != "g"
	sc = rscores(lead$lead[who], z, st, method="U788")
	pvals <- c(pvals, sen2sample(sc, z, alternative="greater", gamma=gamma )$pval)

	## Moderate vs poor hygine
	who <- lead$exposed_gp==1 & lead$job_exposure == "h" & lead$hygiene != 'g'
	st = rep(1, sum(who))
	z = lead[who, 'hygiene'] == "p"
	sc = rscores(lead$lead[who], z, st, method="U788")
	pvals <- c(pvals, sen2sample(sc, z, alternative="greater", gamma=gamma )$pval)

	#pvals
	res_original[which(gamma==gamSeq),] = pvals

	n = length(pvals)

}

#round(res, 6)
rownames(res_original) = gamSeq
#round(res_original, 6)

## Table 4
tab = round(t(apply(res_original[,c(1,2,3,4,5)], 1, function(p) c(
			(5)*min(tail(sort(p), 5)/(1:(5))),
			sort(p)[1]*(5),
			pchisq(-2*log(prod(tail(sort(p), 5))), 2*(5), lower.tail=FALSE),
			truncatedP(p, trunc=.2), 
			getstruncatedPval1(p, xi=foo(xiVIpf, trunc=.2), trunc=.2, N=20))
			)), 5)
rownames(tab) = gamSeq
colnames(tab) = c('Simes', 'Holm', 'Fisher', 'truncatedP', 'struncatedP')
tab


##########################################################################
### Supplementary functions
rscores <- function (y, z, st = NULL, tau = 0, method=c("wilcoxon", "cs1", "cs2", "savage", "cs4", "sh", "U545", "U788")) {
    method = match.arg(method)
    stopifnot(length(tau) == 1)
    stopifnot(length(y) == length(z))
    if (is.null(st)) 
        st <- rep(1, length(y))
    stopifnot(length(st) == length(y))
    ust <- unique(st)
    nst <- length(ust)
    if (tau != 0) 
        y <- y - z * tau
	
	sc <- rep(NA, length(y))

    for (i in 1:nst) {
        who <- st == ust[i]
        yi <- y[who]
        ni <- length(yi)
	  if (ni == 1) {
            sc[who] <- 0
        } else if (ni >= 2) {
		sc[who] <- score(rank(yi), ni, method=method)
		
        }
    }

   sc
}

score <- function(j, Ns, method=c("W", "wilcoxon", "cs1", "cs2", "savage", "cs4", "sh", "U545", "U788"), a=5){
	
	method = match.arg(method)

	sj <- switch(method, W = j, wilcoxon = j/(Ns+1),
				cs1 = sapply(j, function(j) prod( seq(j/(Ns+1), (j+a-2)/(Ns+1), by=1/(Ns+1)))),
				cs2 = (j/(Ns+1))^(a-1),
				savage = sapply(j, function(j) sum(1/(1-1:j/(Ns+1)))),
				cs4 = -log(1-j/(Ns+1)), 
				sh = sapply(j, function(j) ifelse(j<a, 0, prod(seq((j-a+1)/(Ns+1),(j-1)/Ns+1), by=1/(Ns+1) )) ),
				U545 = sapply(j, function(j) ifelse(j<5, 0, (Ns*choose(Ns, 5)^(-1))*sum( choose((j-1), (4-1):(5-1))*choose((Ns-j), (5-4):(5-5)) ))),
				U788 = sapply(j, function(j) ifelse(j<8, 0, (Ns*choose(Ns, 8)^(-1))*sum( choose((j-1), (7-1):(8-1))*choose((Ns-j), (8-7):(8-8)) )))				)
	sj
}


