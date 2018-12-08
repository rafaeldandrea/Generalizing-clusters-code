## Indicator variable, = 1 if running on clusters, 0 if running on personal computer.
flux=as.numeric(1*(Sys.info()[1]!='Windows'))
if(flux) setwd('~/Metrics_Project/Runs/')
if(flux) ind=as.numeric(commandArgs(TRUE))[1] else ind=15

library(plyr); library(untb)

## Saves dtf as a file under the name filename, then attaches header to top and saves again.
WriteHeader=function(dtf,filename,header=NULL,append=FALSE){
	if(!append | !file.exists(filename)){
		write.table(dtf,filename,quote=FALSE,row.names=FALSE)
		header=paste0('## Created on ',date(),'\n',header,'\n')
		txt=c(header,readLines(filename))
		write(txt,filename)
	}else write.table(dtf,filename,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)
}


## ----------------------------------------------- List of models/scenarios ----------------------------------------------------------
## Parameters identical to baseline (ML41) unless otherwise noted.
## ML41  -> LV model, rho = 4, r0 = 1, w = .063 (13 niches), m = .08, theta = 50 (baseline)
## ML45  -> LV model, r0 = x(1-x) (0 at the edges, 1 in the middle)
## ML21  -> LV model, rho = 2
## FD    -> LV model, "fitness differences" (peaked r0 and alpha = 1)
## NTRL  -> LV model, completely neutral (r0 = 1 and alpha = 1)
## RM0   -> Rosenzweig-MacArthur no resource depletion
## RM1   -> Rosenzweig-MacArhtur with resource depletion
## SA0   -> Schwilk-Ackerly no dispersal limitation
## SA1   -> Schwilk-Ackerly with dispersal limitation
## HCCM  -> Competition-colonization tradeoff
## IM1-4 -> LV model, same as ML41 with m = .005, .01, .05, .15
## NC1-4 -> LV model, same as ML41 with theta = 5, 10, 20, 30 (species-to-niches ratio ~ 3, 8, 12, 20. Compare with 30 for theta = 50)
## -----------------------------------------------------------------------------------------------------------------------------------



models=c('ML41','ML45','RM0','RM1','SA0','SA1','HCCM','FD','NTRL','ML21',paste0('IM',1:4),paste0('NC',1:4))
scenarios=merge(data.frame(scen=models),data.frame(run=1:100))
scen=scenarios$scen[ind]
run=scenarios$run[ind]


## ------------------------------------------------------ Regional pool ------------------------------------------------------------
if(!scen%in%paste0('NC',1:4)) theta=50 else theta=c(5,10,20,30)[match(scen,paste0('NC',1:4))]
set.seed(run); Nmet=sample(as.numeric(rand.neutral(150e3,theta))); S=length(Nmet)
set.seed(run); traitmet=sort(runif(S))
##----------------------------------------------------------------------------------------------------------------------------------
	
## ----------------------- Parameters --------------------------
## Community size and number of steps
if(!scen%in%c('SA0','SA1')){ J=21e3; numsteps=1e8}

## Immigration rate
if(!scen%in%paste0('IM',1:4)) m=0.08 else m=c(.005,.01,.05,.15)[match(scen,paste0('IM',1:4))]

## Intrinsic growth rates
r0met=rep(1,S)
if(scen%in%c('ML45','FD')) r0met=seq(0,1,l=S)*(1-seq(0,1,l=S))

## Competition kernel
if(scen%in%c('FD','NTRL')) alphamet=matrix(1,S,S)
if(scen=='ML21'){ 
	d=as.matrix(dist(traitmet))
	alphamet=exp(-(d/.1)^2)
	alphamet=alphamet/max(alphamet)
}
if(scen%in%c('ML41','ML45',paste0('IM',1:4),paste0('NC',1:4))){
	w=.063
	d=as.matrix(dist(traitmet))
	alphamet=exp(-(d/w)^4)
	alphamet=alphamet/max(alphamet)
}
if(scen%in%c('RM0','RM1')){
	nR=500; sigma=.03 
	if(scen=='RM0') K=10 else if(scen=='RM1') K=400
	res.traits=seq(0,1,l=nR)
	restraits=resX=sort(runif(nR)); X=X0=rep(5,nR)	
	alphamet=outer(restraits,traitmet,function(s,r) exp(-((s-r)/sigma)^2))
}
if(scen=='HCCM'){
	f=exp(seq(traitmet[1],log(500*traitmet[S]),l=S)) 	## linear f distribution leads to less coexistence when deterministic, cf. Deterministic Scenarios 2.txt.
														## linear f also doesn't help in the stochastic form. 
														## See Note below.
	G0=outer(f,f,function(x,y) 1/2*(1-tanh(.15*(y-x))))
	alphamet=f*G0+t(t(G0)*f)
	r0met=f-min(f)
	
	## Note: I am defining f as a scaled exponential of the species "trait". In other words, we are imagining the trait of the species as the log 
	## of its fecundity modulo a scaling constant.
}

## Dispersal rates
if(scen=='HCCM') fecmig=f else fecmig=1

## -------------------------------------------------------------

	
## ---------------- Initial Community - All Scenarios except Schwilk-Ackerly model --------------------
if(!scen%in%c('SA0','SA1')){
	com=sample(traitmet,size=J,replace=TRUE,prob=Nmet)
	trait=plyr::count(com)$x
	N=plyr::count(com)$freq
	if(scen%in%c('RM0','RM1')){
		alpha=alphamet[,traitmet%in%trait]; g=alpha*X
	}else{ 
		alpha=alphamet[traitmet%in%trait,traitmet%in%trait]; B=r0met[traitmet%in%trait]; D=alpha%*%N
	}
}
##------------------------------------------------------------------------------------------------------

## --------------------------------------------- Dynamics of All Scenarios except Schwilk-Ackerly model ------------------------------------
set.seed(run)
if(!scen%in%c('SA0','SA1','RM0','RM1')){
	time=0; while(time<numsteps){
		if(flux==0 & time%%1e5==0) plot(trait,N,t='h',main=paste0(scen,'\n',time/1e6,'m steps'),las=1)
		time=time+1
		death=sample(seq_along(trait),size=1,prob=D*N)
		N[death]=N[death]-1
		D=D-as.double(alpha[,death])
		if(N[death]==0){
			trait=trait[-death]
			B=B[-death]
			D=D[-death]
			alpha=as.matrix(alpha[-death,-death])
			N=N[-death]
		}
		if(runif(1)>m){
			birth=sample(seq_along(trait),size=1,prob=B*N)
			N[birth]=N[birth]+1
			D=D+as.double(alpha[,birth])
		}else{
			migrant=sample(traitmet,size=1,prob=fecmig*Nmet)
			if(migrant%in%trait){ 
				N[trait==migrant]=N[trait==migrant]+1
				D=D+as.double(alpha[,trait==migrant])
			}else{
				trait=c(trait,migrant)
				N=c(N,1)
				B=c(B,r0met[traitmet==migrant])
				B=B[order(trait)]
				N=N[order(trait)]
				trait=sort(trait)
				alpha=alphamet[traitmet%in%trait,traitmet%in%trait]
				D=alpha%*%N
			} 
		}
	}
}
if(scen%in%c('RM0','RM1')){
	t1=Sys.time()
	if(flux==0) par(mfrow=c(2,1),mar=c(2,4,2,2))
	time=0; while(time<numsteps){		
		time=time+1
		if(flux==0 & time%%1e4==0){ 
			estimatedtime=50e6/time*(Sys.time()-t1)/60
			Xvec=rep(0,length(resX)); Xvec[resX%in%restraits]=X
			plot(resX,Xvec,t='h',las=1,xlim=c(0,1),main=paste('Resources - Est time to finish (hrs)',round(estimatedtime)),ylim=c(0,max(X)))
			plot(trait,N,t='h',las=1,xlim=c(0,1),main=paste(scen,' , ',time/1e6,'million deaths'),ylim=c(0,max(N)))
			legend('top',legend=paste0(round(100*length(X)/length(resX)),'% surviving prey'),bty='n')
		}
		
		## ----------------- Resource Update ------------------
		death=sample(length(restraits),size=1,prob=g%*%N)
		X[death]=X[death]-1
		if(X[death]==0){
			restraits=restraits[-death]
			alpha=alpha[-death,]
			X=X[-death]
		}			
		
		birth=sample(length(restraits),size=1,prob=pmax(0,X*(1-X/K)))
		X[birth]=X[birth]+1		

		g=alpha*X
		## -------------------------------------------------------
		
		## ------------ Consumer Update --------------------------
		death=sample(seq(length(trait)),size=1,prob=N)
		N[death]=N[death]-1
		if(N[death]==0){
			trait=trait[-death]
			alpha=alpha[,-death]
			N=N[-death]
		}
		g=alpha*X
		
		if(runif(1)>m){
			birth=sample(seq(length(trait)),size=1,prob=colSums(g)*N)
			N[birth]=N[birth]+1		
		}else{
			migrant=sample(traitmet,size=1,prob=Nmet)
			if(migrant%in%trait){ 
				N[trait==migrant]=N[trait==migrant]+1
			}else{
				trait=c(trait,migrant)
				N=c(N,1)
				N=N[order(trait)]
				trait=sort(trait)
				alpha=outer(restraits,trait,function(s,r) exp(-((s-r)/sigma)^2))
			} 	
		}
		g=alpha*X
	}
	res=data.frame(trait=resX,X0=X0,X=rep(0,length(X0))); res$X[resX%in%restraits]=X
}
## -----------------------------------------------------------------------------------------------------------------------------------------

## ---------------------------- Save data - All scenarios except Schwilk-Ackerly model -----------------------------------------------------
if(!scen%in%c('SA0','SA1')){
	dat=cbind(trait=traitmet,N=rep(0,S),meta=Nmet)
	dat[dat[,1]%in%trait,2]=N
	dat=data.frame(dat)
		
	if(flux){	
		WriteHeader(dat,paste0(scen,'_run_',run,'.txt'),header=paste0('## m = ',m))
		if(scen%in%c('RM0','RM1')) WriteHeader(res,paste0(scen,'_run_',run,'_resources.txt'),header=paste0('## m = ',m))
	}
}
## ----------------------------------------------------------------------------------------------------------------------------------------------

## ------------------------------------------------------------ Schwilk-Ackerly model --------------------------------------------------------- 
if(scen%in%c('SA0','SA1')){ 		## Schwilk-Ackerly model		
	J=1000							# community size = number of patches. we need to keep it small because the simulation is spatially explicit
	Erange=1; Ex=seq(J)/J*Erange	# environmental value of each patch
	epss=.05						# width of species tolerance to env conditions
	Es=traitmet*Erange				# trait value of each species in the pool
	
	F=outer(seq(S),seq(J),function(s,x) exp(-.5*((Ex[x]-Es[s])/epss)^2))				# Fitness of each species to each patch.
	sigma=ifelse(scen=='SA1',50,1e5); D=exp(-.5*(as.matrix(dist(seq(J)))/sigma)^2) 		# Dispersal ability (flat if scen==SA0)

	## ----------------------------------------- Initial community -----------------------------------------------------
	xsid=sample(S,size=J,replace=TRUE,prob=Nmet)	# In this model the community is represented by this vector, 
													# the species identity of the individual occupying each site x.
	## -----------------------------------------------------------------------------------------------------------------

	simtime=0
	N=sapply(seq(S),function(s) sum(xsid==s))
	while(simtime<=1e6){	#numsteps can be smaller because of reduced community size
		simtime=simtime+1
		x=sample(J,1)
		if(scen=='SA0') localseeds=sapply(seq(S),function(s) sum(xsid[-x]==s)) else localseeds=sapply(seq(S),function(s) sum(D[which(xsid[-x]==s),x]))
		if(flux==0 & simtime%%1e3==0) plot(Es,N,t='h',las=1,xlab='Trait',ylab='Abundance',main=paste0(scen,'\n',simtime/1e6,'m steps'))
		im=plyr::count(sample(S,size=round(m/(1-m)*sum(localseeds)),replace=TRUE,prob=Nmet)); localseeds[im$x]=localseeds[im$x]+im$freq	# *See Note below
		R=localseeds*F[,x]
		xsid[x]=sample(S,1,prob=R)
		N=sapply(seq(S),function(s) sum(xsid==s))
	}
	
	if(flux){ 	
		dat=data.frame(trait=Es,N=N,meta=Nmet)
		WriteHeader(dat,paste0(scen,'_run_',run,'.txt'),header=paste0('## m = ',m))
		WriteHeader(data.frame(site=Ex,species=xsid,trait=Es[xsid]),paste0(scen,'_run_',run,'_raw.txt'),header=paste0('## m = ',m))
	}
	
	## Note on immigration: in this model immigration is implemented differently than in the LV models as here the state variable
	## is the vector of individuals not the vector of species. In order to ensure that a fracion m of recruitments are from outside,
	## we add m/(1-m)*J*n, where n is the number of seeds produced by each individual in the community.
}
## ------------------------------------------------------------------------------------------------------------------------------------------
