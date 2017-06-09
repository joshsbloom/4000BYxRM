# general Mixed model code for variance component analysis

# y named(!) phenotype vector (n)
# B named list of covariance structures (m x m matrices)
# X a incidence matrix for fixed effects (n x p fixed effects)
# Z a incidence matrix for random effects (n x m) 
# Ze a incidence matrix for random effects for error term only
# reps  (T or F, are there replicates)
# alg derivative based algorithm ('ai' = average information, 'fs' = fisher scoring, or 'nr' = newton-rhapson (fastest to slowest), sqrt  )
# conv.val value for convergence
# Var vector of initialization values for variance components

library(Matrix)
# gpu only works for smaller datasets, otherwise vid card doesn't have enough memory
# advise not using
#library(gputools)

calcMM = function(y, B=NULL,X=NULL, Z=NULL, Ze=NULL, reps=FALSE,
                     alg='ai', conv.val=1e-4, Var=NULL, doGPU=FALSE ,doWood=FALSE){
    strain.names=(names(y))
    unique.sn=unique(strain.names)
    n.to.m=match(strain.names, unique.sn)

    strain.ind  = seq_along(strain.names)
    strain.cnt  = length(unique.sn)
    #for constructing Strain Variance component
    Strain      = Matrix(diag(strain.cnt), sparse=T)

    # If B is NULL then add a covariance term that is the identity matix - will then calculate effect of strain (broad-sense heritability)
    if(is.null(B)) { B=list(Strain=Strain);  } else{
        if (reps) { B=c(B, list(Strain=Strain))  } }    
    # If Z is null and there are no replicates this will make Z a diagonal incidence matrix, otherwise this constructs an incidence matrix based on strain names
    if(is.null(Z)) {   Z=Matrix(0, length(y), strain.cnt,sparse=T);   Z[cbind(strain.ind, n.to.m)]=1 }
    # If X is null assume one fixed effect of population mean
    if(is.null(X) ) {  X=model.matrix(y~1)}
    # If Ze is null assume no error structure
    if(is.null(Ze)) {  Ze=Matrix((diag(length(y))),sparse=T) }

    #number of terms in the structured covariance
    VC.names=paste('sigma', c(names(B), 'E'), sep='')
    N.s.c=length(B)
    Vcmp.cnt=N.s.c+1
    # starting values for VC estimates as 1 / (#of VCs including residual error term)
    if(is.null(Var) ) { Var=rep(1/Vcmp.cnt, Vcmp.cnt) }
    I = matrix(0, ncol= Vcmp.cnt, nrow= Vcmp.cnt)
	s = matrix(0, ncol=1, nrow= Vcmp.cnt)
    
    diffs=rep(10,  Vcmp.cnt)
    # second derivatives of V with respect to the variance components (Lynch and Walsh 27.15)

    VV = list()
    for(i in 1:N.s.c) {
        if(doGPU) {  VV[[i]]=gpuTcrossprod(gpuMatMult(Z, B[[i]]), Z)  }
        else{        VV[[i]]=Z %*% tcrossprod(B[[i]],Z)  }
    }
    VV[[ Vcmp.cnt ]]=Ze 

     solveW = function( AAinv, UU, CCinv) {      AAinv - AAinv %*% UU %*% tcrossprod(solve(CCinv + crossprod(UU, AAinv) %*% UU), UU) %*% AAinv    }
     #new woodbury stuff ----------------3/9/14 ...slow , only good for sparse matrices
     if( doWood==TRUE) { 
         Binvs = list()
         for (i in 1:length(B)){ Binvs[[i]]=solve(B[[i]]) }
         Binvs[[Vcmp.cnt]]=Ze 
    }
     #-------------------3/9/14

    i = 0
    # while the differences haven't converged 
    while ( sum(ifelse(diffs<conv.val, TRUE,FALSE)) <  Vcmp.cnt ) { 
		i = i + 1
        V=matrix(0,length(y), length(y))
	    for( vcs in 1:length(VV)) {  V=V+(VV[[vcs]]*Var[vcs]) }
        
    # new Woodbury stuff-------------------------------
        if(doWood==TRUE) {
            print('Inverting V using woodbury identity')
            Vinv = (1/Var[length(VV)])*Binvs[[length(VV)]]
            for(vcs in ((length(VV)-1):1)  )   {
               if(is.finite(1/Var[vcs])) { 
                 Vinv = solveW(Vinv, Z, (1/Var[vcs])*Binvs[[vcs]])
               } else {Vinv=Vinv}
               print(vcs)
            }               
            print('Done inverting V using woodbury identity')
       }
    #---------------------------------------------------
        if (doWood==FALSE) {
             print('Inverting V')
             Vinv = solve(V)
             print('Done inverting V')
        }
		if(doGPU) {	
            tXVinvX=gpuMatMult(gpuMatMult(t(X),Vinv),X)
            inv.tXVinvX =  solve(tXVinvX) 
            itv         =  gpuMatMult(gpuTcrossprod( inv.tXVinvX, X), Vinv)
            P = Vinv - gpuMatMult(gpuMatMult(Vinv, X), itv) }
        if(doGPU==FALSE){  
            tXVinvX=t(X) %*% Vinv %*% X
            inv.tXVinvX = solve(tXVinvX)
            itv = inv.tXVinvX %*% t(X)%*%Vinv
            P = Vinv - Vinv %*% X %*% itv }

        if(alg=='fs') {print("Fisher scoring algorithm: calculating expected VC Hessian") }
        if(alg=='nr') {print("Netwon rhapson algorithm: calculating observed VC Hessian") }
        if(alg=='ai') {print("Average information algorithm: calculating avg of expected and observed VC Hessians") }
        
        for(ii in 1:Vcmp.cnt) {
           for(jj in ii:Vcmp.cnt) {
                 if (alg=='fs') {    I[ii,jj]= 0.5*sum(diag( ((P%*%VV[[ii]]) %*%P )%*%VV[[jj]])) }
                 if (alg=='nr') {    I[ii,jj]=-0.5*sum(diag(P%*%VV[[ii]]%*%P%*%VV[[jj]])) + (t(y)%*%P%*%VV[[ii]]%*%P%*%VV[[jj]]%*%P%*%y)[1,1] }
                 if (alg=='ai') {    
                     if(doGPU)  {    I[ii,jj]= 0.5*gpuMatMult(gpuMatMult(gpuMatMult(gpuMatMult(gpuMatMult(gpuMatMult(t(y),P),VV[[ii]]),P),VV[[jj]]),P),y) }
                     else{           I[ii,jj]= 0.5*( t(y)%*%P%*%VV[[ii]]%*%P%*%VV[[jj]]%*%P%*%y)[1,1]  } 
                 }
                 print(paste(ii, jj))
                 I[jj,ii]=I[ii,jj]
           }
           if(doGPU) { s[ii,1]= -0.5*sum(diag(gpuMatMult(P,VV[[ii]]))) + .5*(gpuMatMult(gpuMatMult(gpuMatMult(gpuMatMult(t(y),P),VV[[ii]]),P),y) ) }
           else{       s[ii,1]= -0.5*sum(diag(P%*%VV[[ii]])) + .5*(t(y)%*%P%*%VV[[ii]]%*%P%*%y )[1,1] }
        }
        invI = solve(I)
        print(invI)
        print(s) 
        newVar = Var + invI%*%s
        #newVar[newVar<0]=2.6e-9

        for(d in 1:length(diffs)) { diffs[d]=abs(Var[d] - newVar[d]) }
		Var = newVar
        
        cat('\n')
        cat("iteration ", i, '\n')
        cat(VC.names, '\n')
        cat(Var, '\n')
        
        Bhat= itv %*% y
        cat("Fixed Effects, Bhat = ", as.matrix(Bhat), '\n')
        det.tXVinvX=determinant(tXVinvX, logarithm=TRUE)
        det.tXVinvX=det.tXVinvX$modulus*det.tXVinvX$sign
        det.V =determinant(V, logarithm=TRUE)
        det.V=det.V$modulus*det.V$sign

        LL = -.5 * (det.tXVinvX + det.V + log(t(y) %*% P %*% y) )
        cat("Log Likelihood = " , as.matrix(LL), '\n')
        cat("VC convergence vals", '\n')
        cat(diffs, '\n')
	}
    	cat('\n')
	return(list(Var=Var, invI=invI, W=Vinv, Bhat=Bhat, llik=LL))
}



