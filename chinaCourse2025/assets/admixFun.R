
getRMSD<-function(Qest,Qtrue){
    require(gtools)

  npop<-nrow(Qtrue)
  perms<-permutations(npop,npop)
  Nperm<-dim(perms)[1]

  theMin<-c()
  for(g in 1:nrow(Qest)){
    cat("N perm",Nperm,"\n")
    RMSEs<-rep(0,Nperm)
    Qnew<-Qest[-g,]
    for(i in 1:Nperm){
      RMSEs[i]<-sqrt(sum((Qtrue-Qnew[perms[i,],])^2))
    }
    theMin[g]<-min(RMSEs)
  }
  
  w<-which.min(theMin)
  Qnew<-Qest[-w,]

  RMSEs<-rep(0,Nperm)
  for(i in 1:Nperm){
    RMSEs[i]<-sqrt(sum((Qtrue-Qnew[perms[i,],])^2))
  }
 #print(RMSEs)
  Qnew<-rbind(Qnew[perms[which.min(RMSEs),],],Qest[w,])
  Qnew
}

## returns the Q matrix with minimum difference to the old Q(that has 1 less population)
getFast<-function(Q,Qold){
    npop<-nrow(Qold)
    res<-c()
    for(g in 1:nrow(Q)){
        w<-rowSums((rep(Q[g,],each=npop)-Qold)^2)
        res<-rbind(res,c(which.min(w),min(w)))
    }
    dub <- duplicated(res[,1])
    dd<-res[dub,1]
    ww<-which.max(res[res[,1]==dd,2])
    res[which(res[,1]==dd)[ww],1]<-npop+1
    Q[order(res[,1]),]
}


colorFun<-function(x,Q){
  o<-Q[,x]
  if(length(x)==1)
     o<-matrix(o,ncol=1)

  K<-nrow(o)
  res<-matrix(0,ncol=ncol(o),nrow=K^2)
  most<-order(rowSums(o),decreasing=T)
  o<-o[,order(o[most[1],],decreasing=T)]
  if(length(x)==1)
    o<-matrix(o,ncol=1)

  for(k in 1:K)
  res[ (k-1)*K+most[k],]<-o[most[k],]
  res
}
####


mkOrd<-function(p,res,popOrd){
  x<-res[,popOrd==p]
  if(is.matrix(x)){
  w<-which.max(rowSums(x))
  m<-(1:length(popOrd))[popOrd==p][1]-1
  order(x[w,])+m
  }else{
  m<-(1:length(popOrd))[popOrd==p][1]-1
  c(1)+m
  }
}







plotMulti<-function(allQ,Kall,pop,reorder=1,fast=FALSE,lwd=2,lty=2){

      n<-length(pop)

     if(reorder==1){
        res<-allQ[[length(allQ)]]
        ord<-order(pop)
        res <-res[,ord]    
        popOrd<-pop[ord]
        u <-unlist(lapply(unique(pop), mkOrd,res=res,popOrd=popOrd))
        ordd <- ord[u]
    }else
        ordd <- order(pop)

    pop<-pop[ordd]
    Qold<-NA

    par(mfrow=c(length(Kall)+1,1))
cat("K=")
    for(K in Kall){
        par(mar=c(.1,5.1,.6,2.1))
        
        cat(K," ")
        Q<-allQ[[K]]
        Q<-Q[,ordd]
        if(K!=Kall[1]){
            if(fast)
                Q<-getFast(Q,Qold)
            else
                Q<-getRMSD(Q,Qold)

        }
        Qold<-Q
        
        Q<-do.call(cbind,tapply(1:length(pop),match(pop,unique(pop)),colorFun,Q=Q))
        
        

        ta<-tapply(pop,pop,length)
        small<- names(ta)[ta==1]
        h<- barplot(Q,border=NA,col=1:K,space=0,ylab="Admixture 
proportion",main=paste("K = ",K,sep=""))
        abline(v=tapply(h,pop,max)+0.5,col="black",lwd=lwd,lty=lty)

        med<-tapply(h,pop,median)
        
        
        
    }
    par(mar=c(2.9,5.1,.1,2.1))
h<- barplot(Q,border=NA,col="transparent",space=0, axes=F)
    text(med,rep(0.5,length(unique(pop))),names(med),xpd=TRUE,srt = 45,cex=1.2,adj=1)
    
    
}

