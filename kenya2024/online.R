winSum<-function(x,step){
    cs<-c(0,cumsum(x))
    cs[-c(1:step)]-cs[-(length(cs) - 1:step+1)]
}

latticeplot<-function(mat,col=rev(heat.colors(50))[-1],at=0:10/10,...)
  lattice::levelplot(mat,scales=list(x=list(rot=45)),col.regions=col,colorkey=list(at=at,label=list(at=at,lab=at)),at=at,...)
 
info<-function(){
cat("xpd=T\n")
cat("par(mar=c(5.1,4.1,4.1,2.1))\n")
}
pdff<-function(...)
pdf("~/public/albrecht/temp.pdf",...)

norm<-function(x)
  x/sum(x)

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

brewer<-function(x="Dark2"){
    
pdf("/dev/null")
  library(RColorBrewer)
  palette(brewer.pal(8, x))
dev.off()

}



color.palette <- function(steps, n.steps.between=NULL, ...){

    if(is.null(n.steps.between)) n.steps.between <- rep(0, (length(steps)-1))
    if(length(n.steps.between) != length(steps)-1) stop("Must have one less n.steps.between value than steps")

    fill.steps <- cumsum(rep(1, length(steps))+c(0,n.steps.between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[,fill.steps] <- col2rgb(steps)

    for(i in which(n.steps.between>0)){
        col.start=RGB[,fill.steps[i]]
        col.end=RGB[,fill.steps[i+1]]
        for(j in seq(3)){
            vals <- seq(col.start[j], col.end[j], length.out=n.steps.between[i]+2)[2:(2+n.steps.between[i]-1)]
            RGB[j,(fill.steps[i]+1):(fill.steps[i+1]-1)] <- vals
        }
    }

    new.steps <- rgb(RGB[1,], RGB[2,], RGB[3,], maxColorValue = 255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

 

myPals<-function(x=0){

    if(x==0)
        pal <- color.palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","red4","darkred","black"), space="rgb")
    else{
        
        require(wesanderson)
        
        aa<-c( "GrandBudapest","Moonrise1", "Royal1", "Moonrise2", "Cavalcanti", "Royal2","GrandBudapest2", "Moonrise3", "Chevalier" , "BottleRocket" ,"darjeeling", "darjeeling2")
        pal <- wes.palette(name = aa[x], type = "continuous")
    }
    
    pal
}


colorAnd<-function(x="colorblind6b"){

    if(x=="line")
        palette(c("#b2182b", "#2166ac", "#4DAF4A", "#FF7F00", "#F781BF","#984EA3"))
    else if(x=="rasmus")
        palette(c("mistyrose","lavender","lightyellow","lightblue","lightgreen","seashell","lightcyan"))
    else if(x=="colorblind6b")
        palette(c("#b2182b","#ef8a62","#fddbc7","#d1e5f0","#67a9cf","#2166ac"))
    else if(x=="anders")
        palette(c("darkgreen","#00A600FF","yellow","#E9BD3AFF","orange","coral4","red4","black"))
    else if(x=="large")
       palette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"))
 else if(x=="ffs")
     palette(c("#56B4E9", "#E69F00", "#009E73",  "#F0E442", "#0072B2", "#D55E00", "#CC79A7",'#444444'))
    
 }

  
logThin <- function(x,thinTo=10000,fracNoThin=0.01){
    N <-length(x)
    nKeep <- ceiling(fracNoThin*N)
    keep<-c(1:nKeep,round(exp( 1:thinTo*(log(N-nKeep)/thinTo) ))+nKeep)
    keep
}

qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
legend("topleft",paste("lambda=",lambda))
}




 qqp <-
function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,thinLarge=TRUE,col=1,...){
    
    if(length(col)>1)
        col<-col[!is.na(x)]
    x<-x[!is.na(x)]
    if(!missing(maxLogP))
        x[x<10^-maxLogP]<-10^-maxLogP
    N<-length(x)
    ord<-order(x)
    x<-x[ord]
    if(length(col)>1)
        col<-col[ord]
     
    if(thinLarge & N>1.1e5){
        #n1 <- round(N/100)
        #keep <- c(1:n1,1:n1*100)

        ##keep top 10000 without thinning
        ## thin to 10000 points
        keep <- logThin(x,thinTo=10000, fracNoThin=1000/N)
        e<- -log((1:N-0.5)[keep]/N,10)
        x <- x[keep]
    }
    else
        e<- -log((1:N-0.5)/N,10)
 
  if(add)
    points(e,-log(x,10),col=col,...)
  else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,col=col,...)
    abline(0,1,col=2,lwd=2)
  }
#legend("topleft",paste("lambda=",lambda))

    if(ci){
        #https://www.tandfonline.com/doi/abs/10.1080/00949658008810388https://www.tandfonline.com/doi/abs/10.1080/00949658008810388
        #concentration bands
         if(thinLarge & N>1.1e5){
             c97.5<-qbeta(0.975,keep,N-(keep)+1)
             c02.5<-qbeta(0.025,keep,N-(keep)+1)
         }
         else{
            c97.5<-qbeta(0.975,1:N,N-(1:N)+1)
            c02.5<-qbeta(0.025,1:N,N-(1:N)+1)
   
        }
        lines(e,-log(c97.5,10))
        lines(e,-log(c02.5,10))
    }
}


qqpOld<-function(x,ci=TRUE,add=FALSE,ylab="Observed log10(p-value)",xlab="Expected log10(p-value)",maxLogP,col=1,...){
    
    if(length(col)>1)
        col<-col[!is.na(x)] # NB ida changed this from col<-col[!is.a(x)]
    x<-x[!is.na(x)]
    if(!missing(maxLogP))
        x[x<10^-maxLogP]<-10^-maxLogP
    N<-length(x)
    ord<-order(x)
    x<-x[ord]
    if(length(col)>1)
        col<-col[ord]
#lambda<-round(median(x)/qchisq(0.5,1),2)
  e<- -log((1:N-0.5)/N,10)
  if(add)
    points(e,-log(x,10),col=col,...)
  else{
    plot(e,-log(x,10),ylab=ylab,xlab=xlab,col=col,...)
    abline(0,1,col=2,lwd=2)
  }
#legend("topleft",paste("lambda=",lambda))

    if(ci){
        #https://www.tandfonline.com/doi/abs/10.1080/00949658008810388https://www.tandfonline.com/doi/abs/10.1080/00949658008810388
        #concentration bands
    c97.5<-qbeta(0.975,1:N,N-(1:N)+1)
    c02.5<-qbeta(0.025,1:N,N-(1:N)+1) 
    lines(e,-log(c97.5,10))
    lines(e,-log(c02.5,10))
  }
}
 #qqp(runif(1000))


bases<-c("A","C","G","T")
GENO<-c("AA","AC","AG","AT","CC","CG","CT","GG","GT","TT")
GENO2<-c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
GENOrev<-c("TT","GT","CT","AT","GG","CG","AG","CC","AC","AA")
GENO2rev<-rev(GENO2)

wiki<-function(x,file,...){
  
  if(class(x)=="data.frame")
    x<-as.matrix(x)
  if(class(x)=="matrix")
    wiki.matrix(x,file,...)
  


}

wiki.matrix<-function(x,file,class="wikitable",border="1",...){
write(paste("{| class=\"",class,"\" border=\"",border,"\"",sep=""),file=file)
temp<-NULL

if(!is.null(rownames(x)))
  temp<-" "
if(!is.null(colnames(x)))
  writeWikiRow(colnames(x),file=file,type="!",row=temp)

for(tal in 1:nrow(x))
  writeWikiRow(x[tal,],file=file,row=rownames(x)[tal])


write("|}",file=file,append=TRUE)
}

writeWikiRow<-function(x,file,type="|",row=NULL){
  write("|-",file=file,append=TRUE)
  if(!is.null(row))
  write(paste("!",row),file=file,append=TRUE)
  for(tal in 1:length(x))
  write(paste(type,x[tal]),file=file,append=TRUE)
}
wikii<-function(x,...)
wiki(x,file="~/public/albrecht/temp.txt",...)


#################################################
swichMe<-function(w,x){
  n<-length(w)
  if(n==nrow(x)){
    x<-t(x)
    warning("transposing x (slower)")
  }
  if(n!=ncol(x)){
    stop("wrong Dim")
  } 
  x[w+(0:(n-1))*nrow(x)]
}

qtrans<-function(x){
    k<-!is.na(x)
    ran<-rank(x[k])
    y<-qnorm((1:sum(k)-0.5)/sum(k))
    x[k]<-y[ran]
    x
}

##################################
plinkOld<-function(plinkFile){
    pl<-snpMatrix::read.plink(plinkFile)
    ind<-rownames(pl)
    snp<-colnames(pl)
    geno<-as.integer(as.integer(pl)-1)
    dim(geno)<-dim(pl)
    geno[geno==-1]<-NA
    rownames(geno)<-ind
    colnames(geno)<-colnames(pl)
    bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
    fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
    list(geno=geno,bim=bim,fam=fam,pl=pl)
}
plink<-function(plinkFile){
    pl <- snpStats::read.plink(plinkFile)
    ind<-pl$fam[,1]
    snp<-pl$map[,2]
    geno<-as.integer(as.integer(pl$genotypes)-1)
    dim(geno)<-c(length(ind),length(snp))
    geno[geno==-1]<-NA
    rownames(geno)<-ind
    colnames(geno)<-colnames(pl)
    bim<-read.table(paste0(plinkFile,".bim"),as.is=T,header=F)
    fam<-read.table(paste0(plinkFile,".fam"),as.is=T,header=F)
    list(geno=geno,bim=bim,fam=fam,pl=pl)
}

###########################3
writePlink<-function(x,info,file,IID,FID,SEX,PHE){
  n<-ncol(x)
  if(missing(IID))
    IID<-1:n
  if(missing(FID))
    FID<-IID
  if(missing(SEX))
    SEX<-rep(0,n)
  if(missing(PHE))
    PHE<-rep(-9,n)
  
  for(g in GENO2)
    x[x==g]<-paste(strsplit(g,"")[[1]],collapse=" ")

  if(any(is.na(info)))
    stop("no missing in info allowed")
  write.table(cbind(info,x),file=file,col=F,row=F,qu=F,na="0 0")
  write.table(cbind(IID,FID,0,0,SEX,PHE),file=sub("tped","tfam",file),col=F,row=F,qu=F)
  
}


##########################################
sameStrand<-function(A1,B1,A2,B2,wild=0){
    a1<-A1==wild | A1 == A2 | A1==B2
    b1<-B1==wild | B1 == A2 | B1==B2

    a2<-A2==wild | A2 == A1 | A2==B1
    b2<-B2==wild | B2 == A1 | B2==B1
    (a1&b1) | (a2&b2)
}


####
fed<-function(x,n=6){
    if(is.matrix(x) | is.data.frame(x))
        return(x[1:min(nrow(x),n),1:min(ncol(x),n)])
    else
        return(head(x,n))
}
    

##################
#make lover case
lower<-function(x){
    c(letters,letters)[match(x,c(LETTERS,letters))]
}


##################
flip<-function(x){
   c(rev(bases),rev(lower(bases)))[match(x,c(bases,lower(bases)))]
}

###
mplot<-function(x,m,col,ylim,type="l",...){
    if(is.data.frame(m))
        m<-as.matrix(m)
    if(missing(col)){
        col<-1:nrow(m)
    }
    if(missing(ylim))
        ylim<-range(m,na.rm=T)
    plot(x,m[1,],col=1,ylim=ylim,type=type,...)
    for(cc in 2:nrow(m))
        lines(x,m[cc,],type=type,col=cc,...)
}
#mplot(m,0:10/100)

##
fold <- function(x){
if(length(x)%%2==0)
    stop("equal number of cats\n")
hal<-(length(x)+1)/2
x[1:hal]+c(rev(x)[1:(hal-1)],0)
}


legend2<-function(x, y = NULL, legend, fill = NULL, col = par("col"), 
    border = "black", lty, lwd, pch, angle = 45, density = NULL, 
    bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"), 
    box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, 
    xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 
        0.5), text.width = NULL, text.col = par("col"), text.font = NULL, 
    merge = do.lines && has.pch, trace = FALSE, plot = TRUE, 
    ncol = 1, horiz = FALSE, title = NULL, inset = 0, xpd, title.col = text.col, 
    title.adj = 0.5, seg.len = 2){

    myList<-as.list( sys.call() )
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0:1, 0:1, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    myList$xpd=TRUE
    do.call("legend",myList)
}
legendAnd<-function(){

    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0:1, 0:1, type = "n", bty = "n", xaxt = "n", yaxt = "n")
   # myList$xpd=TRUE
   # do.call("legend",myList)
}


xget<-function(x,n=1){
sapply(x,function(y) y[n])

}

ieatIn<-function(a,b){
#k<-paste(bim[,1],bim[,4])%in%paste(gl[,1],gl[,2])
    fun<-function(x){
        aa<-a[x,2]
        bb<-b[b[,1]==a[x[1],1],2]
        aa%in%bb
    }
    ta <- tapply(1:nrow(a),factor(a[,1],levels=unique(a[,1])),fun)
    unlist(ta)
}





## change a matrix order or row/colums without changing the lower/upper triangle

changeMat <- function(matOrg,rank){
    #matOrg <- matrix(1:16,4);rank<-c(4,1,3,2)

    ## mat <- matOrg[length(rank):1,]
    mmat <- matOrg
    newMat <- mmat
    newMat[!is.na(mmat)] <- NA
    rownames(newMat) <- rownames(matOrg)[rank]
    colnames(newMat) <- rownames(newMat)
    dimnames(mmat) <- NULL
    m <- reshape2::melt(mmat)
    m$n1 <- order(rank)[m[,1]]
    m$n2 <- order(rank)[m[,2]]

    m2 <- m[order(m$n1),]
    m2 <- m2[order(m2$n2),]
    
    
    for(i in 1:nrow(m)){
        x<-mmat[m2[i,1],m2[i,2]]
        if( (m[i,1]<m[i,2]) != (m2[i,1]<m2[i,2]))
            x <- mmat[m2[i,2],m2[i,1]]
        
        newMat[m[i,1],m[i,2]] <- x
      #  fed(newMat)

    }
    newMat   
}


### change R with in terminal
wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
   
  options(width=as.integer(howWide))
}

#wideScreen()


##x is a N*M matrix
##y is a M vector
## each column is divided by a entry in M
matScale <- function(x,y,type=1){
    if(type==1)
        x %*% diag(1/y)
    else if(type==2)
        x/y[col(x)]

}






manPlot <- function( x, chr, thinLarge=TRUE,thinTo=2e4,pass,cap=1e-30, collar = c( "darkblue", "#67a9cf", "orange", "grey" ), main = "" ,bonfline=TRUE){
     keep <- !is.na( x )
     x <- x[keep]
     chr <- chr[keep]
     N<-length(x)
     bonf <- -log( 0.05/N, base = 10 )
      k<-1:N
      if(thinLarge & N>1.1e5){
          
          keepRank <- logThin(x,thinTo=100000,fracNoThin=10000/N)
          k <- sort(order(x)[keepRank])
          
#          num <- 10
#         if(N>1e6)
#             num <- 25
#         if(N>5e6)
#             num<- 50
#        n1 <- round(N/num)
#        k <- c(1:n1,1:n1*num)
#        k <- sort(order(x)[k])

          x <- x[k]

        chr <- chr[k]
         if ( !missing( pass ) )
            pass <- pass[k]
      }
    
  
     x <- ifelse( x<cap, cap, x ) 
     col <- chr%%2 + 1
     pch <- rep( 16, length( x ) )
     if ( !missing( pass ) ){
         col[!pass[keep]] <- 4
         pch[!pass[keep]] <- 1
     }
#     par( mar = c( 5.1, 5.1, 4.1, 1.1 ) )
     maxy = max( (-log( x, base = 10 ) ) )
     plot((1:N)[k],-log( x, base = 10 ), col = collar[col], ylab = expression(-log[10]( italic( P ) ) ), xlab = "Chromosomes", main = main, cex = 1, lwd = 2, pch = pch, axes = F, cex.lab = 1, cex.main = 2, ylim = c( 0, maxy + 0.05*maxy ) )
     box( )
     axis( 2, las = 1, cex.axis = 1.8 )
     if(bonfline)
      abline( h = bonf, lty = 2, col = 1 )
     ##legend( 0, t-0.5, "Bonferroni correction", lty = 2, bty = "n", col = 1, cex = 2 )
     ##legend( "topright", c( "known", "novel" ), pch = c( 4, 1 ), bty = "n", cex = 2 )


        if(maxy>-log10(cap)*.99){
            text(0,-log10(cap)*1.03,"capped",adj=0)
            abline(h=-log10(cap),lty=2,col="grey")
        }

     
 }



qqPlot <- function( x, cap=1e-30,main = "" ){
        keep <- !is.na( x )
        x <- x[keep]
        x <- ifelse( x<cap, cap, x )
        maxy = max( (-log( x, base = 10 ) ) )
           
        qqp( x, pch = 16, col = "darkblue", main = main, las = 1, cex.lab = 1, cex.main = 2, ylim = c( 0, maxy + 0.05*maxy ),
            xlab = expression( Expected~~-log[10]( italic( P ) ) ), ylab = expression( Observed~~-log[10]( italic( P ) ) ), cex.axis = 1.8 )
       
        if(maxy>-log10(cap)*.99){
            text(0,-log10(cap)*1.03,"capped",adj=0)
            abline(h=-log10(cap),lty=2,col="grey")
        }

}
    
