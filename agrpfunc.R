##usage£º
##input:d1-expression matrix (Genes*Samples) of cancer class£¬d2-expression matrix of normal tissues
##output: a list containing GRP statistics and its p-values for each gene
agrpfunc<-function(d1,d2){
if(class(d1)=="numeric") d1=matrix(d1,nrow=1);
if(class(d2)=="numeric") d2=matrix(d2,nrow=1);
up1=apply(d1,2,function(x){rowSums(x>=d2)/ncol(d2)});
up2=apply(d2,2,function(x){rowSums(x<=d1)/ncol(d1)});
pu=(ncol(d1)*rowMeans(up1)+ncol(d2)*rowMeans(up2))/(ncol(d1)+ncol(d2))
down1=apply(d1,2,function(x){rowSums(x<=d2)/ncol(d2)});
down2=apply(d2,2,function(x){rowSums(x>d1)/ncol(d1)});
pd=(ncol(d1)*rowMeans(down1)+ncol(d2)*rowMeans(down2))/(ncol(d1)+ncol(d2));
pud=pu-pd;
n=ncol(d1);
m=ncol(d2);
p.value <- 2 * (1 - pnorm(abs(pud), 0,
                                 ((n^2+m^2)/(2*n*m*(n+m))/0.75)^(1/2)));
  GRP <- list( pud = pud, p.value = p.value);
  return(GRP);
}