#---------------------------------------------------------------------------------------------
# R code for computing functional diversity in the paper 
# Chiu, C.-H. and Chao, A. (2014). Distance-based functional diversity measures 
# and their decomposition: a framework based on Hill numbers.PLoS ONE 
# 9(7):e100014.
#--------------------------------------------------------------------------------------------

Func2014=function (Dis,abun,q)  
  # input:
  # Dis: species pairwise functional distance matrix.
  # abun: one assemblage species abundance vector data or multi-assemblage
  #species by assemblages matrix data, where row is species and column is 
  # community or site(plot).
  # q: a non-negative value specifying diversity order.
  
  #output:
  # Q: Rao's quadratic entropy for each community.
  # FuncD= functional diversity of each community. 
  # Gamma= functional gamma diversity.
# Alpha= functional alpha diversity.
# Beta= functional beta diversity.
# FunCqN = functional overlap index CqN(similarity index) 
# FunUqN= functional overlap index UqN(similarity index)
{
  if(is.vector(abun)){
    Dis=as.matrix(Dis); 
    n=sum(abun);
    p=abun/n;I=which(p>0);p=p[I];
    Q=c(t(p)%*%Dis[I,I]%*%p);
    temp=p%*%t(p);
    if(q==1){FD= exp(sum(-Dis[I,I]*temp/Q*log(temp/Q)))}
    else{FD=(t(p^q)%*%Dis[I,I]%*%(p/Q)^q)^(1/(1-q));}
    #output=matrix(ncol=2,nrow=1);
    output=c(Q,FD)
    names(output)=c("Q","FuncD")
    return(output);
  }else{
    abun=as.matrix(abun);  
    N=ncol(abun);n=colSums(abun)
    FuncD=numeric(N);Q=numeric(N)
    Dis=as.matrix(Dis); 
    for(i in 1:N){            
      Q[i]=t(abun[,i]/n[i])%*%Dis%*%(abun[,i]/n[i]);
      temp=(abun[,i]/n[i])%*%t(abun[,i]/n[i]);
      I=which(temp>0);
      
      if(q==1){ FuncD[i]= exp(sum(-Dis[I]*temp[I]/Q[i]*log(temp[I]/Q[i]))); }
      else{ FuncD[i]=sum(Dis[I]*(temp[I]/Q[i])^q)^(1/(1-q)); }
    }
    
    gn=sum(abun);
    pop=abun/gn;  
    p=rowSums(pop);gI=which(p>0);p=p[gI];
    gQ=c(t(p)%*%Dis[gI,gI]%*%p);
    
    atemp=0;      
    for(i in 1:N){
      for(j in 1:N){
        pi=pop[,i];pj=pop[,j]                  
        pij=pi%*%t(pj);aI=which(pij>0);
        if(q==1){ atemp=atemp+ sum(-Dis[aI]*pij[aI]/gQ*log(pij[aI]/gQ))}
        else{atemp=atemp+sum(Dis[aI]*(pij[aI]/gQ)^q) }
      }	
    }
    
    if(q==1){
      gFD= exp(sum(-Dis[gI,gI]*(p%*%t(p))/gQ*log( p%*%t(p)/gQ )))
      aFD= exp(atemp)/N^2
      bFD=gFD/aFD
      CqN= 1-log(bFD)/log(N^2);
      UqN=CqN;      
    }else{
      gFD= (t(p^q)%*%Dis[gI,gI]%*%(p/gQ)^q)^(1/(1-q));
      aFD= atemp^(1/(1-q))/N^2
      bFD=gFD/aFD;            
      CqN=1-(bFD^(1-q)-1)/(N^(2-2*q)-1);
      UqN=1-(bFD^(q-1)-1)/(N^(2*q-2)-1);
    }
    output=c(gQ,gFD,aFD,bFD,CqN,UqN)
    names(output)=c("Q","Gamma","Alpha","Beta","FunCqN","FunUqN")
    return( list(Q=Q,FuncD=FuncD,output))               
  }
}


###pairwise functional similarity indices for CqN and UqN similarity indices
#input:
#abun:species by community matrix dataframe, where row is species and
#column is community or site.
# q: order q should be any non-negative value.
#output: pairwise functional similarity matrix 
# CqN: Sorensen-type overlap index (from a local view)(similarity index) 
# UqN: Jaccard-type overlap index (from a regional view)(similarity index)

pairFunc=function(Dis, abun,q){
  N=ncol(abun);
  CqN=matrix(1,ncol=N,nrow=N);UqN=CqN;
  for(i in 1:(N-1)){
    for(j in i:N){
      o=Func2014(Dis,abun[,c(i,j)],q);
      CqN[i,j]=o[[3]][5];CqN[j,i]=CqN[i,j];
      UqN[i,j]=o[[3]][6];UqN[j,i]=UqN[i,j];
    }
  }
  return(list(CqN=CqN,UqN=UqN));
}
