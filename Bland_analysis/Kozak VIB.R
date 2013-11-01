library(RODBC)
temp=odbcConnectAccess("C:/Users/Frank/Desktop/Demeritt WP- Bland.mdb")
trees=sqlFetch(temp,"Trees")
#plots=sqlFetch(temp,"Plots")
close(temp)
names(trees)[names(trees)=="TotHt"] = "HT"

########################### Kozak Model 2 ######################################
# Parameter estimated taken from Li et al 2012
Kozak_02CRVIB=function(DHT,HT,DBHO,HTLB)  
  
{              z=DHT/HT;          p=1.3/HT; Xi=(1-z^(1/3))/(1-p^(1/3)); Qi=1-z^(1/3);CR=(HT-HTLB)/HT
               
               a0=1.049;a1=1.008;a2=-0.046;b1=0.381;b2=-0.86;b3=0.344;b4=4.608;b5=0.112;b6=-0.552;b7=0
               
               d=(a0*(DBHO^a1)*(HT^a2))*Xi^(b1*z^4+b2*(exp(- DBHO/HT))
                                            
                                            +b3*Xi^0.1+b4*(1/DBHO)+b5*HT^Qi+b6*Xi+b7*CR)
               
               return (d)}

# Volume prediction

smalians<-function(r1,r2,len){
  L=(r1/2)^2*pi
  S=(r2/2)^2*pi
  vol=((L+S)/2)*len
  return(round(vol,4))}

LW.KVolCRVIB=function(DBHO,HT,HTLB){
  
  sgmts <- 100; L  <- HT / sgmts;i <- 0; volme <- 0; aa=0
  
  while(i<(sgmts-1)){
    H1  <- (L * i)
    H2 <- (L * (i+1))
    dib1 <- Kozak_02CRVIB(DHT=H1,HT,DBHO,HTLB)
    dib2 <- Kozak_02CRVIB(DHT=H2,HT,DBHO,HTLB)
    volme <-volme+smalians(dib1,dib2,L*100)
    dib = (dib1+dib2)/2
    
    aa = rbind(aa,dib)
    
    i <- i+1}
  
  return(round(volme/1E6,6))}

livetrees=subset(trees,Status==1,select=,)

livetrees$KozakVIB=mapply(LW.KVolCRVIB,DBHO=livetrees$DBH,HT=livetrees$HT,HTLB=livetrees$LLB)
names(livetrees)[names(livetrees)=="HT"]="TotHt"

write.csv(livetrees, file = ("C:/Users/Frank/Desktop/R output/KozakVIB.csv"), row.names = FALSE)

