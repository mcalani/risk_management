  function [cdev,pdev,caydev]= comparison(B,borrlimit,stateimp,B0,OPTB);

global kappat kappan rbar gamma sigma  rbar omega  SET SEN mu 


TSIM=3;

SIMB=zeros(1,TSIM+1);
SIMP=zeros(1,TSIM);
SIMY=zeros(1,TSIM);

SIMConstraint=zeros(1,TSIM);
SIMB(1)=B0;
SIMyt=zeros(1,TSIM);
SIMyn=zeros(1,TSIM);
SIMBU=zeros(1,TSIM);


for t=1:TSIM

    SIMBU(t+1)=interp1(B,OPTB(:,stateimp(t)),SIMB(t));
    SIMLimit(t)=interp1(B,borrlimit(:,stateimp(t)),SIMB(t));

    
    if SIMBU(t+1)<=SIMLimit(t)     
    SIMConstraint(t)=1;
    SIMB(t+1)=SIMLimit(t);
    else      
    SIMB(t+1)= SIMBU(t+1)  ;
    SIMConstraint(t)=0;
    end
    
    SIMB(t+1)=max(B(1),SIMB(t+1));
    SIMB(t+1)=min(SIMB(t+1),B(end));
    
    SIMCT(t)= SET(stateimp(t)) + SIMB(t)*(1+rbar)  -SIMB(t+1)   ;
    SIMyt(t)=SET(stateimp(t));
    SIMyn(t)=SEN(stateimp(t));
    SIMP(t)=(1-omega)/omega*(SIMCT(t)/SIMyn(t))^(1+mu);  
     SIMY(t)= SIMyt(t)  +  SIMP(t)*SIMyn(t); 
     SIMCAY(t)=(SIMB(t+1)-SIMB(t))/SIMY(t)*100;
end

  SIMC =(omega*(SIMCT.^(-mu))+(1-omega)*(SIMyn.^(-mu))).^(-1/mu);

cdev=SIMC(TSIM);
caydev=SIMCAY(TSIM);
pdev=SIMP(TSIM);

