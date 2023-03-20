%% LOADING DATA

NB =800;
NIT=4;
NIN=NIT;
NSS=NIT*NIN;

global SIG BETA mu omega rbar kappat kappan SET SEN


 cd ..; cd fortran; cd  output
 load gridd.txt
 load parameters.txt
 B=gridd;
SIG=parameters(1); BETA=parameters(2) ;mu=parameters(3); omega=parameters(4) ; rbar=parameters(5) ;kappat=parameters(6); kappan=parameters(7);


cd equil
load shockprocess.txt
Prob=shockprocess(1:NSS*NSS); Prob=reshape(Prob',NSS,NSS);
SET=shockprocess(NSS*NSS+1:NSS*NSS+NSS);
SEN=shockprocess(NSS*NSS+NSS+1:end);


if planner==0  
cd ..;cd EQUIL
elseif planner==1
cd ..; cd planner
 load ktight.txt
 load tax.txt
 tax=reshape(tax,NB,NSS);
ktightp=reshape(ktight,NB,NSS);
end


load guesspn.txt
load guessp.txt
load guessmu.txt
load guessc.txt
load valuef.txt
load constraint.txt
load blimit.txt

load probsb.txt 
load suddens.txt
load simul.txt

  
cd ..; 
load state_simTT.txt; state_sim=state_simTT;
cd ..

%% VARIABLES


value=valuef;

OPTB=reshape(guessp,NB,NSS);
CONSTRAINT=reshape(constraint,NB,NSS);
BORRLIMIT=reshape(blimit,NB,NSS);
prob=reshape(probsb,NB,NSS);
valuef=reshape(valuef,NB,NSS);
probsb=reshape(probsb,NB,NSS);


simB=simul(:,1);
simY=simul(:,2);
simCA=simul(:,3);
simCAY=simul(:,4);
simC=simul(:,5);
simRER=simul(:,6);


simSS1=suddens(:,1);
simSS2=suddens(:,2); 


simTB=simCA-simB*rbar;
simTBY=simTB./simY;


h=length(simY);
[f,xi] = ksdensity(simB);
densityb=f/h;


if planner ==1 
probsbp=probsb;  
OPTBp=OPTB;

BORRLIMITp=BORRLIMIT;
CONSTRAINTp=CONSTRAINT(1:end,:);
valuefp=reshape(valuef,NB,NSS);
taxp=tax;

simBp=simB;
simCp=simC;
simCAYp=simCAY;
simYp=simY;
simRERp=simRER;

simTBYp=simTBY;

simSS1p=simSS1;
simSS2p=simSS2;
densitybp=densityb;
xip=xi;
end
 
