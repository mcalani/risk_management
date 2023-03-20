
replace=1;
binding =1;
planner =0;
load_data;
simSS=simSS2;

cd ..; cd matlab
[state_final,index_b] =  crisis(SET,SEN,simCAY,simSS,state_sim,simB);


simuldata=0;
binding =1;
B0=index_b;
statefinal=state_final;

%% Planner

planner=1; 
binding=1;

load_data


cd ..; cd matlab        
[SimCp,SimPp, SimCAYp]=comparison(B,BORRLIMIT,statefinal,B0,OPTB);
                  
   
    SimCpdev=(SimCp-mean(simC))/mean(simCp)*100;
    SimRERp=(omega^(1/(1+mu))+(1-omega)^(1/(1+mu))*SimPp.^(mu/(1+mu))).^((1+mu)/mu);
    SimRERpdev=(SimRERp-mean(simRERp))/mean(simRERp)*100;
   
%% DE

planner =0;
load_data 

cd ..; cd matlab
[SimC,SimP, SimCAY]=comparison(B,BORRLIMIT,statefinal,B0,OPTB);
  

    SimCdev=(SimC-mean(simC))/mean(simC)*100;
    SimRER=(omega^(1/(1+mu))+(1-omega)^(1/(1+mu))*SimP.^(mu/(1+mu))).^((1+mu)/mu);
    SimRERdev=(SimRER-mean(simRER))/mean(simRER)*100;
    
  
