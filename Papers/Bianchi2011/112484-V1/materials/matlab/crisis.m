function [state_final,index_b]= crisis(SET,SEN,simCAY,simSS,state_sim,simB)
                        
cevent= find(simSS>0); %identifies periods with SS
cevent=cevent(1:end-2);
ssevent=cevent;

if (length(cevent)==0)
    disp('no crisis')
    pause
end

CAY3=simCAY(ssevent-0);

pn3=median(CAY3);

Z=[CAY3 cevent ]; 
[x,ind]=sort(Z(:,1),'ascend');
Z=Z(ind,:);
index= find(Z(:,1)>pn3,1);  

inde=index;
index= Z(inde,2) ;
Z(inde,1);


index_b=simB(index-2);

state_final1=state_sim(index-2);
state_final2=state_sim(index-1);
state_final3=state_sim(index-0);

state_final=[state_final1 state_final2 state_final3];

