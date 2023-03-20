%Lecture Notes 7: Dynamic Programming
%Value Function Iteration
tic
sigma= 1.5; %utility parameter
delta=0.1; %depreciation rate
beta=0.95; %discount factor
alpha=0.30; %capital elasticity of output

nbk=1000; %number of data points in the grid
crit=1; %convergence criterion
epsi=1e-6; %convergence parater

ks= ((1-beta*(1-delta))/(alpha*beta))^(1/(alpha-1));

dev=0.9; %max devition from ss
kmin=(1-dev)*ks; % lower bound on the grid
kmax=(1+dev)*ks; % upper bound on the grid

dk=(kmax-kmin)/(nbk-1); %implied increment
kgrid = linspace(kmin, kmax, nbk)'; %builds the grid
v= zeros(nbk,1); %value function
dr = zeros(nbk,1); % decision rule 
tv=zeros(nbk,1);

while crit>epsi
    for i=1:nbk
        tmp=(kgrid(i)^alpha+(1-delta)*kgrid(i)-kmin);
        imax=min(floor(tmp/dk)+1,nbk);
        c=kgrid(i)^alpha+(1-delta)*kgrid(i)-kgrid(1:imax);
        util= (c.^(1-sigma)-1)/(1-sigma);
        [tv(i),dr(i)]=max(util+beta*v(1:imax));
    end;
    crit=max(abs(tv-v)); %compute convergence criterion
    v=tv;%update de value function
end
toc

%final solution
kp=kgrid(dr);
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        