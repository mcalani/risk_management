
load data_arg.txt

year=data_arg(:,1);
pn=data_arg(:,2);  %  real exchange rate (log-diff)
cay=data_arg(:,3);  % current account-GDP ratio
tby=data_arg(:,4);   % tradabe balance-GDP ratio
gdp=data_arg(:,5);  % GDP in units of tradables
cons=data_arg(:,6);  % Consumption in real terms

omega=   0.3020;
mu =   0.2050;
rer=(  omega^(1/(1+mu))  +  (1-omega)^(1/(1+mu)).*pn.^(mu/(1+mu))).^((1+mu)/mu);
rerdev=[nan ;log(rer(2:end))-log(rer(1:end-1))];
clear mu omega

N=length(year); % time frame is 1965-2007
smooth=100; %  hpfilter smoothing parameter

loggdp=log(gdp);
logcons=log(cons);
hp_gdp=hpfilter(loggdp,smooth);
hp_cons=hpfilter(logcons,smooth);
gdpdev=loggdp-hp_gdp;
consdev=logcons-hp_cons;

stddev=nanstd([consdev,rerdev, cay, tby]);

corrgdp=zeros(4,1);
num=corrcoef([gdpdev, consdev],'rows','complete');
corrgdp(1)=num(1,2);
num=corrcoef([gdpdev, rerdev],'rows','complete');
corrgdp(2)=num(1,2);
num=corrcoef([gdpdev, cay],'rows','complete');
corrgdp(3)=num(1,2);
num=corrcoef([gdpdev, tby],'rows','complete');
corrgdp(4)=num(1,2);


