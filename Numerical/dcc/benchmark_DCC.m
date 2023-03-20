clear;
close all;
clc;

%% Read data
filePath = matlab.desktop.editor.getActiveFilename; %fprintf('%s\n',filePath)
filefolderpath = filePath(1:(end-16));
cd(filefolderpath)

Rstock = readtable('Rstockjl.csv', 'ReadVariableNames', false, 'HeaderLines', 1);
CovQ   = readtable('CovQjl.csv', 'ReadVariableNames', false, 'HeaderLines', 1);

%% Set variables
y = table2array(Rstock(:,2:end));

[t , k] =size(y);


iv = table2array(CovQ(:,2:end));
iv = reshape(iv,k,k,t);

%in-sample days
t1 = 1000; 
y = y(1:t1,:);
iv = iv(:,:,1:t1);

%plot data
ii=3;
plot(y(:,ii))
hold on
plot(reshape(iv(ii,ii,:),1,t1))

%% Step 1: Estimate GARCH and get residuals, D = sqrt(diag(H))

params = zeros(3,k);
Ht = zeros(t1,k);
for i =1:k
    resid       = y(:,i)-mean(y(:,i));
    params(:,i) = tarch(resid,1,0,1);
    [params(:,i),~,Ht(:,i),~,~,~,~] = tarch(resid,1,0,1);
end

ii=2;
plot(y(:,ii))
hold on
plot(sqrt(Ht(:,ii)))

Dt = sqrt(Ht); %sigma

err = y./Dt; % GARCH residuals
plot(err)


%% Step 2: Use residuals to estimate DCC and get Qt and Rt objects
options = optimset('TolFun',1e-12,'MaxIter',2000,'Display','iter',...
                    'MaxFunEvals',8000);
                
[par_DCC, ~ ,Ht, ~, ~, ~] = dcc(err,[],1,0,1,[],[],[],[],'2-stage',[],[],options);



%% BEKK, RARCH, HEAVY, RCC, OGARCH, GOGARCH?





5









