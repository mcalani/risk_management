%% Table 2

tex=str2mat( 'Consumption','Current Account-GDP','Real exchange rate');

table('Table 2: Severity of Crises (calculated for median crisis in current simulation)',strvcat('VARIABLE','DE','SP'  ),str2mat(tex),....
  [        SimCdev             SimCpdev;....
          SimCAY             SimCAYp  ;...
           SimRERdev       SimRERpdev ]  ,12,8,1);

%% Table 3

tex=str2mat( 'Consumption','RER','cay','tby');

 table('Standard deviation',strvcat('VARIABLE','DE','SP','DATA'),str2mat(tex),....
[                std(simC)/mean(simC)*100        std(simCp)/mean(simCp)*100          stddev(1)*100 ;...
                 std(simRER)/mean(simRER)*100   std(simRERp)/mean(simRERp)*100     stddev(2)*100    ;...       
                std(simCAY)*100                  std(simCAYp)*100                  stddev(3);
                std(simTBY)*100                std(simTBYp)*100                  stddev(4);  ],10,8,1)



tex=str2mat( 'Consumption','RER','cay','tby');

 table('Correlation with output in units of tradables',strvcat('VARIABLE','DE','SP','DATA'),str2mat(tex),....
 [                    rhoyc rhoycp  corrgdp(1)        ;
                rhoyrer rhoyrerp  corrgdp(2)  ;...
                rhoycay rhoycayp corrgdp(3)  ;
                 rhoytby rhoytbyp  corrgdp(4)    ],10,8,2)

             