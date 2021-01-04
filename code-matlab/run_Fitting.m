clc;
clear all;
%==============================================
global PC1 Age

strBase = 'AR25.Sigma0p70.rep01.k3';
strInfile = ['./AgePC1_Simulated/AgePC1.' strBase '.csv'];


data = csvread( strInfile , 0 , 0 );
PC1 = data(:,2);
Age = data(:,1);
%Age = Age - 0.5;
%===========Search Routine Settings =================
options = optimset('Display','iter','TolFun',1e-8,'MaxIter', 10000, 'MaxFunEvals',10000);
iter = 50;

LB      = [0 0 0 0];
UB      = [30 30 30 5];
% initial = random('unif',LB,UB)

%====================================================

 recordPrms = [];
 record_llh = 10000000000;
    
 
 for it = 1:iter
     ini = random('unif',LB,UB);
     while isnan(LogLikelihoodPC1(ini))
        ini = random('unif',LB,UB);
     end
        
     [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
     if fval < record_llh
        record_llh = fval;
        recordPrms = x;
     end
end

dlmwrite('AttackRate.txt',[recordPrms record_llh],'delimiter', '\t','precision', 6);
 




