clc;
clear all;
%==============================================
global PC1 Age

load ../var_for_attackrate_fitting.mat
locnames = {'HC', 'KH', 'HU', 'DL', 'KG', 'ST', 'AG', 'DT', 'QN', 'BD' };

subtype = 'H1N1';
%subtype = 'H3N2';

%===========Search Routine Settings =================
options = optimset('Display','iter','TolFun',1e-8,'MaxIter', 10000, 'MaxFunEvals',10000);
iter = 50;

LB      = [0 0 0 0];
UB      = [30 30 30 5];
% initial = random('unif',LB,UB)

%====================================================

 if strcmp(subtype,'H1N1')
    D  = SCOREH1( LOC==1 | LOC==2 | LOC==3, : );
 else
    assert(strcmp(subtype,'H3N2'))
    D  = SCOREH3(  LOC==1 | LOC==2 | LOC==3, : );
 end
 
 PC1 = D(:,1);
 Age = AGE( LOC==1 | LOC==2 | LOC==3 );
 
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

dlmwrite(['./' subtype '/AR_Estimate_All_loc.txt'],[recordPrms record_llh],'delimiter', '\t','precision', 6) 
 
 



%     
%     
% for lc=1:3
%    
%     
%     
%     recordPrms = [];
%     record_llh = 10000000000;
%     
%     for it = 1:iter
%         ini = random('unif',LB,UB);
%         while isnan(LogLikelihoodPC1(ini))
%             ini = random('unif',LB,UB);
%         end
%         
%         [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
%         if fval < record_llh
%             record_llh = fval;
%             recordPrms = x;
%         end
%             
%     end
%     dlmwrite(['./' subtype '/AR_Estimate_All_loc.txt'],[recordPrms record_llh],'delimiter', '\t','precision', 6)
% end
% 
% 
%     
% 
