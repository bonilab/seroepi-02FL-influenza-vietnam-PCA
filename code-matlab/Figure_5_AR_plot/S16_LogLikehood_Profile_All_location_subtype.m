clc;
clear all;
%==============================================
global PC1 Age

load ../var_for_attackrate_fitting.mat
locnames = {'HC', 'KH', 'HU', 'DL', 'KG', 'ST', 'AG', 'DT', 'QN', 'BD' };

%subtype = 'H1N1';
subtype = 'H3N2';

%===========Search Routine Settings =================
options = optimset('Display','iter','TolFun',1e-8,'MaxIter', 10000, 'MaxFunEvals',10000);
iter = 50;

if strcmp(subtype,'H1N1')
    D  = SCOREH1( LOC==1 | LOC==2 | LOC==3, : );
else
    assert(strcmp(subtype,'H3N2'))
    D  = SCOREH3(  LOC==1 | LOC==2 | LOC==3, : );
end
 
PC1 = D(:,1);
Age = AGE( LOC==1 | LOC==2 | LOC==3 );
 
record = [];
mle = dlmread(['./' subtype '/AR_Estimate_All_loc.txt']);
record = [record; mle];
 
LB      = [0 0 0 0];
UB      = [30 30 30 5];

for i=1:13
    ini = mle(1:4);
    ini(4) = mle(4) - i*0.005;
    UB(4) = ini(4);
    LB(4) = ini(4);
    [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
    record = [record; [x fval]];
end
    
for i=1:13
    ini = mle(1:4);
    ini(4) = mle(4) + i*0.005;
    UB(4) = ini(4);
    LB(4) = ini(4);
    [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
    record = [record; [x fval]];
end
    
dlmwrite(['./' subtype '/LLHProfile_all_locations.txt'],record, 'delimiter','\t','precision', 6)
    

for lc=1:3
    if strcmp(subtype,'H1N1')
        D  = SCOREH1( LOC==lc, : );
    else
        assert(strcmp(subtype,'H3N2'))
        D  = SCOREH3( LOC==lc, : );
    end
    PC1 = D(:,1);
    Age = AGE( LOC==lc );
    
   
    record = [];
    mle = dlmread(['./' subtype '/AR_Estimate_' locnames{lc} '.txt']);
    record = [record; mle];
    
    for i=1:13
        ini = mle(1:4);
        ini(4) = mle(4) - i*0.005;
        UB(4) = ini(4);
        LB(4) = ini(4);
        [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
        record = [record; [x fval]];
    end
    
    for i=1:13
        ini = mle(1:4);
        ini(4) = mle(4) + i*0.005;
        UB(4) = ini(4);
        LB(4) = ini(4);
        [x,fval,exitflag,output]=fminsearchbnd(@LogLikelihoodPC1,ini,LB,UB,options);
        record = [record; [x fval]];
    end
    
    dlmwrite(['./' subtype '/LLHProfile_' locnames{lc} '.txt'],record, 'delimiter','\t','precision', 6)
end
