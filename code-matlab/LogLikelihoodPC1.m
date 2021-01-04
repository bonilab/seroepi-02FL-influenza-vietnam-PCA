function llh = LogLikelihoodPC1(prms)

    global PC1 Age;
    assert(length(PC1) == length(Age));
    
    B = prms(1);
    C = prms(2);
    standarDeviation = prms(3);
    lambda = prms(4);
    llh = 0;

    for i=1:length(PC1)
        age = Age(i) - 0.5;
        obs = PC1(i);
    
       Expected = B - C*exp(-lambda*age);     
       Residual = obs - Expected;
       templlh = - log(standarDeviation) - Residual^2/(2*standarDeviation^2);
       
       llh = llh + templlh;

    end

    llh = -1*llh;

end