function [est, lowerbound, upperbound] = Find95CI(prms,llh)
   [prms_sorted, order] = sort(prms);
   llh_sorted = llh(order);
   
   mle =  min(llh_sorted);
   
  % ============find estimate================== 
   est = prms_sorted(find(llh_sorted == mle));
   est = est(1)
  
  % ============find lowerbound================ 
    idx = find(llh_sorted < mle + 1.92);
    idx = min(idx);
    x1 = prms_sorted(idx-1);
    x2 = prms_sorted(idx);
    llh1 = llh_sorted(idx-1);
    llh2 = llh_sorted(idx);
    targetllh = mle + 1.92;
    
    assert(llh1 > mle + 1.92);
    assert(llh2 < mle + 1.92);

    lowerbound = x1 + (x2-x1)*(targetllh-llh1)/(llh2-llh1);
    
   %==========find upperbound==================
   
    idx = find(llh_sorted < mle + 1.92);
    idx = max(idx);
    x1 = prms_sorted(idx);
    x2 = prms_sorted(idx+1);
    llh1 = llh_sorted(idx);
    llh2 = llh_sorted(idx+1);
    targetllh = mle + 1.92;
    
    assert(llh1 < mle + 1.92);
    assert(llh2 > mle + 1.92);

    upperbound = x1 + (x2-x1)*(targetllh-llh1)/(llh2-llh1);
   
   

end