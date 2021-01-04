function y = besttiter( x )

    fourFoldFailed = 0;
    %y = NaN;
    y = -99;
    
    if x(3) == 1 % meaning that the four-fold is valid
        
        %if isnan( x(1) )
        if x(1) < -50
            fourFoldFailed = 1;
        end
        
        %if ~isnan( x(1) )
        if x(1) >= 0
            y = x(1);
        end
        
    end
    
    %if x(3) == 0 || isnan( x(3) ) || fourFoldFailed
    if x(3) == 0 || x(3)<-50  || fourFoldFailed
        
        %if x(4) == 1 && ~isnan(x(2))
        if x(4) == 1 && x(2) >= 0
            y = x(2);
        end
        
    end


