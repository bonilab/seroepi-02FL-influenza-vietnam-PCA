% start_w is the start position of the left-most window
%
% end_w is the end position of the right-most window
%
% ww is the window width
%
% xx_sorted is the x values sorted in order
%
% yy is the corresponding y values
%
% yy_smooth is the smoothed y values
%
% qtile .... put .975 here if you want a 95% prediction interval

function AA = get_sigmas_in_windows(  start_w, end_w, ww , xx_sorted, yy, yy_smooth , qtile )

        % number of windows 
        numw = (end_w- start_w)/ww;

        sigma_estimates = zeros(numw,4);

        % fill in column 1 with the mid-x-value for each window
        for r=1:numw
            sigma_estimates(r,1) = start_w + (r-1)*ww + ww/2;
        end

        for r=1:numw

                % all the indices for the x-values in this window
                ind = xx_sorted > sigma_estimates(r,1) - ww/2 & xx_sorted < sigma_estimates(r,1) + ww/2;

                % the x-values that we care about
                x=xx_sorted(ind);

                % the y values that we care about
                y = yy(ind);

                % the smoothed values
                yyy = yy_smooth(ind);
                sigma_estimates(r,2) = mean(yyy);

                % number of points
                n = size(x,1);

                sm=0;
                for i=1:n
                    sm = sm + (yyy(i)-y(i))*(yyy(i)-y(i));
                end

                sigma_estimator = sqrt(sm/n);

                sigma_estimates(r,3) = sigma_estimator;

                % now record the x-value that gets you the .975 quantile of
                % a normal distribution with this sigma
                pred_max_dev = norminv(qtile,0.0,sigma_estimator);

                sigma_estimates(r,4) = pred_max_dev;


        end

        AA = sigma_estimates;

