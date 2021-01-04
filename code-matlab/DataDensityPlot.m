function [ f ] = DataDensityPlot( x, y, levels, mymap, bnds )
%DATADENSITYPLOT Plot the data density 
%   Makes a contour map of data density
%   x, y - data x and y coordinates
%   levels - number of contours to show
%
% By Malcolm Mclean
%
    map = dataDensity(x, y, bnds, 256, 256 ); % original function call looks like this
    %map = dataDensity(x, y, bnds, 256, 256, bnds, 0);
    map = map - min(min(map));
    map = floor(map ./ max(max(map)) * (levels-1));
    %f = figure();
    
    imagesc(map);
    %colormap(jet(levels));
    colormap(mymap);
    set(gca, 'XTick', [1 256]);
    set(gca, 'XTickLabel', {'', ''});
    set(gca, 'YTick', [1 256]);
    set(gca, 'YTickLabel', {'', ''});
    %uiwait;
    f=map;  
end

