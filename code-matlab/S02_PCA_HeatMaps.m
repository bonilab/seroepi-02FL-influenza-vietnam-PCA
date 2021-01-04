%exported as 14" wide, 7" tall, 600dpi, fixed fonts at 7
clear all;

load mfb_workspace_titermatrix.mat;

%%


lvls=128;
ind=lvls:-1:1;

%mymap = colormap(1-pink(256));  black=0; white=1;
%mymap = pink(lvls);
mymap = bone(lvls);
%mymap = copper(lvls);
%mymap = hot(lvls);
mymap2 = mymap(ind,:); % reverses the map
mymap2 = mymap2(1:(lvls),:);
%mymap2 = [1 1 1; mymap2];

pc1max = max(SCORE(:,1)); % this is 11.7059
pc1min = min(SCORE(:,1)); % this is -13.1364
pc2max = max(SCORE(:,2)); % this is 10.0543
pc2min = min(SCORE(:,2)); % this is -9.9370

xp=1.13;

bnds = [pc1min pc1max pc2min pc2max];

for spi=1:12

    subplot(3,6,spi)
    D=SCORE( B(:,12) >= spi-1 & B(:,12) < spi, : );    % col 12 is age
    
    %DataDensityPlot(D(1:180,1), D(1:180,2), lvls, mymap2, bnds);
    DataDensityPlot(D(:,1), D(:,2), lvls, mymap2, bnds); hold on;
    
    plot( [136 136], [0 256], '-', 'Color', [0 0 0 0.3]); hold on;
    plot( [0 256], [127 127], '-', 'Color', [0 0 0 0.3]);
    if spi==1
        text(205,235,sprintf('%1.1f-%d',spi-0.5,spi));
    else
        text(205,235,sprintf('%d-%d',spi-1,spi));
    end
    
    axis xy;
    box off;
    axis([-5 256 1 256]);
    
    if spi==6 || spi==12
        %colorbar;
    end

    if spi==1 || spi==7
        %ylabel('PC2');
        op=get(gca,'OuterPosition');
        set(gca,'OuterPosition',op.*[1 1 xp xp]);
        text(-45,117,'PC2','Rotation',90);
    else
        op=get(gca,'OuterPosition');
        set(gca,'OuterPosition',op.*[1 1 xp xp]);
    end

end

z=[12 14 16 18 30 50 100];
for spi=13:18

    subplot(3,6,spi)
    D=SCORE( B(:,12) >= z(spi-12) & B(:,12) < z(spi-11), : );

    %DataDensityPlot(D(1:180,1), D(1:180,2), lvls, mymap2, bnds);
    DataDensityPlot(D(:,1), D(:,2), lvls, mymap2, bnds);  hold on;
    plot( [136 136], [0 256], '-', 'Color', [0 0 0 0.3]); hold on;
    plot( [0 256], [127 127], '-', 'Color', [0 0 0 0.3]);
    if spi==18
        text(205,235,sprintf('>%d',z(spi-12)));
    else
        text(205,235,sprintf('%d-%d',z(spi-12),z(spi-11)));
    end
    
    if spi==13 
        %ylabel('PC2');
        text(-45,117,'PC2','Rotation',90);
    end
    text(118,-20,'PC1');
    %xlabel('PC1');
    
    axis xy;
    box off;
    
    axis([-5 256 1 256]);
    op=get(gca,'OuterPosition');
    set(gca,'OuterPosition',op.*[1 1 xp xp]);
end

%colorbar;


%D=SCORE( B(:,12) >= 0.0 & B(:,12) < 1.0, : );
%ff=DataDensityPlot(D(:,1), D(:,2), lvls, mymap2, bnds);
%D=SCORE( B(:,12) >= 6.0 & B(:,12) < 7.0, : );
%nr=size(D,1);



