function [fig] = plotCorrelationPatientData(patientdata,workmetricPSM,titlelabel,xtitle,ytitle,ymin,ymax,xmin,xmax)

%numplots = 2;

fig = figure; set(gcf, 'Position',[100 300 400 400])
for pat = 1:8
    if pat == 1 || pat ==  3 ||pat ==  5 ||pat ==  6
        markeredge = 'k';
        markerfill = 'w';
    elseif pat == 2 ||pat ==  4 ||pat ==  8
        markeredge = [0.5 0.5 0.5];
        markerfill = [0.5 0.5 0.5];
    elseif pat == 7
        markeredge = 'k';
        markerfill = 'k';
    end
    if pat == 1 || pat == 6 || pat == 7
        markershape = 'o';
    elseif pat == 3 || pat == 8
        markershape = '^';
    elseif pat == 2 || pat == 4 || pat == 8
        markershape = 'diamond';
    end

    hold all
%    subplot(1,numplots,1); hold all
    plot(workmetricPSM(pat),patientdata(pat),markershape,'MarkerFaceColor',markerfill,...
        'MarkerEdgeColor',markeredge)
    
    if pat == 8
        %xRange = max(workmetricPSM) - min(workmetricPSM);
        minXVal = xmin;%0;%round((min(workmetricPSM)-0.4*xRange)*10)/10;
	    maxXVal = xmax;%0.8;%round((max(workmetricPSM)+0.4*xRange)*10)/10;
        xRange = maxXVal - minXVal;
        xInterval = (maxXVal - minXVal)/20;
        ylabel(ytitle)
        xlabel(xtitle)
        R2Val = 1-det(corrcoef(workmetricPSM,patientdata));
        xtext = minXVal + 0.08.*xRange;
        ytext = ymax - 0.08.*(ymax-ymin);
        text(xtext,ytext,['R^2 = ',num2str(R2Val,'%.2f')])
        p=polyfit(workmetricPSM,patientdata,1);
	    yfit = polyval(p,minXVal:xInterval:maxXVal);
        plot(minXVal:xInterval:maxXVal,yfit,'Color',[129./255 208./255 200./255]);
        ylim([ymin ymax])
        title(titlelabel)
        xlim([minXVal maxXVal])
        
    end
    
    % 
    % subplot(1,numplots,2); hold all
    % plot(workmetricCT(pat),patientdata(pat),'o','MarkerFaceColor',markerfill,'MarkerEdgeColor',markeredge)
    % if pat == 8
    %     xRange = max(workmetricCT) - min(workmetricCT);
    %     minXVal = round((min(workmetricCT)-0.4*xRange)*10)/10;
	%     maxXVal = round((max(workmetricCT)+0.4*xRange)*10)/10;
    %     xInterval = (maxXVal - minXVal)/20;
    %     ylabel('\DeltaESV (mL)')
    %     xlabel([location,' Negative Work Fraction, CT'])
    %     R2Val = 1-det(corrcoef(workmetricCT,patientdata));
    %     xtext = minXVal + 0.15.*xRange;
    %     text(xtext,80,['R^2 = ',num2str(R2Val,'%.2f')])
    %      p=polyfit(workmetricCT,patientdata,1);
	%     yfit = polyval(p,minXVal:xInterval:maxXVal);
    %     plot(minXVal:xInterval:maxXVal,yfit);
    %     ylim([-100 100])
    % end

end

 set(findall(gcf,'-property','FontSize'),'FontSize',20)
 set(findall(gcf,'-property','LineWidth'),'LineWidth',1)
 set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)

end