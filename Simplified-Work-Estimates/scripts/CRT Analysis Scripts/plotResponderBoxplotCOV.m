function [fig] = plotResponderBoxplotCOV(resultPSM,resultCT,location,Volume,Segments)

CRTrespondersPSM = [resultPSM{6}, resultPSM{3}, resultPSM{5}, resultPSM{1}];
CRTnonrespondersPSM = [resultPSM{8}, resultPSM{2}, resultPSM{4}, resultPSM{7}];
fig = figure; 
set(fig, 'Position',[100 300 1600 600])
s1 = subplot(1,3,1);
boxplot([CRTrespondersPSM',CRTnonrespondersPSM'],'Labels',{'Responders','Non-Responders'})
title('PSM Result')
axis('square')
%ylim([1 1.8])
[~,p,~] = ttest2(CRTrespondersPSM',CRTnonrespondersPSM');
%Xrel = interp1([0 1], xlim(),0.05);
Yrel = interp1([0 1], ylim(),0.9);
text(2,Yrel,['p = ',num2str(p,'%.2f')],'HorizontalAlignment','center') %2, 0.25
if Volume == 1
    ylabel(['COV of ',location, ' PSM Work Density'])
% elseif Volume == 1 && Segments == 1
%     ylabel(['COV of Segmental',location, ' PSM Work Density'])
elseif Volume == 0 %whole chamber, just negative work locations
    ylabel(['COV of ',location, ' PSM Work'])
% elseif Volume == 0 && Segments == 1
%     ylabel(['COV of Segmental',location, ' PSM Work'])
end
s1pos = get(s1,'position');

s2 = subplot(1,3,2);
CRTrespondersCT = [resultCT{6}, resultCT{3}, resultCT{5}, resultCT{1}];
CRTnonrespondersCT = [resultCT{8}, resultCT{2}, resultCT{4}, resultCT{7}];
boxplot([CRTrespondersCT',CRTnonrespondersCT'],'Labels',{'Responders','Non-Responders'})
title('CT Result')
%ylabel([{'Volume Fraction of '},{[location, ' Negative Work']}])
axis('square')
%ylim([1 1.8])
[~,p,~] = ttest2(CRTrespondersCT',CRTnonrespondersCT');
%Xrel = interp1([0 1], xlim(),0.05);
Yrel = interp1([0 1], ylim(),0.9);
text(2,Yrel,['p = ',num2str(p,'%.2f')],'HorizontalAlignment','center') %2,0.25
if Volume == 1 && Segments == 0
    ylabel(['COV of ',location, ' CT Work Density'])
elseif Volume == 1 && Segments == 1
     ylabel(['COV of Segmental ',location, ' CT Work Density'])
elseif Volume == 0 && Segments == 0 %whole chamber, just negative work locations
    ylabel(['COV of ',location, ' CT Work'])
elseif Volume == 0 && Segments == 1
    ylabel(['COV of Segmental ',location, ' CT Work'])
end
s2pos = get(s2,'position');
s2pos(3:4) = s1pos(3:4);
set(s2,'position',s2pos)


s3 = subplot(1,3,3);
plot(CRTrespondersPSM,CRTrespondersCT,'.b'); hold on;
plot(CRTnonrespondersPSM,CRTnonrespondersCT,'.r'); hold on;
xlimit = xlim();
ylimit = ylim();
% range = max(xlimit(2),ylimit(2)) - min(xlimit(1),ylimit(1));
% minVal = round((min(xlimit(1),ylimit(1))-0.15*range)*10)/10;
% maxVal = round((max(xlimit(2),ylimit(2))+0.15*range)*10)/10;
% limits = [minVal,maxVal];
limits = [min(xlimit(1),ylimit(1)) - 0.05, max(xlimit(2),ylimit(2)) + 0.05];
plot(limits,limits,'k')
axis('square')
xlabel('PSM Result')
ylabel('CT Result')
title('PSM and CT Result Agreement')
legend('Responders','Nonresponders','location','northwest')
%r = corr([CRTrespondersPSM,CRTnonrespondersPSM]',[CRTrespondersCT,CRTnonrespondersCT]','type','Spearman');
R2Val = 1-det(corrcoef([CRTrespondersPSM,CRTnonrespondersPSM]',[CRTrespondersCT,CRTnonrespondersCT]'));
Xrel = interp1([0 1], xlim(),0.1);
Yrel = interp1([0 1], ylim(),0.8);
text(Xrel,Yrel,['r^2 = ',num2str(R2Val,'%.2f')])

%Patient # Labels
xlimits = xlim();
xInterval = (xlimits(2) - xlimits(1))/20;
labelOffset = xInterval/2;
CRTresponder_label = [6,3,5,1];
CRTnonresponder_label = [8,2,4,7];
patlabel = [CRTresponder_label,CRTnonresponder_label];
data_x = [CRTrespondersPSM,CRTnonrespondersPSM];
data_y = [CRTrespondersCT,CRTnonrespondersCT];
for k = 1:8
    text(data_x(k)+labelOffset,data_y(k), ['BiV',num2str(patlabel(k))],...
				'horizontal','left','vertical','middle')
end

s3pos = get(s3,'position');
s3pos(1:2) = s2pos(1:2) + (s2pos(1:2) - s1pos(1:2));
s3pos(3:4) = s1pos(3:4);
set(s3,'position',s3pos)

s1pos(1:2) = s1pos(1:2) - ((s2pos(1:2) - s1pos(1:2)).*0.1);
s2pos(1:2) = s2pos(1:2) - ((s2pos(1:2) - s1pos(1:2)).*0.05);
set(s1,'position',s1pos)
set(s2,'position',s2pos)

if Volume == 1 && Segments == 0
    sgtitle(['COV of ',location,' Work Density in CRT Responders vs Non-responders'])
elseif Volume == 1 && Segments == 1
    sgtitle(['COV of Segmental ',location,' Work Density in CRT Responders vs Non-responders'])
elseif Volume == 0 && Segments == 0
    sgtitle(['COV of ',location,' Work in CRT Responders vs Non-responders'])
elseif Volume == 0 && Segments == 1
    sgtitle(['COV of Segmental ',location,' Work in CRT Responders vs Non-responders'])
end

set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf,'-property','LineWidth'),'LineWidth',3)
set(findall(gcf,'-property','MarkerSize'),'MarkerSize',35)


end