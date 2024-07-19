function [fig,mean_results] = plotResponderBoxplotAllEsimates(resultPSM,resultCT,ymin,ymax,ytitle)


CRTrespondersPSM = [resultPSM(6); resultPSM(3); resultPSM(5); resultPSM(1)];
CRTnonrespondersPSM = [resultPSM(8); resultPSM(2); resultPSM(4); resultPSM(7)];
for i = 1:5
    CRTrespondersCT(:,i) = [resultCT(i,6); resultCT(i,3); resultCT(i,5); resultCT(i,1)];
    CRTnonrespondersCT(:,i) = [resultCT(i,8); resultCT(i,2); resultCT(i,4); resultCT(i,7)];
end
input = [CRTrespondersPSM,CRTnonrespondersPSM,CRTrespondersCT(:,1),...
    CRTnonrespondersCT(:,1),CRTrespondersCT(:,2),CRTnonrespondersCT(:,2),...
    CRTrespondersCT(:,3),CRTnonrespondersCT(:,3),CRTrespondersCT(:,4),...
    CRTnonrespondersCT(:,4),CRTrespondersCT(:,5),CRTnonrespondersCT(:,5)];
labels = {'R','NonR','R','NonR','R','NonR','R','NonR','R','NonR','R','NonR'};
%labels = {'Resp','Non-Resp','Resp','Non-Resp','Resp','Non-Resp','Resp','Non-Resp','Resp','Non-Resp','Resp','Non-Resp'};
%labels = {'Responders','Non-Responders','Responders','Non-Responders','Responders','Non-Responders','Responders','Non-Responders','Responders','Non-Responders','Responders','Non-Responders'};
fig = figure; hold all
set(fig, 'Position',[100 300 1800 600])
boxplot(input,'Labels',labels,'Colors','k');
h=findobj('LineStyle','--'); set(h, 'LineStyle','-');
ylabel(ytitle)
ylim([ymin ymax])
xlim([0.25 12.75])
ytext = ymax + 0.043.*(ymax-ymin);
text(1.5,ytext,'SSA_P_S_M','HorizontalAlignment','center')
text(3.5,ytext,'P_L_H_CSA','HorizontalAlignment','center')
text(5.5,ytext,'WS_E_DSA','HorizontalAlignment','center')
text(7.5,ytext,'WS_T_VSA','HorizontalAlignment','center')
text(9.5,ytext,'P_g_e_nSA','HorizontalAlignment','center')
text(11.5,ytext,'P_g_e_n_,_s_c_a_l_e_dSA','HorizontalAlignment','center')
xline(2.5,'--','Color',[0.8 0.8 0.8])
xline(4.5,'--','Color',[0.8 0.8 0.8])
xline(6.5,'--','Color',[0.8 0.8 0.8])
xline(8.5,'--','Color',[0.8 0.8 0.8])
xline(10.5,'--','Color',[0.8 0.8 0.8])

yline_start = ymax - 0.1.*(ymax-ymin);
yline_edge = ymax - 0.12.*(ymax-ymin);
text_line = ymax - 0.0714.*(ymax-ymin);
for i = 1:2:length(input)
    [~,p,~] = ttest2(input(:,i),input(:,i+1));
    if p < 0.05
        line([i i+1],[yline_start yline_start],'Color','k')
        line([i i],[yline_start+0.001 yline_edge],'Color','k')
        line([i+1 i+1],[yline_start+0.001 yline_edge],'Color','k')
        text(i+0.5,text_line,'*','HorizontalAlignment','center')
    end

end



 set(findall(gcf,'-property','FontSize'),'FontSize',20)
 set(findall(gcf,'-property','LineWidth'),'LineWidth',2)
 set(findall(gcf,'-property','MarkerSize'),'MarkerSize',10)
for z = 1:length(input)
    mean_results(z,:) = [mean(input(:,z)) std(input(:,z))];
end
 
end