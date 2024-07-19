function [fig] = plotHistogram(measure,title_label)

count = 1;
fig = figure;
for pat = [6,8,3,2,5,4,1,7]%pat = [6,3,5,1,8,2,4,7]

    subplot(4,2,count)
    histogram(measure(:,pat),'Normalization','probability')
    cov = std(measure(:,pat))./mean(measure(:,pat));
    text(0,0.8,['COV = ',num2str(cov,'%.2f')])
    if count == 1
        title('Responders')
    elseif count == 2
        title('Nonresponders')
    end
    count = count + 1;
    ylim([0 1])
end
sgtitle(title_label)
end