function [result_vector] = data_analysis(group1,group2,param_type,fit_indx)


if param_type==1 % Continuous Varaible
    
    % Check for normality
    % Report mean+-std or median with IQR
    % Compare values across groups
    
   %nt=normalitytest(group1');
   nt(7,3) = 0;
    
    if nt(7,3)==0
        
        % Report median and IQR
        val1a=nanmedian(group1);
        val1b=prctile(group1,25);
        val1c=prctile(group1,75);
        
        val2a=nanmedian(group2);
        val2b=prctile(group2,25);
        val2c=prctile(group2,75);        
        
        % Do a non-parametric comparison
        if numel(find(~isnan(group1)))*numel(find(~isnan(group2))) > 0
            %grp_vec=[zeros(numel(group1),1); ones(numel(group2),1); 2*ones(numel(group3),1)];
            for i = 1:(length(group1)+length(group2))
                if i <= length(group1)
                    grp_vec{i} = '1';
                else
                    grp_vec{i} = '2';
                end
            end
            gen_vec=[group1; group2];

            if fit_indx == 1
                [r,p] = corr(group1,group2,'Type','Spearman'); %Spearman correlation
                fit = polyfit(group1,group2,1); %Linear fit
            elseif fit_indx == 2
                %mdlr = fitlm(group1,group2,'RobustOpts','on');
                mdlr = fitlm(group1,group2);
                fit = mdlr.Coefficients.Estimate;
                r = sqrt(mdlr.Rsquared.Ordinary);
                p = mdlr.Coefficients.pValue(2);
            end
          
 %           [comp1,~,stats]=ranksum(group1,group2); %wilcoxon rank sum test for non-parametric groups
            
        else
            comp1=NaN;
        end
        
    elseif nt(7,3)==1
        
        % Report mean and std
        val1a=mean(group1);
        val1b=std(group1);
        val1c=NaN;
        
        val2a=mean(group2);
        val2b=std(group2);
        val2c=NaN;
        
        % Do a unpaired t-test
        disp('P')
        %grp_vec=[zeros(numel(group1),1); ones(numel(group2),1); 2*ones(numel(group3),1)];
        for i = 1:(length(group1)+length(group2))
            if i <= length(group1)
                grp_vec{i} = 'no dyssync';
            else
                grp_vec{i} = 'dyssync';
            end
        end
   
        gen_vec=[group1; group2];
        if fit_indx == 1
            [r,p] = corr(group1,group2,"type","Pearson"); %Pearson correlation
            fit = polyfit(group1,group2,1); %Linear fit

        elseif fit_indx == 2
            mdlr = fitlm(group1,group2,'RobustOpts','on');
            fit = mdlr.Coefficients.Estimate;
            r = sqrt(mdlr.Rsquared.Ordinary);
            p = mdlr.Coefficients.pValue(2);
        end
        
        

 %       [~,comp1,~,stats]=ttest(group1,group2);

    end
    
  
    
    
    result_vector=[nt(7,3) val1a val1b val1c val2a val2b val2c r p fit(1) fit(2)];
    
elseif param_type==2 % Categorical
    
    [result_vector] = prop_table_analysis(group1,group2);
    
end