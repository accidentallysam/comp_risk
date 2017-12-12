%%%
% MAIN OUTPUT FUNCTION
%

function [] = output_survivor(X,param_hat,cov_labels,fval,exitflag,output,s,H,predict,pSE)
    % store some basic dimensionality information
    J = length(param_hat);
    T = arrayfun(@(x) length(x.bhaz),param_hat);
    K = arrayfun(@(x) length(x.b),param_hat);

    % convert the parameters into matrix form so that we can easily
    % calculate test statistics
    param_ary = param2mtx(param_hat);
    
    % calculate the asymptotic variance (times 1/n), standard errors, and
    % t-statistics
    Hinv = inv(H);
    Avar = diag(Hinv);
    SE = sqrt(Avar);
    t = param_ary./SE;
    
    % convert statistics into structure (rather than matrix) form
    SE = mtx2param(SE,J,T,K);
    t = mtx2param(t,J,T,K);
    
    % plot the survivor functions and confidence intervals
    output_survivorplot(X,param_hat,cov_labels,SE,predict,pSE);
end

%%%
% OTHER FORMATTING AND OUTPUT FUNCTIONS
%

function [] =  output_survivorplot(X,param_hat,cov_labels,SE,predict,pSE)
    % store some basic dimensionality information
    J = length(param_hat);
    T = max(arrayfun(@(x) length(x.bhaz),param_hat));

    % now do plots just for the average predicted survivor function
    for j = 1:J
        ymax = 0;
        N = size(X,1);
        K = size(X{1,j},2);
        % calculate the median income
        Xmat = cell2mat(X(:,j));
        icov = min(find(cellfun(@(x) ~isempty(strfind(x, 'linc')),...
            cov_labels(:,j))));
        med = median(Xmat(:,icov));
        index = Xmat(:,icov)>=med;
        for k = 1:2
            p = mean(reshape(predict(index,j),T,[]),2);
            pstde = sqrt((1/N)*sum(reshape(pSE(index,j).^2,T,[]),2));
            pmarg = abs(norminv((1-0.95)/2,0,1)).*pstde;
            h = plot_conf(100+j,p,pmarg,[0 0 0]);
            ymax = max(ymax, max(ylim(gca)));
            xlabel('Time (days)');
            ylabel('-log(-log(S(t)))');
            set(gca,'XLim',[0 log(length(p))]);
            set(gca,'XTick',log(1:length(p)));
            set(gca,'XTickLabel',1:length(p));
            index = ~index;
            hold on
        end
        set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        % Import the STATA file
        newData1 = importdata('data/Sinc_illn.csv');

        % Break the data up into a new structure with one field per column.
        colheaders = genvarname(newData1.colheaders);
        for i = 1:length(colheaders)
            dataByColumn1.(colheaders{i}) = newData1.data(:, i);
        end

        % Create new variables in the base workspace from those fields.
        vars = fieldnames(dataByColumn1);
        for i = 1:length(vars)
            eval([vars{i} '= dataByColumn1.(vars{i});']);
        end
        t = 1:T;
        x = log(t);
        if j==1
            line(x,llS1(find(hinc)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS1(find(~hinc)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        else
            line(x,llS2(find(hinc)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS2(find(~hinc)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        end
        hold off
        legend('Predicted','Empirical','Location','NorthEast');
    end
    % now do plots just for the average predicted survivor function
    for j = 1:J
        ymax = 0;
        N = size(X,1);
        K = size(X{1,j},2);
        % calculate the median distance
        Xmat = cell2mat(X(:,j));
        icov = min(find(cellfun(@(x) ~isempty(strfind(x, 'dist')),...
            cov_labels(:,j))));
        med = median(Xmat(:,icov));
        index = Xmat(:,icov)>=med;
        for k = 1:2
            p = mean(reshape(predict(index,j),T,[]),2);
            pstde = sqrt((1/N)*sum(reshape(pSE(index,j).^2,T,[]),2));
            pmarg = abs(norminv((1-0.95)/2,0,1)).*pstde;
            h = plot_conf(200+j,p,pmarg,[0 0 0]);
            ymax = max(ymax, max(ylim(gca)));
            xlabel('Time (days)');
            ylabel('-log(-log(S(t)))');
            set(gca,'XLim',[0 log(length(p))]);
            set(gca,'XTick',log(1:length(p)));
            set(gca,'XTickLabel',1:length(p));
            index = ~index;
            hold on
        end
        set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        % Import the STATA file
        newData1 = importdata('data/Sdist_illn.csv');

        % Break the data up into a new structure with one field per column.
        colheaders = genvarname(newData1.colheaders);
        for i = 1:length(colheaders)
            dataByColumn1.(colheaders{i}) = newData1.data(:, i);
        end

        % Create new variables in the base workspace from those fields.
        vars = fieldnames(dataByColumn1);
        for i = 1:length(vars)
            eval([vars{i} '= dataByColumn1.(vars{i});']);
        end
        t = 1:T;
        x = log(t);
        if j==1
            line(x,llS1(find(hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS1(find(~hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        else
            line(x,llS2(find(hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS2(find(~hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        end
        hold off
        legend('Predicted','Empirical','Location','NorthEast');
    end
    % now do plots just for the average predicted survivor function
    for j = 1:J
        ymax = 0;
        N = size(X,1);
        K = size(X{1,j},2);
        % calculate the median distance
        Xmat = cell2mat(X(:,j));
        icov = min(find(cellfun(@(x) ~isempty(strfind(x, 'kno')),...
            cov_labels(:,j))));
        med = median(Xmat(:,icov));
        index = Xmat(:,icov)>=med;
        for k = 1:2
            p = mean(reshape(predict(index,j),T,[]),2);
            pstde = sqrt((1/N)*sum(reshape(pSE(index,j).^2,T,[]),2));
            pmarg = abs(norminv((1-0.95)/2,0,1)).*pstde;
            h = plot_conf(300+j,p,pmarg,[0 0 0]);
            ymax = max(ymax, max(ylim(gca)));
            xlabel('Time (days)');
            ylabel('-log(-log(S(t)))');
            set(gca,'XLim',[0 log(length(p))]);
            set(gca,'XTick',log(1:length(p)));
            set(gca,'XTickLabel',1:length(p));
            index = ~index;
            hold on
        end
        set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        % Import the STATA file
        newData1 = importdata('data/Skno_illn.csv');

        % Break the data up into a new structure with one field per column.
        colheaders = genvarname(newData1.colheaders);
        for i = 1:length(colheaders)
            dataByColumn1.(colheaders{i}) = newData1.data(:, i);
        end

        % Create new variables in the base workspace from those fields.
        vars = fieldnames(dataByColumn1);
        for i = 1:length(vars)
            eval([vars{i} '= dataByColumn1.(vars{i});']);
        end
        t = 1:T;
        x = log(t);
        if j==1
            line(x,llS1(find(hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS1(find(~hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        else
            line(x,llS2(find(hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS2(find(~hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        end
        hold off
        legend('Predicted','Empirical','Location','NorthEast');
    end
    % now do plots just for the average predicted survivor function
    for j = 1:J
        ymax = 0;
        N = size(X,1);
        K = size(X{1,j},2);
        % calculate the median distance
        Xmat = cell2mat(X(:,j));
        icov = min(find(cellfun(@(x) ~isempty(strfind(x, 'male')),...
            cov_labels(:,j))));
        med = median(Xmat(:,icov));
        index = Xmat(:,icov)>=med;
        for k = 1:2
            p = mean(reshape(predict(index,j),T,[]),2);
            pstde = sqrt((1/N)*sum(reshape(pSE(index,j).^2,T,[]),2));
            pmarg = abs(norminv((1-0.95)/2,0,1)).*pstde;
            h = plot_conf(400+j,p,pmarg,[0 0 0]);
            ymax = max(ymax, max(ylim(gca)));
            xlabel('Time (days)');
            ylabel('-log(-log(S(t)))');
            set(gca,'XLim',[0 log(length(p))]);
            set(gca,'XTick',log(1:length(p)));
            set(gca,'XTickLabel',1:length(p));
            index = ~index;
            hold on
        end
        set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off');
        % Import the STATA file
        newData1 = importdata('data/Smale_illn.csv');

        % Break the data up into a new structure with one field per column.
        colheaders = genvarname(newData1.colheaders);
        for i = 1:length(colheaders)
            dataByColumn1.(colheaders{i}) = newData1.data(:, i);
        end

        % Create new variables in the base workspace from those fields.
        vars = fieldnames(dataByColumn1);
        for i = 1:length(vars)
            eval([vars{i} '= dataByColumn1.(vars{i});']);
        end
        t = 1:T;
        x = log(t);
        if j==1
            line(x,llS1(find(hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS1(find(~hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        else
            line(x,llS2(find(hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
            line(x,llS2(find(~hdist)),'LineWidth',2,'Color',[0.7 0.7 0.7],'LineStyle','--');
        end
        hold off
        legend('Predicted','Empirical','Location','NorthEast');
    end
end

function [h] = plot_conf(j,val,pmarg,color)
    figure(j);
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
    % confidence invervals
    t = 1:length(val);
%     x(t*2)=t; x(t*2-1)=t-1;
%     y(t*2)=upper; y(t*2-1)=upper;
%     z(t*2)=lower; z(t*2-1)=lower;
    x = log(t);
    % fill([x fliplr(x)],[upper' fliplr(lower')],color);
    % errorbar(x,val,pmarg,'LineWidth',1,'Color',color,'LineStyle','-');
    h = line(x,val,'LineWidth',2,'Color',color,'LineStyle','-');
end
