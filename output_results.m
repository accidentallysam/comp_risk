%%%
% MAIN OUTPUT FUNCTION
%

function [] = output_results(X,param_hat,cov_labels,fval,exitflag,output,s,H,predict,pSE)
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
    
    % output the baseline hazard
    output_newline();
    output_heading('BASELINE HAZARD');
    output_bhaz(param_hat,t,SE)
    output_divline();

    % output the covariate coefficients and stardard errors
    output_newline();
    output_heading('COVARIATES');
    output_coeff(param_hat,cov_labels,t,SE);
    output_divline();
    output_wald(param_hat,cov_labels,H);
    output_divline();
    output_newline();

    % plot the hazard functions and confidence intervals
    output_hazardplot(X,param_hat,SE,predict,pSE);
end

%%%
% OTHER FORMATTING AND OUTPUT FUNCTIONS
%

function [] = output_wald(param_hat,cov_labels,H)
    % store some basic dimensionality information
    J = length(param_hat);
    T = arrayfun(@(x) length(x.bhaz),param_hat);
    K = arrayfun(@(x) length(x.b),param_hat);

    % convert the parameters to matrix form
    param_ary = param2mtx(param_hat);

    param_index = mtx2param(zeros(size(param_ary)),J,T,K);
    for j = 1:J
        param_index(j).b(:) = ...
            cellfun(@(x) ~isempty(strfind(x, '_t')), cov_labels(:,j));
    end
    nindex = find(param2mtx(param_index));
    
    % output wald stats for all covariates
    for j = 1:J
        param_index = mtx2param(zeros(size(param_ary)),J,T,K);
        param_index(j).b(:) = 1;
        index = setdiff(find(param2mtx(param_index)),nindex);
        fprintf('W(%02d)=%-10s| ',length(index),...
            [get_wald(index,param_ary,H) ' (all)  ']);
    end
    output_newline();

%     % output wald stats for age
%     for j = 1:J
%         param_index = mtx2param(zeros(size(param_ary)),J,T,K);
%         param_index(j).b(:) = ...
%             cellfun(@(x) ~isempty(strfind(x, 'age')), cov_labels(:,j));
%         index = setdiff(find(param2mtx(param_index)),nindex);
%         fprintf('W(%02d)=%-10s| ',length(index),...
%             [get_wald(index,param_ary,H) ' (age)  ']);
%     end
%     output_newline();
% 
%     % output wald stats for quantile
%     for j = 1:J
%         param_index = mtx2param(zeros(size(param_ary)),J,T,K);
%         param_index(j).b(:) = ...
%             cellfun(@(x) ~isempty(strfind(x, 'quant')), cov_labels(:,j));
%         index = setdiff(find(param2mtx(param_index)),nindex);
%         if ~isempty(index)
%             fprintf('W(%02d)=%-10s| ',length(index),...
%                 [get_wald(index,param_ary,H) ' (qua)  ']);
%             if j==J, output_newline(); end
%         end
%     end
% 
%     % output wald stats for region
%     for j = 1:J
%         param_index = mtx2param(zeros(size(param_ary)),J,T,K);
%         param_index(j).b(:) = ...
%             cellfun(@(x) ~isempty(strfind(x, 'region')), cov_labels(:,j));
%         index = setdiff(find(param2mtx(param_index)),nindex);
%         if ~isempty(index)
%             fprintf('W(%02d)=%-10s| ',length(index),...
%                 [get_wald(index,param_ary,H) ' (reg)  ']);
%             if j==J, output_newline(); end
%         end
%     end
    
    % output wald stats for prop haz test
    for j = 1:J
        param_index = mtx2param(zeros(size(param_ary)),J,T,K);
        param_index(j).b(:) = ...
            cellfun(@(x) ~isempty(strfind(x, '_t')), cov_labels(:,j));
        index = find(param2mtx(param_index));
        if ~isempty(index)
            a = param_ary(index);
            A = zeros(length(param_ary),length(index));
            A(index+((1:length(index))'-1)*length(param_ary)) = 1;
            % now calculate the wald test statistic
            W = a'*inv(A'*inv(H)*A)*a;
            % construct critical values for the t-tests
            conf = chi2inv([0.99 0.95 0.90],length(a));
            % assign stars according to significance
            x = sum(W>conf);
            stars = [repmat('*',1,x) repmat(' ',1,3-x)];
            str = sprintf(['%05.2f' stars], W);

            fprintf('W(%02d)=%s,p=%03.2f | ',length(index),str,...
                1-chi2cdf(W,length(a)));
            if j==J, output_newline(); end
        end
    end
    
    output_divline();

    % output wald stats for all covariates
    param_index = mtx2param(zeros(size(param_ary)),J,T,K);
    for j = 1:J, param_index(j).b(:) = 1; end
    index = setdiff(find(param2mtx(param_index)),nindex);
    fprintf('W(%02d)=%-10s| ',length(index),...
        [get_wald(index,param_ary,H) ' (all)  ']);
    output_newline();

%     % output wald stats for sex
%     param_index = mtx2param(zeros(size(param_ary)),J,T,K);
%     for j = 1:J, param_index(j).b(:) = ...
%             cellfun(@(x) ~isempty(strfind(x, 'male')), cov_labels(:,j)); end
%     index = setdiff(find(param2mtx(param_index)),nindex);
%     fprintf('W(%02d)=%-10s| ',length(index),...
%         [get_wald(index,param_ary,H) ' (sex)  ']);
%     output_newline();
% 
%     % output wald stats for age
%     param_index = mtx2param(zeros(size(param_ary)),J,T,K);
%     for j = 1:J, param_index(j).b(:) = ...
%             cellfun(@(x) ~isempty(strfind(x, 'age')), cov_labels(:,j)); end
%     index = setdiff(find(param2mtx(param_index)),nindex);
%     fprintf('W(%02d)=%-10s| ',length(index),...
%         [get_wald(index,param_ary,H) ' (age)  ']);
%     output_newline();
% 
%     % output wald stats for quantile
%     param_index = mtx2param(zeros(size(param_ary)),J,T,K);
%     for j = 1:J, param_index(j).b(:) = ...
%             cellfun(@(x) ~isempty(strfind(x, 'quant')), cov_labels(:,j)); end
%     index = setdiff(find(param2mtx(param_index)),nindex);
%     if ~isempty(index)
%         fprintf('W(%02d)=%-10s| ',length(index),...
%             [get_wald(index,param_ary,H) ' (qua)  ']);
%         output_newline();
%     end
% 
%     % output wald stats for region
%     param_index = mtx2param(zeros(size(param_ary)),J,T,K);
%     for j = 1:J, param_index(j).b(:) = ...
%             cellfun(@(x) ~isempty(strfind(x, 'region')), cov_labels(:,j)); end
%     index = setdiff(find(param2mtx(param_index)),nindex);
%     if ~isempty(index)
%         fprintf('W(%02d)=%-10s| ',length(index),...
%             [get_wald(index,param_ary,H) ' (reg)  ']);
%     end
%     output_newline();


    % output wald stat for the proportional hazards test
    output_divline();
    fprintf('PROP. HAZARDS TEST\n');
    output_divline();
    a = param_ary(nindex);
    A = zeros(length(param_ary),length(nindex));
    A(nindex+((1:length(nindex))'-1)*length(param_ary)) = 1;
    % now calculate the wald test statistic
    W = a'*inv(A'*inv(H)*A)*a;
    % construct critical values for the t-tests
    conf = chi2inv([0.99 0.95 0.90],length(a));
    % assign stars according to significance
    x = sum(W>conf);
    stars = repmat('*',1,x);
    str = sprintf(['%05.2f' stars ', p=%03.2f' ], W,1-chi2cdf(W,length(a)));
    fprintf('W(%02d)=%-10s   | ',length(a),str);
    output_newline();
end

function [str] = get_wald(index,param_ary,H)
    a = param_ary(index);
    A = zeros(length(param_ary),length(index));
    A(index+((1:length(index))'-1)*length(param_ary)) = 1;
    % now calculate the wald test statistic
    W = a'*inv(A'*inv(H)*A)*a;
    % construct critical values for the t-tests
    conf = chi2inv([0.99 0.95 0.90],length(a));
    % assign stars according to significance
    x = sum(W>conf);
    stars = [repmat('*',1,x) repmat(' ',1,3-x)];
    str = sprintf(['%05.2f' stars], W);
end

function [] =  output_hazardplot(X,param_hat,SE,predict,pSE)
    % store some basic dimensionality information
    J = length(param_hat);
    T = max(arrayfun(@(x) length(x.bhaz),param_hat));

    % output the hazard plots
    ymax = 0;
    for j = 1:J
        bhaz = param_hat(j).bhaz';
        stde = SE(j).bhaz';
        marg = abs(norminv((1-0.95)/2,0,1)).*stde;
        hc = plot_conf(j,bhaz,max(bhaz-marg,0),bhaz+marg,[0.9 0.9 0.9]);
        hp = plot_pred(j,param_hat,predict);
        ymax = max(ymax, max(ylim(gca)));

        xlabel('Time (days)');
        ylabel('Hazard');
        legend([hc, hp],'Baseline hazard','Mean predicted hazard');
    end
    % now do plots just for the average predicted hazard
    for j = 1:J
        N = size(X,1);
        p = mean(reshape(predict(:,j),T,[]),2);
        pstde = sqrt((1/N)*sum(reshape(pSE(:,j).^2,T,[]),2));
        pmarg = abs(norminv((1-0.95)/2,0,1)).*pstde;
        plot_conf(j+J,p,max(p-pmarg,0),p+pmarg,[0.9 0.9 0.9]);
        ymax = max(ymax, max(ylim(gca)));
        xlabel('Time (days)');
        ylabel('Mean predicted hazard');
    end
    name = {'Hazard 1', 'Hazard 2'};
    for j = 1:J
        figure(j);
        set(gcf, 'Name',name{j},'NumberTitle','off');
        ylim(gca, [0 ymax]);
        figure(j+J);
        set(gcf, 'Name',name{j},'NumberTitle','off');
        ylim(gca, [0 ymax]);
    end
end

function [h] = plot_conf(j,val,lower,upper,color)
    figure(j);
    % confidence invervals
    t = 1:length(lower);
    x(t*2)=t; x(t*2-1)=t-1;
    y(t*2)=upper; y(t*2-1)=upper;
    z(t*2)=lower; z(t*2-1)=lower;
    set(fill([x,x(end:-1:1)],...
        [y,z(end:-1:1)],color),'EdgeColor',color);
    line(x,y,'LineWidth',1,'Color',[0 0 0]);
    line(x,z,'LineWidth',1,'Color',[0 0 0]);
    line([x(1) x(1)],[y(1) z(1)],'LineWidth',1,'Color',[0 0 0]);
    line([x(end) x(end)],[y(end) z(end)],'LineWidth',1,'Color',[0 0 0]);

    % estimated values
    w(t*2)=val; w(t*2-1)=val;
    h = line(x,w,'LineWidth',2,'Color',[0 0 0]);
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
    set(gca, 'XTick',[1 2 3 4 5 6 7]);
    xlim([0 5]);
end

function [h] = plot_pred(j,param_hat,pred)
    % store some basic dimensionality information
    J = length(param_hat);
    T = max(arrayfun(@(x) length(x.bhaz),param_hat));
    t = 1:T;

    p = mean(reshape(pred(:,j),T,[]),2);
    x(t*2)=t; x(t*2-1)=t-1;
    y(t*2)=p; y(t*2-1)=p;
    figure(j);
    hold on
    h = line(x,y,'LineStyle','--','LineWidth',1,'Color',[0.2 0.2 0.2]);
    set(gca, 'FontSize', 20, 'FontName', 'Times New Roman');
    set(gca, 'XTick',[1 2 3 4 5 6 7]);
    xlim([0 5]);
    hold off
end

function [stars] = get_stars(t_ary)
    % construct critical values for the t-tests
    conf = norminv((1-[0.99 0.95 0.90])./2,0,1);
    % assign stars according to significance
    stars = arrayfun(@(x) [repmat('*',1,x) repmat(' ',1,3-x)],... 
        arrayfun(@(x) sum(abs(x)>abs(conf)),t_ary),...
        'UniformOutput',false)';
end

function [str] = format_coeff(b, t)
    str = arrayfun(@(x) sprintf(' %7.4f',x),b,'UniformOutput',false);
    str = cellfun(@(x,y) [x y],str,get_stars(t),'UniformOutput',false);
end

function [str] = format_stderr(SE)
    str = arrayfun(@(x) sprintf('(%7.4f)  ',x),SE,'UniformOutput',false);
end

function [] = output_bhaz(param_hat,t,SE)
    % store some basic dimensionality information
    J = length(param_hat);
    T = max(arrayfun(@(x) length(x.bhaz),param_hat));
    
    % construct the coefficient label and value strings
    bstr1 = format_coeff(param_hat(1).bhaz',t(1).bhaz);
    bstr2 = format_coeff(param_hat(2).bhaz',t(2).bhaz);
    bout = [bstr1; bstr2];

    % construct the standard error strings
    vstr1 = format_stderr(SE(1).bhaz');
    vstr2 = format_stderr(SE(2).bhaz');
    vout = [vstr1; vstr2];

    % construct the spacing strings
    sout = [cellstr(repmat('           f',length(SE(1).bhaz),1))';...
        cellstr(repmat('           f',length(SE(2).bhaz),1))'];
    sout = cellfun(@(x) x(1:end-1),sout,'UniformOutput',false);

    bind = cumsum(ones(1,size(bout,2)))*3-2;
    vind = bind+1;
    sind = bind+2;

    % construct the output strings
    allout(:,vind) = vout;
    allout(:,bind) = bout;
    allout(:,sind) = sout;
    
    tind = (1:T)*3-2;
    sind = setdiff(1:T*3,tind);
    tout(sind) = {'   '};
    tout(tind) = arrayfun(@(x) sprintf('%3d',x),1:T,...
        'UniformOutput',false);
    allout = [tout; allout];

    % output the coefficient information
    fprintf([' %3s | ' repmat(' %-10d | ',1,J) '\n'],'t',cumsum(ones(1,J)));
    output_divline();
    fprintf('     |             |             | \n');
    fprintf([' %3s | ' repmat('%9s | ',1,J) '\n'],allout{:,:});
end

function [] = output_coeff(param_hat,cov_labels,t,SE)
    % store some basic dimensionality information
    J = length(param_hat);

    % construct the coefficient label and value strings
    bstr1 = format_coeff(param_hat(1).b',t(1).b);
    bstr2 = format_coeff(param_hat(2).b',t(2).b);
    bout = [cov_labels(:,1)'; bstr1; cov_labels(:,2)'; bstr2];

    % construct the standard error strings
    vstr1 = format_stderr(SE(1).b');
    vstr2 = format_stderr(SE(2).b');
    vout = [cellstr(repmat(' ',length(SE(1).b),1))'; vstr1;...
        cellstr(repmat(' ',length(SE(2).b),1))'; vstr2];

    % construct the spacing strings
    sout = [cellstr(repmat(' f',length(SE(1).b),1))';...
        cellstr(repmat('           f',length(SE(1).b),1))';...
        cellstr(repmat(' f',length(SE(2).b),1))';...
        cellstr(repmat('           f',length(SE(2).b),1))'];
    sout = cellfun(@(x) x(1:end-1),sout,'UniformOutput',false);

    bind = cumsum(ones(1,size(bout,2)))*3-2;
    vind = bind+1;
    sind = bind+2;

    % construct the output strings
    allout(:,vind) = vout;
    allout(:,bind) = bout;
    allout(:,sind) = sout;

    % output the coefficient information
    fprintf([repmat(['%-9d ' repmat(' ',1,11) ' | '] ,1,J) '\n'],...
        cumsum(ones(1,J)));
    output_divline();
    fprintf('                      |                       | \n');
    fprintf([repmat('%-9s %9s | ',1,J) '\n'],allout{:,:});
end
