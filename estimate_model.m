%%%
% MAIN ESTIMATION FUNCTION
%

function [param_hat,fval,exitflag,output,s,H] = ...
        estimate_model(exit_state, spell_len, X)
    % get the dimensionality of the parameter space
    pdim = cellfun(@(x) size(x), X(1,:), 'UniformOutput', false);
    % set the initial paramter values
    param0 = struct(...
        'bhaz',{0.001*ones(pdim{1}(1),1),...
            0.001*ones(pdim{2}(1),1)},...
        'b',{0.001*ones(pdim{1}(2),1),...
            0.001*ones(pdim{2}(2),1)}...
        );
    % estimate the model
    [param_hat,fval,exitflag,output,s,H] = ...
        max_likelihood(exit_state,spell_len,X,param0);
end


%%%
% MAXIMUM LIKELIHOOD FUNCTIONS
%

function [param_hat,fval,exitflag,output,s,H] = ...
        max_likelihood(j,t,X,param0)
    global iter
    % store the number of exit states
    J = size(X,2);
    % store the maximum number of time periods observed
    % (this is the number of periods for which we 
    % estimate a hazard)
    T = zeros(1,J);
    for i = 1:J T(i) = max(t); end % = max(t(j==i)); end
    % store the number of covariates for each exit state
    K = cellfun(@(x) size(x,2),X(1,:));
    % set the optimization options
    % options = optimoptions('fminunc','GradObj','on');
    options = optimset('GradObj','on',...
        'Hessian','on',...
        'LargeScale','on',...
        'TolX',1e-10...
        );
    % get the coefficient estimates
    iter = [];
    [param_hat,fval,exitflag,output,s,H] = fminunc(...
        @(param) nlog_likelihood(j,t,X,mtx2param(param,J,T,K)),...
            param2mtx(param0),options...
        );
    param_hat = mtx2param(param_hat,J,T,K);
end

function [tau] = bhaz(t,param)
    tau = param.bhaz(t);
end

function [phi] = mult(X,param)
   phi = exp(X*param.b);
end

function [z] = inthaz(t,X,param)
    J = size(X,2);
    z = arrayfun(@(x) bhaz(1:t,param(x))'*mult(X{x}(1:t,:),...
            param(x)),1:J);
end

function [h] = haz(t,X,param)
    J = size(X,2);
    h = cell2mat(arrayfun(@(x) bhaz(t,param(x)).*...
            mult(X{x}(t,:),param(x)),1:J,'UniformOutput',false));
end

function [a] = alpha(t,X,param)
   a = log(inthaz(t,X,param));
end

function [y] = f(x)
   y = exp(x).*exp(-exp(x));
end

function [y] = F(x)
   y = -exp(-exp(x));
end

function [l] = L(j,t,X,param)
    % store the number of exit states
    J = size(X,2);
    if (j<=J) % uncensored
        % calculate the integrated hazards
        inthaz_tm1 = inthaz(t-1,X,param);
        inthaz_t   = inthaz(t,X,param);
        % get the censored exit state
        k = setdiff(1:length(inthaz_tm1),j);
        % calcuate the likelihood
        l = (exp(-inthaz_tm1(k)-inthaz_tm1(j))-...
                exp(-inthaz_t(j)-inthaz_t(k))).*...
                    (-inthaz_t(j)+inthaz_tm1(j))./...
                        (-inthaz_t(k)+inthaz_tm1(k)...
                            -inthaz_t(j)+inthaz_tm1(j));
    else % censored
        % calculate the integrated hazard
        inthaz_t = inthaz(t,X,param);
        % calcuate the likelihood
        l = exp(-sum(inthaz_t));
    end % end if
end

function [s] = score(j,t,X,param)
    % store the number of exit states and number of time periods
    J = size(X,2);
    % get the censored exit state(s)
    k = setdiff(1:J,j);
    % construct the structure used to store the score contributions
    sparam = param;
    % zero out score contributions to start
    for l=1:J
        sparam(l).bhaz = zeros(size(sparam(l).bhaz));
        sparam(l).b    = zeros(size(sparam(l).b));
    end % end for
    if (j<=J) % uncensored
        T = size(X{j},1);
        % calculate the hazard and the integrated hazard
        myhaz    = haz(1:T,X,param);
        % calculate constants in the score function
        A = myhaz(t,k)+myhaz(t,j);
        B = exp(A);
        C = -1/A+B/(B-1)+1/myhaz(t,j);
        m = arrayfun(@(x) mult(X{j}(x,:),param(j)),1:T)';
        % calculate the score for the baseline hazard coeffs first
        sparam(j).bhaz(t)=C*m(t);
        sparam(j).bhaz = sparam(j).bhaz-...
            [m(1:t);zeros(length(sparam(j).bhaz)-t,1)];
        % calculate the score for the regression coeffs
        sparam(j).b = C*X{j}(t,:)'*myhaz(t,j);
        sparam(j).b = sparam(j).b-...
            arrayfun(@(x) (bhaz(1:t,param(j)).*X{j}(1:t,x))'*...
                mult(X{j}(1:t,:),param(j)),1:size(X{j},2))';
    else % censored
        for l=1:J
            T = size(X{l},1);
            % calculate constants in the score function
            m = arrayfun(@(x) mult(X{l}(x,:),param(l)),1:T)';
            % calculate the score for the baseline hazard coeffs
            % first
            sparam(l).bhaz = ...
                -[m(1:t);zeros(length(sparam(l).bhaz)-t,1)];
            % calculate the score for the regression coeffs
            sparam(l).b = -arrayfun(@(x) ...
                (bhaz(1:t,param(l)).*X{l}(1:t,x))'*...
                    mult(X{l}(1:t,:),param(l)),1:size(X{l},2))';
        end % end for
    end % end if
    s = param2mtx(sparam);
end

function [logl grad H] = nlog_likelihood(j,t,X,param)
    global iter
    if isempty(iter), iter=0; end
    iter = iter+1;
    n = size(X,1);
    logl = 0;
    s = [];
    grad = zeros(size(param2mtx(param)));
    H = zeros(length(grad));
    for m = 1:n
        logl = logl + log(L(j(m),t(m),X(m,:),param));
        s= score(j(m),t(m),X(m,:),param);
        grad = grad + s;
        H = H + s*s';
    end % end for
    fprintf('iteration %03d, log likelihood = % 5.10f\n',iter,logl);
    logl = -logl;
    grad = -grad;
end

