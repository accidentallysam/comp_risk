function [predict SE] = predicted_S(X,param_hat,H)
    fprintf('generating predicted Survivor function and standard errors...');
    J = size(X,2);
    T = max(arrayfun(@(x) length(x.bhaz),param_hat));
    K = arrayfun(@(x) length(x.b),param_hat);

    param_ary = param2mtx(param_hat);
    Hinv = inv(H);
    sum_part = zeros(size(X(:,1),1)*T,K(1),J);
    
    for j = 1:J
        Xmat = cell2mat(X(:,j));
        n = size(Xmat,1);
        % first get the predicted value
        for day = T:-1:1
            predict(day:T:n,j) = param_hat(j).bhaz(day)*...
                exp(Xmat(day:T:n,:)*param_hat(j).b);
        end
        % now calculate the sum of the predicted values (the integrated
        % hazard)
        sum_predict(:,j) = reshape(cumsum(reshape(predict(:,j),T,[]),1),size(predict(:,j)));
        % now calculate the negative log of the integrated hazard
        % (this is the same as -log(-log(S(t)))
        predict(:,j) = -log(sum_predict(:,j));
        % now get the multipliers and partials (for the standard errors)
        for day = T:-1:1
            param_index = mtx2param(zeros(size(param_ary)),J,arrayfun(@(x) length(x.bhaz),param_hat),K);
            param_index(j).b(:)=1;
            index = find(param2mtx(param_index));
            mult(day:T:n,j) = exp(Xmat(day:T:n,:)*param_hat(j).b);
            partial(day:T:n,:,j) = [Xmat(day:T:n,:)*param_hat(j).bhaz(day).*repmat(mult(day:T:n,j),1,length(index))];
        end
        % now get the cumulative sum of the partials
        sum_part(:,:,j) = reshape(cumsum(reshape(partial(:,:,j),T,[],K(j)),1),size(partial(:,:,j)));
        % now that we have the predicted values, multipliers, and partials,
        % calculate the A matrix and standard errors
        for day = T:-1:1
            param_index = mtx2param(zeros(size(param_ary)),J,arrayfun(@(x) length(x.bhaz),param_hat),K);
            param_index(j).bhaz(day)=1;
            param_index(j).b(:)=1;
            index = find(param2mtx(param_index));
            A = zeros(n/T,length(param_ary));
            A(:,index) = [mult(day:T:n,j) sum_part(day:T:n,:,j)]./(-repmat(sum_predict(day:T:end,j),1,length(index)));
            SE(day:T:n,j) = arrayfun(@(x) sqrt(A(x,:)*Hinv*A(x,:)'),1:n/T);
        end
    end
    fprintf('done.\n');
end
