function [predict SE] = predicted_values(X,param_hat,H)
    fprintf('generating predicted values and standard errors...');
    J = size(X,2);
    T = max(arrayfun(@(x) length(x.bhaz),param_hat));
    K = arrayfun(@(x) length(x.b),param_hat);

    param_ary = param2mtx(param_hat);
    Hinv = inv(H);
    
    for j = 1:J
        Xmat = cell2mat(X(:,j));
        n = size(Xmat,1);
        for day = T:-1:1
            % first get the predicted value
            predict(day:T:n,j) = param_hat(j).bhaz(day)*...
                exp(Xmat(day:T:n,:)*param_hat(j).b);
            % second get its standard error
            param_index = mtx2param(zeros(size(param_ary)),J,arrayfun(@(x) length(x.bhaz),param_hat),K);
            param_index(j).bhaz(day)=1;
            param_index(j).b(:)=1;
            index = find(param2mtx(param_index));
            mult = exp(Xmat(day:T:n,:)*param_hat(j).b);
            A = zeros(n/T,length(param_ary));
            A(:,index) = [mult Xmat(day:T:n,:)*param_hat(j).bhaz(day).*repmat(mult,1,length(index)-1)];
            SE(day:T:n,j) = arrayfun(@(x) sqrt(A(x,:)*Hinv*A(x,:)'),1:n/T);
        end
    end
    fprintf('done.\n');
end
