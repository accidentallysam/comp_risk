function [] = output_goodnessoffit(X,param_hat,cov_labels,fval,exitflag,output,s,H,predict,pSE,predict_S,pSE_S,exit_state,spell_len)
    % store some basic dimensionality information
    J = length(param_hat);
    T = max(arrayfun(@(x) length(x.bhaz),param_hat));
    N = length(spell_len);

%     S1 = [1; exp(-exp(-mean(reshape(predict_S(:,1),T,[]),2)))];
%     S2 = [1; exp(-exp(-mean(reshape(predict_S(:,2),T,[]),2)))];
%     h1 = [0; mean(reshape(predict(:,1),T,[]),2)];
%     h2 = [0; mean(reshape(predict(:,2),T,[]),2)];
    
    S1 = exp(-exp(-predict_S(:,1)));
    S2 = exp(-exp(-predict_S(:,2)));
    h1 = predict(:,1);
    h2 = predict(:,2);
    
    S1 = [ones(1, size(S1,1)/T); reshape(S1,T,[])];
    S2 = [ones(1, size(S2,1)/T); reshape(S2,T,[])];
    h1 = [zeros(1, size(h1,1)/T); reshape(h1,T,[])];
    h2 = [zeros(1, size(h2,1)/T); reshape(h2,T,[])];

    P1hat = zeros(T,1);
    P2hat = zeros(T,1);
    P1 = zeros(T,1);
    P2 = zeros(T,1);
    
    %mean(reshape(XXX,T,[]),2);
    
    for t = 2:T+1
        P1hat(t-1) = mean(h1(t,:)./(h1(t,:)+h2(t,:)).*(S1(t-1,:).*S2(t-1,:)-S1(t,:).*S2(t,:)));
        P2hat(t-1) = mean(h2(t,:)./(h1(t,:)+h2(t,:)).*(S1(t-1,:).*S2(t-1,:)-S1(t,:).*S2(t,:)));
        P1(t-1) = sum(exit_state==1 & spell_len==(t-1))/N;
        P2(t-1) = sum(exit_state==2 & spell_len==(t-1))/N;
    end
    output_newline();
    output_divline();
    fprintf('t  \tP1hat\tP1\tP2hat\tP2\tR_hat\tRatio\n');
    output_divline();
    fprintf('%02d\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n',[(1:T)',P1hat,P1,P2hat,P2,P2hat./P1hat,P2./P1]');
    output_newline();
    fprintf('Total\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n',[sum(P1hat(1:3,:)),sum(P1(1:3,:)),sum(P2hat(1:3,:)),sum(P2(1:3,:)),sum(P2hat(1:3,:))./sum(P1hat(1:3,:)),sum(P2(1:3,:))./sum(P1(1:3,:))]');
    output_newline();

    h1 = h1*1.3;
    S1 = exp(1.3*log(S1));
    for t = 2:T+1
        P1hat(t-1) = mean(h1(t,:)./(h1(t,:)+h2(t,:)).*(S1(t-1,:).*S2(t-1,:)-S1(t,:).*S2(t,:)));
        P2hat(t-1) = mean(h2(t,:)./(h1(t,:)+h2(t,:)).*(S1(t-1,:).*S2(t-1,:)-S1(t,:).*S2(t,:)));
        P1(t-1) = sum(exit_state==1 & spell_len==(t-1))/N;
        P2(t-1) = sum(exit_state==2 & spell_len==(t-1))/N;
    end
    output_newline();
    output_divline();
    fprintf('t  \tP1hat\tP1\tP2hat\tP2\tR_hat\tRatio\n');
    output_divline();
    fprintf('%02d\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n',[(1:T)',P1hat,P1,P2hat,P2,P2hat./P1hat,P2./P1]');
    output_newline();
    fprintf('Total\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n',[sum(P1hat(1:3,:)),sum(P1(1:3,:)),sum(P2hat(1:3,:)),sum(P2(1:3,:)),sum(P2hat(1:3,:))./sum(P1hat(1:3,:)),sum(P2(1:3,:))./sum(P1(1:3,:))]');
    output_newline();
end
