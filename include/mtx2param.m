function [param] = mtx2param(mtx,J,T,K)
    % index the baseline hazards (start and end indices for each exit
    % state)
    ibhaz0 = 1+[0 cumsum(T)]';
    ibhaz1 = circshift(ibhaz0,-1)-1;
    % index the regression coefficients (start and end indices for each
    % exit state)
    ib0 = ibhaz1(J)+1+[0 cumsum(K)]';
    ib1 = circshift(ib0,-1)-1;
    % construct an array of parameter structures (one struct for each exit
    % state)
    param = struct(...
        'bhaz',arrayfun(@(x) mtx(ibhaz0(x):ibhaz1(x)),1:J,...
            'UniformOutput',false),...
        'b',arrayfun(@(x) mtx(ib0(x):ib1(x)),1:J,...
            'UniformOutput',false) ...
    );
end
