function [mtx] = param2mtx(param)
    % stack the baseline hazard coefficients (in order of exit state) on
    % top of the regression parameters (in order of exit state)
    mtx = [vertcat(param.bhaz); vertcat(param.b)];
end
