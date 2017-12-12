%%%
% EXECUTE AND OUTPUT ALL EMPIRICS
%

% set the current directory
cd '/Volumes/SVBROWN 1/Hazard Model/MATLAB'
% add the "include" folder to the path
addpath '/Volumes/SVBROWN 1/Hazard Model/MATLAB/include'

output_heading('HAZARD MODEL');

% load  the data
[exit_state spell_len X cov_labels] = load_data('data/data.csv');

% create the data
% param = struct(...
%     'bhaz',{[0.04 0.03 0.02 0.01]', [0.02 0.03 0.04 0.05]'},...
%     'b',{[.1 .05]',[.01]'}...
%     );
% [exit_state spell_len X] = create_data(2000,param);

% estimate the model
[param_hat,fval,exitflag,output,s,H] = estimate_model(exit_state, spell_len, X);

% get the predicted values
[predict pSE] = predicted_values(X_better,param_hat,H);

close all

% output the results
output_results(X_better,param_hat,cov_labels,fval,exitflag,output,s,H,predict,pSE);

% now get the predicted negative log of the survivor function
[predict_S pSE_S] = predicted_S(X_better,param_hat,H);

% output the survivor plots
output_survivor(X_better,param_hat,cov_labels,fval,exitflag,output,s,H,predict_S,pSE_S);

% output goodness of fit
output_goodnessoffit(X_better,param_hat,cov_labels,fval,exitflag,output,s,H,predict,pSE,predict_S,pSE_S,exit_state,spell_len);
