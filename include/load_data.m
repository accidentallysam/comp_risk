function [exit_state spell_len X cov_labels] = load_data(fname)
    fprintf('loading data...');
    newData1 = importdata(fname);

    % Break the data up into a new structure with one field per column.
    colheaders = genvarname(newData1.colheaders);
    for i = 1:length(colheaders)
        dataByColumn1.(colheaders{i}) = newData1.data(:, i);
    end
    % store the variables
    vars = fieldnames(dataByColumn1);
    % first retrieve 'id', 'spell_num', 'spell_len', 'exit_state', and 't'
    %from the file (the first five columns of the file)
    for i = 1:5
        eval([vars{i} ' = dataByColumn1.(vars{i});']);
    end
    % second retrieve the matrix of covariates for each exit state (the
    % remaining columns of the file)
    J = 0; % number of exit states
    while 1
        xindex = cellfun(@(x) strcmp(x(end),num2str(J+1)),...
            newData1.colheaders);
        xindex(1:5) = zeros(1,5);
        if sum(xindex)==0 break; end;
        J = J+1;
        % get the covariates
        covars{J} = newData1.data(:,xindex);
        cov_labels(:,J) = vars(xindex);
    end % end while
    % remove exit state numbers from the end of the covariate labels
    cov_labels = cellfun(@(x) x(1:end-1),cov_labels,'UniformOutput',false);
    
    % index unique entries in the order that they appear in the file
    [C,index,ic] = unique([id spell_num],'rows');
    [index,ia,ic] = unique(index);

    nspell = length(index); % total number of spells

    % store the spell length and exit state for each spell
    exit_state = exit_state(index);
    spell_len = spell_len(index);

    % store the time-varying covariates for each spell (and each exit
    % state)
    for i = 1:nspell
        for j = 1:J
            % index the covarites for this individual spell
            cindex = id == id(index(i)) & spell_num == spell_num(index(i));
            % get the covariates
            x = covars{j}(cindex,:);
            % order covariets by time and store them
            X{i,j}(t(cindex),:) = x;
        end % end for
    end % end for
    fprintf('done.\n');
end % end load_data()

