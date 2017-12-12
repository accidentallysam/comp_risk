clc

stepsize = 50;

fileToRead1 = '/Volumes/SVBROWN 1/Hazard Model/MATLAB/data/gps.csv';

% Import the file
newData1 = importdata(fileToRead1);

% Break the data up into a new structure with one field per column.
colheaders = genvarname(newData1.colheaders);
for i = 1:length(colheaders)
    dataByColumn1.(colheaders{i}) = newData1.data(:, i);
end

% Create new variables in the base workspace from those fields.
vars = fieldnames(dataByColumn1);
for i = 1:length(vars)
    assignin('base', vars{i}, dataByColumn1.(vars{i}));
end

clearvars -except stepsize CO latitude longitude

n = length(CO);
i = 1;

for k = 0:stepsize:n-1
    str = sprintf('http://maps.googleapis.com/maps/api/distancematrix/xml?origins=');
    index = (k+1:min(k+stepsize,n))';
    m = length(index);
    for j = 1:m
        str = strcat(str,sprintf('%.10f,%.10f',latitude(index(j)),longitude(index(j))));
        if j<m
            str = strcat(str,sprintf('|'));
        end
    end
    str = strcat(str,sprintf('&destinations=12.68041,-7.975924&mode=walking&units=metric&sensor=false'));
    urlwrite(str,sprintf('gps_xml%03d.xml',i));
    i = i+1;
    pause(5);
    % fprintf('http://maps.googleapis.com/maps/api/distancematrix/xml?origins=12.69127,-7.970567|12.6797,-7.97465&destinations=12.68041,-7.975924&mode=walking&units=metric&sensor=false');
end
