alldist = [];
for i = 1:14
    xDoc = xmlread(sprintf('gps_xml%03d.xml',i));
    allDist = xDoc.getElementsByTagName('distance');

    n = allDist.getLength;

    dist = zeros(n,1);

    for k = 1:n
        thisDist = allDist.item(k-1);
        dist(k) = str2num(thisDist.getElementsByTagName('value').item(0).getFirstChild.getData);
    end
    alldist(end+1:end+n,1) = dist;
end

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

clearvars -except alldist CO

csvwrite('data/gps_walk.csv',[CO alldist]);