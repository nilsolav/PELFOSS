function [pos] = surveycn_transectendpoints(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [pos.STRATUM,TRANSECT,LAT,LON,vessel] = IMPORTFILE(FILENAME) Reads data
%   from text file FILENAME for the default selection.
%
%   [STRATUM,TRANSECT,LAT,LON,vessel] = IMPORTFILE(FILENAME, STARTROW,
%   ENDROW) Reads data from rows STARTROW through ENDROW of text file
%   FILENAME.
%
% Example:
%   [stratum,transect,lat,lon,vessel] = importfile('transect_endpoints_IESNS2017.txt',2, 121);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/04/26 09:19:46

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
% pos.stratum = dataArray{:, 1};
% pos.transect = dataArray{:, 2};
% pos.lat = dataArray{:, 3};
% pos.lon = dataArray{:, 4};
% pos.vessel = dataArray{:, 5};

for i=1:(size(dataArray{1},1)/2)
    pos(i).stratum  = dataArray{1}(i*2);
    pos(i).transect = dataArray{2}(i*2);
    pos(i).vessel  = dataArray{5}{i*2};
    % lat
    pos(i).latstart = dataArray{3}(i*2-1);
    pos(i).latstop = dataArray{3}(i*2);
    % lon
    pos(i).lonstart = dataArray{4}(i*2-1);
    pos(i).lonstop = dataArray{4}(i*2);
end


