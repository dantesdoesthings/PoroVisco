
function [files] = fdir(directory, filetype)
%FILLFULLDIR Summary of this function goes here
%   Detailed explanation goes here

% clear
% directory = '/Users/DanStrange/Dropbox/Eng Work/Research/Instron/Gelatin15/Gelatin15-4.is_crelax_RawData';
% filetype =  'Spec*.csv';

%% check whether / is included at the end. If it is delete it.

if directory(length(directory)) == '/'
    directory(length(directory)) = '';
end

%% Create cell of files in the directory

filecell = struct2cell(dir([directory '/' filetype])); %get cell of file info
filenames = filecell(1,:); % cell file names
files = {};
for i = length(filenames):-1:1 
    files{i} = [directory '/' filenames{i}];  % append directory onto name
end


%% Determine if there are folders in the directory

celldir  = struct2cell(dir(directory)); % get listing of full directory as a cell
findex = cell2mat (celldir(4, :));      % finds the indexs of cells which are folders
folders =  celldir(1, findex);          % extracts these names into another cell

%% Remove spurious folders eg. '.' and '..'

for i = 1:size(folders,2)
    if (strcmp(folders{i},'.')||strcmp(folders{i},'..'))
        folders{i} = '';
    end
end

folders = folders(~cellfun('isempty', folders));

%% search folders for files and append to file list using fDir function
for i = 1:size(folders,2)
    folderpath = [directory '/' folders{i}];
    folderfiles = fdir(folderpath, filetype);
    files = [files folderfiles];
end


end

