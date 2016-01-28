function path = datapath()

path= mfilename('fullpath'); %path of the current file
path= fileparts(path); %remove filename to get the folder
path= fullfile(path,'Data');

addpath(genpath(fileparts(path)))

end