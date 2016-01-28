function [s] = load_mvnx2(varargin)
% mvnx = load_mvnx(filename)
% load mvnx file
% filename is name of file
% mvnx is struct containing all data of mvnx file
% Ben van Oeveren, 7 dec. 2015
% Use mvn2traj2(s) to load the file into the VU 3d model.

%% force filename
filename = varargin{1};
[~,f,~] = fileparts(filename);
filename = fullfile([f,'.mvnx']); %force mvnx name
filename = which(filename);

%% 2. Open as text file.
R = importfile1(filename);
s = getfields(R,[]);

end


function s = getfields(r,s)
fields = cellfun(@(x) regexpi(x,'<(?!/)(\w*)','tokens','once'),r,'uni',0);
ind2 = ~cellfun(@isempty,fields);
fields = fields(ind2);
R3 = r(ind2);
[uniqfields,~,g] = unique([fields{:}],'stable');

for i_field = unique(g)'
    temp_ind = i_field ==g;
    fieldname = uniqfields{i_field};
    data = R3(temp_ind,:);
    
    %% check first to find indices of data
    data1 = data{1,:};
    if strcmpi(fieldname,'points');
        frame_start = findFirstinCell(r, '<frames');
        n = regexp(r(1:frame_start),'points');
        ind = find(~cellfun(@isempty,n));
        n = cell2mat(n(ind));
        L = ind(n==3) - ind(n==2);
        s.segmentlabelnum = (L-1)/3; %van point to point kost 3 regels
    end
    %check string
    type = checkFirst(data1);
    
    switch type
        case 1
            new = forStrings(data);
        case 2
            [~,a,b] = regexp(data1,'>(.*)<','match');
            finish = length(data1)-b(end);
            matrix = cell(size(data));
            for i_row = 1:length(data)
                d = data{i_row};
                matrix{i_row,:} =  d(a+1:end-finish-1);
            end
            new = matrix;
        case 3
            new = forData(data);
    end
    s.(fieldname) = new;
end
end
function new = forStrings(data)
data1 = data{1,:};
subfields = regexp(data1,'(\w*)="','tokens');
subfields = [subfields{:}];

for i_subfields = 1:length(subfields)
    matrix = cell(size(data));
    subfield = subfields{i_subfields};
    for i_datarow = 1:size(data,1)
        d = data{i_datarow};
        [~,x,z] = regexp(d,'="(.*?)"','tokens');
        rng = x(i_subfields)+2:z(i_subfields)-1;
        matrix{i_datarow,:} = d(rng);
    end
    if numel(matrix)==1; matrix = char(matrix{1}); end
    new.(subfield) = matrix;
end
end
function matrix = forData(data)
data1 = data{1,:};
data2 = regexp(data1,'\d*[.]{0,1}\d*','match');
[~,a,b] = regexp(data1,'>(.*)<','match');
finish = length(data1)-b(end);

i = [a(1),finish];
sz = [size(data,1),length(data2)];

if isempty(data); matrix = []; return; end
i = i+[1 -1];
matrix = nan(sz);
for i_row = 1:sz(1)
    d = data{i_row,:};
    d = d(i(1):end-i(2));
    matrix(i_row,:) = sscanf(d, '%f');
end
end
function type = checkFirst(data1);
if ~isempty(strfind(data1,'='))
    type=1;
elseif ~isempty(regexp(data1,'>(.*)<','match'))
    A = regexp(data1,'>(.*)<','match');
    if length(A{1})<=2;
        type = 0; return;
    end
    B = isstrprop(A{:}, 'alpha');
    if ~any(B)
        type=3; %this is data
    else
        type=2; %this is another string
    end
else
    type =0;
end
end
function i =findFirstinCell(R,expr)
found =false;
i = 1;
while ~found && i<length(R)
    found = ~isempty(strfind(R{i}, expr));
    if ~found
         i = i+1;
    end
end
end


function [dataArray] = importfile1(filename, startRow, endRow)
if nargin<=2;  startRow = 2;   endRow = inf; end

delimiter = '';
formatSpec = '%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
dataArray = dataArray{1,:};
%% Close the text file.
fclose(fileID);
end