function x = xlsreadsettings(fullfilename,tabname,headerrow,varargin)
% read metadata from xlsfile
% input (1): filepath
% input(2): cell with tabnames
% output: struct containing the fields according to the headers of the
% metadata. First two columns contain the Metadata[label,data].
% Ben van Oeveren, benvoeveren@gmail.com
% 06-05-2015: Handles empty cells and date formats better.

if isempty(headerrow);
    headerrow = 1;
end
% varargin, pairs with label and number of columns e.g. {'Condition',2}

if ~iscell(tabname)
     error('Tabname should be a cell')
end


for j_sheet = 1:length(tabname)
    [num,txt,both] = xlsread(fullfilename,tabname{j_sheet});
    %%metadata
    header = txt(headerrow,:);
    %first two columns contain the metadata
    %     %% Metadata from both
    %     for i_row = headerrow+1:size(txt,1)
    %         value = both{i_row,2};
    %         field = deblank(both{i_row,1});
    %         x.(field) = value;   %labels +value
    %     end
    
    %% Built struct from varargin
    x.empty = [];
    for i_var = 1:2:length(varargin);
        columnNr= find(ismember(header,varargin{i_var}));
        index = varargin{i_var+1};
        indexdata = both(headerrow+1:end,columnNr+index);
        try timingbothi = ~any(cellfun(@(x) ...
                isnan(x)|strcmpi(x,'#VALUE!')|strcmpi(x,'#REF!'),indexdata),2);
        catch
            timingbothi = ones(length(indexdata),1);
        end
        
        if all(cellfun(@isnumeric,indexdata));
            indexdata = cell2mat(indexdata(timingbothi,:));
        end
        
        condition = both(headerrow+find(timingbothi),columnNr);
        
        %remove nans
        c_nan = cellfun(@(x) all(isnan(x)),condition,'uni',1);condition(c_nan)={''};
        c_word = cellfun(@(x) ~any(regexp(x,'[a-zA-Z0-9]')),condition,'uni',1);
        condition(c_word)={''};
        %replace space for _
        condition = cellfun(@(x) regexprep(x,' ','_'),condition,'uni',0);
        
        
        for i_row=1:length(condition)
            value = indexdata(i_row,:);
            if iscell(value); value = cat(1,value); end
            if isnan(condition{i_row}); continue; end
            
            time_i= ~cellfun(@isempty,regexpi(condition(i_row),'.*time|.*tijd')); 
 
            field = deblank(condition{i_row});
            field = regexp(field,'\w*','match','once');
            
            if time_i
                value = datetime(value,'ConvertFrom','excel');
                [~, ~, ~, H, MN, S] = datevec(value);
                value = H*3600+MN*60+S;
            end
            
            if ~isnan(field);
                if ~isfield(x,field)
                    x.(field) = value;
                else %add to existing field
                     x.(field) =  [x.(field); value];
                end
            end
        end
                
        fields = fieldnames(x);
        datei= ~cellfun(@isempty,regexpi(fields,'.*date|.*datum|.*day|.*dag'));
        for i_date = find(datei)'
            v = x.(fields{i_date});
            if iscell(v); v = v{1}; end
            if isdatetime(v);
                x.(fields{i_date}) = v;
            elseif ~isnan(v)
                if isnumeric(v);
                    value = forcedateformat(x2mdate(v));
                    %value = datetime(v,'ConvertFrom','excel');
                else
                    value = forcedateformat(v,'dd-mm-yyyy');
                end
            x.(fields{i_date}) = value;
            end
        end 
    end
    
end




end

function dateout = forcedateformat(datein,varargin)
% Force date format string based on the current date.
%
%recognizes: 
%    strings: '19-Aug-2014' or '19-08-2014' or '19/08/2014'
%    vectors: (now)/ 7.3583e+05
if isnan(datein)
    return;
end

A = regexp(datein,'\d*|\w*','match');

%asume ymd,mdy,dym,dmy,myd,ydm   
allperms = perms(1:3); %all possible permutations
%1 =y, 2=m, 3 = d

%% Check month or day string
% put together as date
if length(A)<3 %likely only yy and month
    A= [A,0];
    remove_size = allperms(:,3)~=3;
else
    remove_size = zeros(size(allperms,1),1);    
end
%1 =y, 2=m, 3 = d
nums = cellfun(@(x) str2double(x),A);
remove_string = allperms(:,isnan(nums))==1; %string can be a day or a month (3 or 2), not 1.
if isempty(remove_string);remove_string=zeros(size(remove_size));end
remove_year1 = ~any(allperms(:,nums>9)==1,2);
remove_year = allperms(:,nums>30)~=1; %only year can be larger than 30
if isempty(remove_year);remove_year=zeros(size(remove_size));end
remove_month = any(allperms(:,nums>12)==2,2);

try remove_tot = remove_size|remove_string|remove_year1|remove_year|remove_month;
catch % invalid date construction
    dateout = [];
    return
end
if all(remove_tot) % invalid date construction
    dateout = [];
    return
end
allperms(remove_tot,:)=[];

C = A(allperms); %make strings
for i_row = 1:size(C,1)
    x = strjoin(C(i_row,:),'-');
    try N{i_row} = etime(datevec(now),datevec(x,'yy-mm-dd'));
    catch 
        N{i_row} = nan;
    end
end
[~,ind] = nanmin(abs(cell2mat(N)));

x = strjoin(C(ind,:),'-');
dateout = datenum(x,'yy-mm-dd');


%when the format is given
if ~isempty(varargin)
    formatOut = varargin{1};
    dateout = datestr(dateout,formatOut);
end

end


