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




