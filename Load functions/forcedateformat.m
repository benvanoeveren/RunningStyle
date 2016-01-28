function dateout = forcedateformat(datein,varargin)
% Force date format string based on the current date.
%
%recognizes: 
%    strings: '19-Aug-2014' or '19-08-2014' or '19/08/2014'
%    vectors: (now)/ 7.3583e+05
if isnan(datein)
    return;
end

if isnumeric(datein)
    dateout = datestr(datein,'yy-mm-dd');
    return
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