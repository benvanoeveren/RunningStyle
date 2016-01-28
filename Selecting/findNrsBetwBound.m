function [index,NR] = findNrsBetwBound(Nrs,boundries)
% [V,Index] = FINDNRSBETWBOUND(Nrs,boundries) Finds numbers within boundries. 
% Nrs: vector of nrs
% boundries: double column array [start,finish]
% V = found nrs
% index = found indixes
% Ben van Oeveren, 13-03-2015
index = [];
NR = [];
if isempty(Nrs)|| isempty(boundries); return; end

%remove nans
nonan = find(~any(isnan(boundries),2));
if any(size(Nrs)==1) 
    Nrs=Nrs(:);
end

% index = [];
L = size(boundries,1);
index = cell(L,1); 
NR  = cell(L,1); 
for i_event = nonan';
  ind = find(any(Nrs>boundries(i_event,1),2) & any(Nrs<=boundries(i_event,2),2));
  V  = Nrs(ind,:);
  index{i_event} = ind;
  NR{i_event} = V;
%   if isempty(V)
%       warning('warnig')
%   end
end

end
