function [traj]= trimDistance(traj,distance,startingPoint,varargin)
% [traj]= trimDistance() trims the traj based on a given distance and starting point

%% parse inputs
ix = cellfun(@(x) strcmpi(x,'fs'),varargin);
if isfield(traj,'fs'); fs = traj.fs;
elseif any(ix); fs = varargin{ix+1};
else fs = 100; end
if ~exist('startingPoint','var') || isempty(startingPoint); startingPoint =1; end

% Use data
if isfield(traj,'centerOfMass'); pos = traj.centerOfMass(startingPoint:end,:);
else pos = traj.segment(1).origin(startingPoint:end,1:3); end

x = pos(:,1); y = pos(:,2);
travelled = cumtrapz(x, y)./fs;
if ~travelled(end)>distance; 
    finish = find(travelled>distance,1,'first')+startingPoint;
else
    finish = length(distance)+startingPoint;
end

traj = trim_traj(traj_mvn,startingPoint:7000);
end