function [UP,DOWN,varargout] = findintercept(Yin,varargin)
% [up,down] = findintercept(Yin), detects the intercept with signal Yin
% Options: 'Threshold','MinGradient','InterceptWidthRange',...
% 'MinInterceptDistance','MinInterceptAvarage','MinPeakHeight'
% followed by a number and Method: 'Linear','Nearest'
%
% 'Threshold': Gives the intercept with Yin
% 'MinSlope': The slope calculated from the first neighbour peak to
%  neighbour valley. Supports both single value and range.
% 'InterceptWidthRange': Minimal and maximal width between intercepting points.
%  Supports both single value and range. Starts with first up or down. To
%  force a start use 'StartPair'.
% 'MinInterceptDistance': Distance between intercept points of either downward
%  or upward intercepts.
% 'MinInterceptAvarage': Avarage over points above or under the intercept.
%  Supports both single value and two values for averages above and below threshold.
% 'StartPair': 'Up','Down','First'_ forces starting with an
%  upward or downward slope, by default it will start with the 'first'
%  'NIntercepts': specifies the maximum number of intercepts
%
%
% Optional output:
%   no output =  figure with results
%   varargout{1} = Av; %Average between intercept points
%   varargout{2} = slope; %Slope of intercepting point
%   varargout{3} = x; %x after interpolation
%   varargout{4} = Yin; %Y after interpolation
%
% Ben van Oeveren, 15-12-2016

% extract the parameters from the input argument list
[x,threshold,minSlope,intW,minD,peakAv,method,startPair,Nintc,minPeak] ...
    = parse_inputs(Yin,varargin{:});
Yin = Yin(:);

% Get all intercepts
[up,down] = getAllIntercepts(x,Yin,threshold);
[up,down,x,Yin] = methodForIntercept(x,up,down,Yin,threshold,method);
[up,down] = makePair(up,down,startPair);
[up,down] = interceptWidth(x,up,down,intW);
[up,down,Av,peakSum] = interceptAverage(x,Yin,up,down,peakAv,threshold,nargout,minD,minPeak);
[up,down,slope] = interceptSlope(x,Yin,up,down,minSlope,nargout,minD);
[up,down,Av,slope,peakSum] = minInterceptDist(x,peakSum,Av,slope,up,down,minD);
[up,down,Av,slope,peakSum] = keepAtMostNpPeaks(x,up,down,Nintc,Av,slope,peakSum);

UP = x(up);
DOWN = x(down);

if nargout>2
    varargout{1} = Av; %Average between intercept points
    varargout{2} = slope; %Slope of intercepting point
    varargout{3} = x; %when method is linear, than extra points are made
    varargout{4} = Yin; %when method is linear, than extra points are made
end

if nargout==0
    figure; hold on;
    grid on
    plot(x,Yin,'.-')
    line([xlim],[threshold threshold],'linestyle','--')
    plot(x(up),Yin(up),'or')
    plot(x(down),Yin(down),'og')
end

end %fcn end
function [x,threshold,minSlope,intW,minD,peakAv,method,startPair,Nintc,minPeak] = parse_inputs(Yin,varargin);

% Validate input signal
validateattributes(Yin,{'numeric'},{'nonempty','real','vector'},...
    'findpeaks','Y');
yIsRow = isrow(Yin);
y = Yin(:);

% copy over orientation of y to x.
xIsRow = yIsRow;

% indicate if the user specified an Fs or X
hasX = ~isempty(varargin) && isnumeric(varargin{1});

if hasX
    startArg = 2;
    if isscalar(varargin{1})
        % Fs
        Fs = varargin{1};
        validateattributes(Fs,{'double'},{'real','finite','positive'},'findpeaks','Fs');
        x = (0:numel(y)-1).'/Fs;
    else
        % X
        Xin = varargin{1};
        validateattributes(Xin,{'double'},{'real','finite','vector','increasing'},'findpeaks','X');
        if numel(Xin) ~= numel(Yin)
            if coder.target('MATLAB')
                throwAsCaller(MException(message('signal:findpeaks:mismatchYX')));
            else
                coder.internal.errorIf(true,'signal:findpeaks:mismatchYX');
            end
        end
        xIsRow = isrow(Xin);
        x = Xin(:);
    end
else
    startArg = 1;
    % unspecified, use index vector
    x = (1:numel(y)).';
end


if coder.target('MATLAB')
    p = inputParser;
    % Options: 'Threshold','MinGradient','MinInterceptWidth','MaxInterceptWidth',...
    % 'MinInterceptDistance','MinInterceptAvarage'
    % followed by a number and Method: 'Linear','Nearest'
    
    %put in defaults
    addParameter(p,'Threshold',0);
    addParameter(p,'MinSlope',[]);
    addParameter(p,'MinInterceptAvarage',[-inf inf]);
    addParameter(p,'InterceptWidthRange',[-inf,inf]);
    addParameter(p,'MinInterceptDistance',0);
    addParameter(p,'MinPeakHeight',[0 0]);
    addParameter(p,'Method','nearest');
    addParameter(p,'StartPair','First');
    addParameter(p,'NIntercepts',inf);
    parse(p,varargin{startArg:end}); %place varargin in struct p
else
    warning('no Matlab coder?')
end

%output
threshold = p.Results.Threshold;
minSlope = p.Results.MinSlope;
intW = p.Results.InterceptWidthRange;
minD = p.Results.MinInterceptDistance;
peakAv = p.Results.MinInterceptAvarage;
method = p.Results.Method;
startPair = p.Results.StartPair;
Nintc = p.Results.NIntercepts;
minPeak = p.Results.MinPeakHeight;

end
function [up,down] = getAllIntercepts(x,y,threshold)
t = y>threshold;
ud = diff(t);
up = find(ud ==1);
down = find(ud ==-1);
if isempty(up) && isempty(down) %there are no intercepting points
    return;
elseif isempty(down)
    down = [1;length(y)];
elseif isempty(up)
    up = [1;length(y)];
else
    %% include extremes
    if up(1)<down(1)
        if down(1)~=1; down = [1;down]; end
    else %start with down
        if up(1)~=1; up = [1;up]; end
    end
    
    L = length(y);
    if up(end)<down(end) %last is up, than make it a down
        if down(end)~=L; down = [down;L]; end
    else %start with down
        if up(end)~=L; up = [up;L]; end
    end
end
end
function [up,down,x,Yin] = methodForIntercept(x,up,down,Yin,threshold,method)
if ~any(up)||~any(down); return; end
if ~isempty(method)
    switch lower(method)
        case{'nearest'} %Get index closest to threshold
            up = nearest(up,Yin,threshold);
            down = nearest(down,Yin,threshold);
        case{'linear'}
            upx =interpolate(([up,up+1]),x,Yin,threshold);
            downx =interpolate(([down,down+1]),x,Yin,threshold);
            %new data points should be added to Yin and up and down
            [x,v] = sort([x;upx;downx]);
            newY = threshold*ones(length(up)+length(down),1);
            Yin2 = [Yin;newY]; %add newY
            Yin = Yin2(v); %new sorted Yin
            up = find(ismember(x,upx));
            down = find(ismember(x,downx));
        case{'first'}
            %             warning('sdf')
    end
end

    function x_out =interpolate(xx,x,y,threshold)
        y_0= [y(xx(:,1)),y(xx(:,2))];
        x_0 = [x(xx(:,1)),x(xx(:,2))];
        opp = diff(y_0,[],2);
        adj = diff(x_0,[],2);
        alp = atan(opp./adj);
        x_out = (threshold-y_0(:,1))./tan(alp)+x_0(:,1);
    end

    function x0 = nearest(x,Yin,threshold)
        x0 = x;
        if x(end) == length(Yin); x(end) = x(end)-1; end
        if numel(threshold)==1
            [~,i] = min(abs([Yin(x),Yin(x+1)]-threshold),[],2);
        else
            [~,i] = min(abs([Yin(x),Yin(x+1)]-[threshold(x) threshold(x)]),[],2);
        end
        x0=x+i-1;
    end

end
function [up,down] = makePair(up,down,startPair)
if ~any(up)||~any(down); return; end

switch lower(startPair)
    case{'first'}
        [L,~] = min([length(up),length(down)]);
        up = up(1:L);
        down = down(1:L);
    case{'up'}
        %Start with up and end with down
        start = find(down>up(1),1,'first');
        down = down(start:end);
        [L,~] = min([length(up),length(down)]);
        up = up(1:L);
        
        if length(down)~=length(up)
            warning('Length does not match')
        end
        if any((down-up)<0)
            warning('order does not match')
        end
    case{'down'}
        %Start with down and end with up
        start = find(up>down(1),1,'first');
        up = up(start:end);
        down = down(1:length(up));
        
        if length(up)~=length(down)
            warning('Length does not match')
        end
        if any((up-down)<0)
            warning('order does not match')
        end
end
ix = abs(up-down)<1;
up(ix)=[];
down(ix)=[];
end
function [up,down] = interceptWidth(x,up,down,intW)
if ~any(up)||~any(down); return; end

if numel(intW)<2
    intW = [intW, inf];
end
if all(isinf(intW)); return; end

if up(1)<down(1) %force starting with up
    %Start with up and end with down
    start = find(down>up(1),1,'first');
    down = down(start:end);
    up = up(1:length(down));
end

A = [x(up) x(down)];
A = diff(A,[],2);
ind = A>intW(1) & A<intW(2);
up = up(ind);
down = down(ind);

% copy peak values and locations to a temporary place
% firstDown = up(1)>down(1); %first is true, it starts with up.
% if firstDown; s1 = x(down); s2 = x(up); else s1 = x(down); s2 = x(up); end

% V = false(size(up));
% V(1) = true;
% keep =1;
% while keep
%     test = s2(keep:end)-s1(keep);
%     if all(test>intW(2)); keep = false; end
%     keep0 = find(test>intW(1),1,'first');
%     V(keep0+keep-1) = true;
%     keep = keep0+keep-1;
% end
% ind = find(V);

%
% if firstDown
%     down = down(ind);
%     up = up(ind+1);
% else
%     up = up(ind);
%     down = down(ind+1);
% end
end
function [up,down,aV,peakSum] = interceptAverage(x,Yin,up,down,peakAv,threshold,nOut,minD,minPeak)
if ~any(up)||~any(down); aV=[]; peakSum=[]; return; end
if all(isinf(peakAv)) && nOut<3 && isempty(minD); aV=[]; peakSum=[]; return;end
%Can be positive and negative!
% if numel(peakAv)<2
%     if peakAv-threshold<0; %down
%         peakAv = [threshold,peakAv];
%     else %peakAv-threshold>0
%         peakAv = [peakAv,threshold];
%     end
% end
% if numel(minPeak)<2
%     minPeak = [minPeak,-minPeak];
% end

Yin= Yin-threshold;
t = sort([up;down]); %indices of up and down

L = size(t,1);
ampl=nan(L,1);
mpeak=nan(L,2);
for i = 1:L-1
    rng = t(i):t(i+1);
    time = x(rng);
    dtime = diff(time([1,end]));
    if dtime<=0
        ampl(i) = 0;
        mpeak(i,:) = [0 0];
    else
        ampl(i) = trapz(time, Yin(rng)) ./ diff(time([1,end]));
        mpeak(i,:) = [max(Yin(rng)),min(Yin(rng))];
    end
end
if up(1)<down(1);
    upA = ampl(1:2:end);
    downA = ampl(2:2:end);
    
    upP = mpeak(1:2:end,1);
    downP = mpeak(2:2:end,2);
else %start with down
    upA = ampl(2:2:end);
    downA = ampl(1:2:end);
    
    upP = mpeak(2:2:end,1);
    downP = mpeak(1:2:end,2);
end

%% remove the next indices
if numel(peakAv)==2 && all(~isinf(peakAv))
    index1 = any([upA<peakAv(1),downA>=peakAv(2)],2);
elseif all(~isinf(peakAv))
    index1 = upA<peakAv(1);
else
    index1 = false(size(upA));
end

if numel(minPeak)==2
    index2 = any([upP<minPeak(1),downP>=minPeak(2)],2);
else
    index2 = upA<minPeak(1);
end
index = index1|index2;
%remove indexes below peakAv
up(index)=[];
down(index)=[];
aV = [upA(~index),downA(~index)];

if size(aV,2)==1; aV = [aV aV]; end

dist = abs([x(up)-x(down)]);
peakSum = aV.*[dist dist];

end
function [up,down,slope] = interceptSlope(x,Yin,up,down,minSlope,nOut,minD)
if ~any(up)||~any(down); slope = []; return; end
if isempty(minSlope) && nOut<4 && isempty(minD); slope = []; return; end
% pair = [up,down];
% based on the first local min and max before or after the intercept
% with down: %peak before,intercept,valley after
% with up: %valley before,intercept,peak after
[valley_y,valley_x] = findpeaks(-Yin,x);
[peaks_y,peaks_x] = findpeaks(Yin,x);

if numel(minSlope)==2 || isempty(minSlope)
    minSlope = sort(minSlope,'descend');
    if isempty(minSlope); minSlope = [-inf inf]; end
    option = 'slopeUpDown'; %gradUp,gradDown
elseif minSlope>0;
    option = 'slopeUp'; %gradUp,gradDown
else %minSlope<0;
    option = 'slopeDown'; %gradUp,gradDown
    minSlope = [inf,minSlope];
end

slope = nan(length(up),2);
switch option
    case{'slopeUp','slopeUpDown'}
        u =x(up);
        if ~isempty(u);
            for i = 1:length(u)
                n1 = find(valley_x<u(i),1,'last');
                n2 = find(peaks_x>u(i),1,'first');
                
                if ~isempty(n1) && ~isempty(n2);
                    changeY = diff([valley_y(n1),peaks_y(n2)]);
                    changeX = diff([valley_x(n1),peaks_x(n2)]);
                    slope(i,1) = changeY/changeX;
                end
                
            end
        end
end
switch option
    case{'slopeDown','slopeUpDown'}
        d =x(down);
        if ~isempty(d)
            for i = 1:length(d)
                n1 = find(peaks_x<d(i),1,'last');
                n2 = find(valley_x>d(i),1,'first');
                if ~isempty(n1) && ~isempty(n2);
                    changeY = diff([peaks_y(n1),valley_y(n2)]);
                    changeX = diff([peaks_x(n1),valley_x(n2)]);
                    slope(i,2) = changeY/changeX;
                end
            end
        end
end
index = any([slope(:,1)<minSlope(1),slope(:,2)>minSlope(2)],2);
up(index) = [];
down(index) = [];

end
function [up,down,aV,slope,peakSum] = minInterceptDist(x,peakSum,aV,slope,up,down,minD)
if isempty(peakSum); return; end
idx = findPeaksSeparatedByMoreThanMinPeakDistance(x,peakSum,[up down],minD);
up = up(idx);
down = down(idx);
try aV = aV(idx,:); end
try slope = slope(idx,:); end
try peakSum = peakSum(idx,:); end

    function idx = findPeaksSeparatedByMoreThanMinPeakDistance(x,y,iPk,minD)
        % Start with the larger peaks to make sure we don't accidentally keep a
        % small peak and remove a large peak in its neighborhood.
        
        if isempty(iPk) || minD==0
            idx = [1:size(iPk,1)]';
            return
        end
        
        % copy peak values and locations to a temporary place
        pks = abs(y); %use max?
        [~,col] = max(nansum(pks,1)); pks = pks(:,col);
        locs = x(iPk(:,col));
        
        % Order peaks from large to small
        [~, sortIdx] = sort(pks,'descend');
        locs_temp = locs(sortIdx);
        
        idelete = ones(size(locs_temp))<0;
        for i = 1:length(locs_temp)
            if ~idelete(i)
                % If the peak is not in the neighborhood of a larger peak, find
                % secondary peaks to eliminate.
                idelete = idelete | (locs_temp>=locs_temp(i)-minD)&(locs_temp<=locs_temp(i)+minD);
                idelete(i) = 0; % Keep current peak
            end
        end
        
        % report back indices in consecutive order
        idx = sort(sortIdx(~idelete));
    end
end
function [up,down,Av,slope,peakSum] = keepAtMostNpPeaks(x,up,down,Np,Av,slope,peakSum)
%based on the maximal surface
if isinf(Np);return; end
S = peakSum(:,1);
[~,sinx] = sort(S,'descend');
if length(up)>Np
    sinx =sort(sinx(1:Np));
    up = up(sinx);
    down = down(sinx);
    Av =  Av(sinx);
    slope = slope(sinx);
    peakSum = peakSum(sinx);
end

end
