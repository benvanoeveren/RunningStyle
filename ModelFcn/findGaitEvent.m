function event = findGaitEvent(traj,traj_raw,varargin)
% Finds events from the kinetic model traj and the raw accelerations from
% the xsens sensors.
%
% 1. First remove the heading
% 2. Based on the found steps, it will search for the other events.
% IC_ initial contact, based on accelerations
% TO_ toe-off, based on gyroscope
% MST_ mid stance, lowest point of the COM
%
% Ben van Oeveren, 21-12-2015


%% inputs
if isfield(traj,'fs'); fs = traj.fs;
elseif any(ix); fs = varargin{ix+1};
else fs = 100; end

%% functions
segNr = @(x) find(strcmpi([traj.segment(:).name],x));

%% settings
SF = fliplr(60./[30 200].*fs); %duration of a step in samples

%% 1. Remove heading
[traj_rot, old_R_new]= rotate2subjHeading(traj); %x:ant-post, y:medio-late, z:vertical
raw_rot = transform_traj(traj_raw,old_R_new,[0 0 0]);

%% Use signals
%all signals are based on the foot
pos_R = traj_rot.segment(segNr('rightfoot')).joint(:,:);
pos_L = traj_rot.segment(segNr('leftfoot')).joint(:,:);
acc_R = raw_rot.segment(segNr('rightfoot')).sens_acc;
acc_L = raw_rot.segment(segNr('leftfoot')).sens_acc;
gyr_R = raw_rot.segment(segNr('rightfoot')).sens_ang_vel;
gyr_L = raw_rot.segment(segNr('leftfoot')).sens_ang_vel;

% COM
% if isfield(traj_rot,'CenterOfMass'); COM = traj_rot.CenterOfMass(:,1:3);
% else COM = traj_rot.segment(segNr('Pelvis')).origin(:,1:3); end
COM = traj_rot.segment(segNr('Pelvis')).com(:,1:3);

%% Get steps
%based on minimal speed in x+z.
[IC_R_temp,TO_R_temp] = getSteps(pos_R(:,4:6),fs,SF);
[IC_L_temp,TO_L_temp] = getSteps(pos_L(:,4:6),fs,SF);

%% Initial contact
%based on acceleration peaks, near initial IC
IC_R = initialContact(acc_R,IC_R_temp,fs,SF);
IC_L = initialContact(acc_L,IC_L_temp,fs,SF);

%% Toe-off
%toe-off based on the gyroscope, after IC.
TO_R = toeOff(acc_R,gyr_R,IC_R,TO_R_temp,fs,SF);
TO_L = toeOff(acc_L,gyr_L,IC_L,TO_L_temp,fs,SF);

%% midstance
% 50% of IC to TO
MST_R = round(TO_R+(TO_R-IC_R)./2);
MST_L = round(TO_L+(TO_L-IC_L)./2);
% MST_R = midstance(IC_R,TO_R,fs,pos_R,COM);
% MST_L = midstance(IC_L,TO_L,fs,pos_L,COM);

%% vertical oscillation
event.vertOsci_R = verticalOscillation(COM,IC_R,fs);
event.vertOsci_L = verticalOscillation(COM,IC_L,fs);

%% stepFreq
event.strideFreqR = 60*fs./diff(IC_R);
event.strideFreqL = 60*fs./diff(IC_L);

show = false;
if show
    % Plot results for right
    foot_R = pos_R(:,3);
    foot_L = pos_L(:,3);
    gyr_R = sqrt(sum(gyr_R.^2,2));
    gyr_R = lowpassfilterzerolag([1,10],fs,5,gyr_R);
    
    figure; hold on;
    plot(event.vertOsci_R);
    plot(foot_R); plot(foot_L);
    plot(zscore(acc_R(:,end))./10,'g')
    plot(zscore(gyr_R)./10,'r')
    vertlines(MST_R,'color','b')
    vertlines(IC_R,'color','g');
    vertlines(TO_R,'color','r');
    legend({'vertOsci','foot_R','foot_L','acc_R','gyr_R'},'interpreter','none')
    vertlines(TO_R_temp,'color','m')
    xlim([0 200])
    
    [h1,h2,h3] = plot_3d2(traj_rot,[],[],[],fs,[],[],COM(:,3));
    vertlines(IC_R,'color','g');
    vertlines(MST_R,'color','b');
    vertlines(TO_R,'color','r');
end

%% Phases
R = [IC_R,MST_R,TO_R];
L = [IC_L,MST_L,TO_L];

%% Stance times
event.ST_R = diff(R(:,[1,3]),[],2)./fs;
event.ST_L = diff(L(:,[1,3]),[],2)./fs;

%% flight times
event.FL_R = diff([R(1:end-1,3),R(2:end,1)],[],2)./fs;
event.FL_L = diff([L(1:end-1,3),L(2:end,1)],[],2)./fs;

%% dubble support
suppR = R(:,[1,3]);
suppL = L(:,[1,3]);
event.DS_R = doublesupport(suppL,suppR);
event.DS_L = doublesupport(suppR,suppL);

%% output
event.IC_R = IC_R;
event.TO_R = TO_R;
event.MST_R = MST_R;
event.IC_L = IC_L;
event.TO_L = TO_L;
event.MST_L = MST_L;

end

function [IC,TO] = getSteps(pos,fs,SF)
%% Filter signal
Vel = sum(diff(pos(:,[1,3]),[],1),2);
Vel = lowpassfilterzerolag(5,fs,5,Vel);
Vel(Vel<0)=0;
[TO,IC] = findintercept(Vel,...
    'method','down',...
    'Threshold',0.005,...
    'MinInterceptDistance',min(SF));
if diff(Vel(IC(1):IC(1)+1))>0 || ...
        ~(Vel(IC(1))<0 && Vel(IC(1)+1)>0);
    IC(1) = [];
end
if TO(1)<IC(1); TO(1) = []; end


end
function IC = initialContact(Raw,IC_temp,fs,SF)
%% get the ones with the minimal distance
Raw_v = sqrt(sum(Raw.^2,2));
Raw_lp = lowpassfilterzerolag([1,15],fs,5,Raw_v);
[~,newSteps1] = findintercept(Raw_lp);
IC_temp2 = findFirstAfter(IC_temp,newSteps1);

[~,newSteps2] = findpeaks(Raw_v);
IC = nan(size(IC_temp));
for i_st = 1:length(IC_temp)
    x = IC_temp2(i_st);
    a = newSteps2([find(newSteps2<=x,1,'last');find(newSteps2>x,1,'first')]);
    if ~isempty(a)
        [~,i_small] = max(Raw_v(a));
        IC(i_st) = a(i_small);
    end
end

end
function TO = toeOff(acc,gyr,IC,TO_temp,fs,SF)
% a valley in the gyroscope is corresponding to a toe-off
% the initial search range is limited from stride to stride
% the TO closest to TO_temp is used.

gyr_lp = lowpassfilterzerolag([1 5],fs,5,gyr());
[~,X] = findpeaks(-gyr_lp(:,2),'MinPeakDistance',min(SF),'MinPeakHeight',0);
rng = [0.1 1]*fs;

[index,NR] = findNrsBetwBound(X,[IC+rng(1) [IC(2:end);length(gyr)]]); %after IC and before p
TO = nan(size(IC));
for i_step =1 :length(index);
    data = NR{i_step};
    if ~isempty(data)
        [~,i]= min(abs(data-TO_temp(i_step)));
        TO(i_step) = data(i);
    end
end

end
function MST = midstance(IC,TO,fs,foot,COM)
% Defined as the lowest point of the pelvic (COM)
signal = lowpassfilterzerolag([1 5],fs,5,COM(:,3));
MST = nan(size(IC));
for i_step =1:length(IC)
    rng = IC(i_step):TO(i_step);
%     n = round(length(rng)*.18);
%     n(n<1)=1;
%     rng = rng(n:end-n);
    if numel(rng)<2 || isempty(rng)
        continue;
    end
    data = signal(rng);
    [~,i_max] = findpeaks(abs(data));
%     [~,i_max] = max(abs(data));
if numel(i_max)==1
    MST(i_step) = i_max+IC(i_step);
elseif ~isempty(i_max)
    [~,i_min] = min(abs(i_max - length(rng)/2));
    MST(i_step) = i_max(i_min)+IC(i_step);
%     try close(100); end
%     figure(100); 
%     plot(data)
%     vertlines(i_max)
end
end

%% alternative approach,
% when the com accelerates forward/upward
% In this case the angle of acceleration is between 0 and 90.
% x = signal(:,1); z = signal(:,3);
% ang = atan2d(x,z);
% ang0 = find(ang>0); ang90 = find(ang>90);
% MST = findFirstAfter(IC,[ang0;ang90]);
% figure; hold on; plot(ang)
% horzlines(0,'color','k');horzlines(90,'color','k')
% vertlines(IC,'color','g'); vertlines(TO,'color','r');


end
function DS = doublesupport(S1,S2)
indexR = findNrsBetwBound(S1(:,1),S2); %after IC and before p
ixR = ~cellfun(@isempty,indexR);
ix2x = cellfun(@(x) length(x)>1,indexR);
if any(ix2x); indexR(ix2x) = {indexR{ix2x}(end)}; end
DS = double(ixR); DS(~ixR) = nan;
DS(ixR) = S1([indexR{ixR}],1);
end
function VO = verticalOscillation(COM,IC_R,fs)
vertOsci = lowpassfilterzerolag([1 15],fs,5,COM(:,3));
extremes = @(x) [min(x),max(x)];
VO = nan(size(IC_R,1),2);
for i_step = 1:length(IC_R)-1
    try VO(i_step,:) = extremes(vertOsci(IC_R(i_step):IC_R(i_step+1))); end
end
end
%% subfunctions
function I = findLastBefore(p,X)
I = nan(size(p));
for i_st = 1:length(p)
    x_p = p(i_st);
    try I(i_st) = X(find(X<x_p,1,'last')); end
end
end
function I = findFirstAfter(steps,X)
I = nan(size(steps));
for i_st = 1:length(steps)
    xStep = steps(i_st);
    try I(i_st) = X(find(X>xStep,1,'first')); end
end
end
function smoothsignal=lowpassfilterzerolag(lpfreq,sampfreq,butorder,thesignal)
for i = 1: size(thesignal,2)
    [num,den]=butter(butorder,lpfreq./(sampfreq/2));%5 Hz lowpass filter 100 hz filter design, butterworth, 4th order
    smoothsignal(:,i) =(filtfilt(num,den,thesignal(:,i)));
end
end