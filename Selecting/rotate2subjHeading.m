function [traj_new old_R_new]= rotate2subjHeading(traj)
warning('rotate2subjHeading kan beter')

segNr = @(x) find(strcmpi([traj.segment(:).name],x));
fs = traj.fs;

if ~exist('signal','var')
    % use Pelvis position as default
    %     segm = {'Pelvis';'L5';'L3';'T12';'T8';'Neck';'Head'};
    segm = {'Pelvis'};
    segmNr = cellfun(@(x) segNr(x),segm);
    poss = [traj.segment(segmNr).com];
    heading = [nanmedian(poss(:,1:3:end),2), nanmedian(poss(:,2:3:end),2), nanmedian(poss(:,3:3:end),2)];
    head_POS_f = lowpassfilterzerolag(5,fs,5,heading);
    
    %     figure; hold on;
    %     hold on; plot3(heading(:,1),heading(:,2),heading(:,3))
    %     figure; hold on;
    %     plot(heading(:,3));
    %     sm = smooth(heading(:,3),sf*3,'loess');
    %     plot(sm)
    
%     head_POS_f(:,1) = smooth(heading(:,1),fs*3,'loess');
%     head_POS_f(:,2) = smooth(heading(:,2),fs*3,'loess');
%     head_POS_f(:,3) = smooth(heading(:,3),fs*3,'loess');
    
    if size(heading,1)>100
        segm = {'RightFoot','LeftFoot'};
        segmNr = cellfun(@(x) segNr(x),segm);
        posFeet = [traj.segment(segmNr).com];
        [~,~,Rsteps] = getSteps(posFeet(:,1:3),fs);
        [~,~,Lsteps] = getSteps(posFeet(:,4:6),fs);
        steps = sort([Rsteps;Lsteps]);
    else
        head_POS_f = heading; steps = [];
    end
    
    if isempty(steps) || length(steps)<2;
        %         warning('No steps detected to determine heading')
        xsteps =[];
    elseif length(steps)<4 && length(steps)>1
        xsteps = steps(1:2:end);
        meanPos = head_POS_f(xsteps,:);
    else
        steps = steps(1:floor(length(steps)/2)*2);
        xsteps1 = steps(1:2:end);
        xsteps2 = steps(2:2:end);
        xsteps = round(mean([xsteps1 xsteps2],2));

        meanPos = [head_POS_f(steps(1:2:end),:)+head_POS_f(steps(2:2:end),:)]./2;
    end
    %     figure; hold on;
    %     plot3(heading(:,1),heading(:,2),heading(:,3),'r')
    %     hold on; plot3(head_POS_f(:,1),head_POS_f(:,2),head_POS_f(:,3))
    %     hold on; plot3(head_POS_f(xsteps,1),head_POS_f(xsteps,2),head_POS_f(xsteps,3),'*-')
    
    xnew = 1:length(head_POS_f);
    if ~isempty(xsteps)
        x_axis_temp = interp1(xsteps,meanPos,xnew,'spline','extrap');
        %         figure; hold on; plot(head_POS_f-x_axis_temp);
    else
        x_axis_temp = head_POS_f;
    end
end
% build coordinate system
old_R_new = coordinateSystem(x_axis_temp);

% rotate trajectory
traj_new = transform_traj(traj,old_R_new,[0 0 0]);
%     set_ax = [85 130 -3 1 -2 2];
%     plot_3d2(traj_new,[],[],[],100,[],[],traj_new.segment(1).origin)  ;hline([0 10],'k:');
end

function old_R_new = coordinateSystem(x_axis_temp)
% build coordinate system
z_axis = repmat([0 0 1],size(x_axis_temp,1),1);
y_axis = cross(z_axis,x_axis_temp);
x_axis = cross(y_axis,z_axis);

x_axis   = x_axis./norm_col(x_axis);
y_axis   = y_axis./norm_col(y_axis);
z_axis   = z_axis./norm_col(z_axis);
old_R_new = [x_axis y_axis z_axis];
end

function [IC,TO,mstance] = getSteps(pos,fs)
SF = fliplr(60./[30 200].*fs); %duration of a step in samples
pos = lowpassfilterzerolag(15,fs,5,pos);
if mean(pos(:,1))<0;
    pos(:,1) = -pos(:,1);
end
if mean(pos(:,2))<0;
    pos(:,2) = -pos(:,2);
end

%% Filter signal
Vel = sum(pos(:,[1,3]),2);
Vel = lowpassfilterzerolag([1 5],fs,5,Vel);
% Vel(Vel<0)=0;
[TO,IC] = findintercept(Vel,...
    'method','down',...
    'MinInterceptDistance',min(SF));
if diff(Vel(IC(1):IC(1)+1))>0 || ...
        ~(Vel(IC(1))<0 && Vel(IC(1)+1)>0);
    IC(1) = [];
end
%start with IC
if TO(1)<IC(1); 
    TO(1:find(IC(1)>TO,1,'last'))=[]; 
end
n= min([length(IC),length(TO)]);
IC = IC(1:n);TO = TO(1:n);

mstance = round(mean([TO,IC],2));
end