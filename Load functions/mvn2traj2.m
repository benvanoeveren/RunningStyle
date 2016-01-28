function [traj, traj_mvn_acc] = mvn2traj2(s)
segment_names={
    'Pelvis'%1
    'L5' %2
    'L3' %3
    'T12' %4
    'T8' %5
    'Neck'%6
    'Head' %7
    'RightShoulder'% 8
    'RightUpperArm' %9
    'RightForeArm'
    'RightHand'
    'LeftShoulder'
    'LeftUpperArm'
    'LeftForeArm'
    'LeftHand'
    'RightUpperLeg'
    'RightLowerLeg'
    'RightFoot'
    'RightToe'
    'LeftUpperLeg'
    'LeftLowerLeg'
    'LeftFoot'
    'LeftToe'};

%% sizes
[rows, cols] = size(s.position);

%% create dummy trajectory
traj=create_dummy_traj(length(segment_names),rows,0);
[~,i_sens_list] = ismember(segment_names,s.sensor.label);

%% fill dummy
%% first add rows to fields with less samples
%normal = rows, otherwise rows-2
sfields = fieldnames(s);
[z,~] = structfun(@size,s);
changefields = sfields(z == rows-2); %search for fields with lower number of samples
for i_field = 1:length(changefields)
    data = s.(changefields{i_field});
    s.(changefields{i_field}) = [zeros(2,size(data,2));data];
end

%% start filling from struct s
A= reshape([1:cols],3,cols/3);
for i_s = 1:length(segment_names)
    traj.segment(i_s).name=segment_names(i_s);
%     traj_mvn.segment(i_s).position = s.position(:,A(:,i_s));
%     traj_mvn.segment(i_s).acceleration = s.acceleration(:,A(:,i_s));
    traj.segment(i_s).origin = s.position(:,A(:,i_s));
    traj.segment(i_s).com   = traj.segment(i_s).origin*nan;
    traj.segment(i_s).bracemarkers = traj.segment(i_s).bracemarkers*nan;
    traj.segment(i_s).joint = traj.segment(i_s).origin*nan;
    traj.segment(i_s).blm = traj.segment(i_s).blm*nan;
    
    q=s.orientation(:,i_s*4-3:i_s*4);
    traj.segment(i_s).R=quaternion2orientation_matrix(q);
    
    i_sens=i_sens_list(i_s);
    if i_sens~=0
        traj_mvn_acc.segment(i_s).name=segment_names(i_s);
        q = s.sensorOrientation(:,i_sens*4-3:i_sens*4);
        traj_mvn_acc.segment(i_s).sens_acc =  s.sensorAcceleration(:,A(:,i_sens));
        traj_mvn_acc.segment(i_s).sens_ang_vel =  s.sensorAngularVelocity(:,A(:,i_sens));
        traj_mvn_acc.segment(i_s).Rsensor = quaternion2orientation_matrix(q);
        traj_mvn_acc.segment(i_s).origin_pos = s.position(:,A(:,i_sens));
        traj_mvn_acc.segment(i_s).origin_acc = s.acceleration(:,A(:,i_sens));
        traj_mvn_acc.segment(i_s).ang_vel = s.angularVelocity(:,A(:,i_sens));
        traj_mvn_acc.segment(i_s).ang_acc = s.angularAcceleration(:,A(:,i_sens));
    end
end

%% put blm into struct
p = [0;cumsum(s.segmentlabelnum)];
for i_s=1:length(segment_names);
    gRseg = traj.segment(i_s).R;
    origin = traj.segment(i_s).origin;
    inx = p(i_s)+1:p(i_s+1);
    
    for i_point = 1:length(inx)
        gP = s.pos_s(inx(i_point),:);
        if ismember(i_s,[9:11 13:15]);
            segRg = transpose_col(gRseg(1,:));
            segP = prod_col_vect(segRg,gP);
        else
            segP=gP;
        end
        gP=origin+prod_col_vect(gRseg, segP);
        traj.segment(i_s).blm(:,i_point*3-2:i_point*3)= gP;
        traj.segment(i_s).blm_name = s.point.label(inx);
    end
end

subjFields= fieldnames(s.subject);
for i_field =1:length(subjFields)
    stringData = s.subject.(subjFields{i_field});
    B = isstrprop(stringData, 'alpha');
    if all(~B)
        s.subject.(subjFields{i_field}) = str2double(stringData);
    end
end
% S = s.subject;

%% additional fields
traj.fs= s.subject.frameRate;
traj.recDate = s.subject.recDate;
traj.originalFilename = s.subject.originalFilename;
traj.centerOfMass = s.centerOfMass;

return
%% pelvis
traj_pelvis=create_dummy_traj(2,length(frame_range),0);
i_s=1;
for i_frame=frame_range;
    if ismember(i_frame,[1:2, frame_range(end)]);
        traj_pelvis.segment(1).R(i_frame,:) = [nan nan nan nan nan nan nan nan nan];
        traj_pelvis.segment(2).R(i_frame,:) = [nan nan nan nan nan nan nan nan nan];
    else
        traj_pelvis.segment(1).name = 'pelvis segment';
        q = s.orientation(:,i_s*4-3:i_s*4);
        traj_pelvis.segment(1).R(i_frame,:) = quaternion2orientation_matrix(q);
        
        traj_pelvis.segment(2).name = 'pelvis sensor';
        q = s.sensorOrientation(:,i_s*4-3:i_s*4);
        traj_pelvis.segment(2).R(i_frame,:) = quaternion2orientation_matrix(q);
    end
end
end

function R=quaternion2orientation_matrix(q)
q0=q(:,1);
q1=q(:,2);
q2=q(:,3);
q3=q(:,4);

R=[1-2*q2.^2-2*q3.^2  2*q1.*q2+2*q3.*q0   2*q1.*q3-2*q2.*q0...
    2*q1.*q2-2*q3.*q0    1-2*q1.^2-2*q3.^2    2*q2.*q3+2*q1.*q0...           
    2*q1.*q3+2*q2.*q0   2*q2.*q3-2*q1.*q0  1-2*q1.^2-2*q2.^2];% see Xsens manual and http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/index.htm
end

%clear tree
