clear; clc; close all;
%which('=VU 3D model=')
sbj_weight = 67;

% filename = 'Ben-001';
fPath = '/Users/benvanoeveren/Documents/VU/Projecten/running style/Pilot 10-12-2015';
addpath(genpath(fPath))
 filename = 'Ben-001';
%filename = 'Bart_001-008';
%filename = 'subject6-008';
try 
    load(filename); 
catch
    mvnx = load_mvnx2(filename);
    [traj_mvn,traj_mvn_acc]= mvn2traj2(mvnx);
    save(filename,'traj_mvn','traj_mvn_acc')
end
traj = build_bodymodel(traj_mvn,sbj_weight);
nFrames = size(traj.segment(1).origin,1);
rng = 4000:7000;
traj_mvn_acc = trim_traj(traj_mvn_acc,rng);
traj = trim_traj(traj,rng);
% [traj]= trimDistance(traj,10,20);
event = findGaitEvent(traj,traj_mvn_acc);
% traj = trim_traj(s,10:20);

[combined_com,combined_mass] = calc_combined_com(traj,1:length(traj_mvn.segment));
fs = traj.fs;?

traj = rotate2subjHeading(traj);
%% XCOM
% COM = traj.centerOfMass;
% XCOM = calcXcom(COM,fs);
% traj.single.marker = COM; 
% traj.single(2).marker = combined_com; 

[h1,h2,h3] = plot_3d2(traj,[],[],[],fs);
