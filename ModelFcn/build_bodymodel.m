function traj_mvn = build_bodymodel(traj_mvn,sbj_weight)

traj_mvn.segment(1).joint=[traj_mvn.segment(16).origin traj_mvn.segment(2).origin traj_mvn.segment(20).origin];%pelvis
traj_mvn.segment(2).joint=[traj_mvn.segment(2).origin traj_mvn.segment(3).origin ];% 2'L5' 
traj_mvn.segment(3).joint=[traj_mvn.segment(3).origin traj_mvn.segment(4).origin ];% 3'L3' 
traj_mvn.segment(4).joint=[traj_mvn.segment(4).origin traj_mvn.segment(5).origin ];% 4'T12' 
traj_mvn.segment(5).joint=[traj_mvn.segment(5).origin traj_mvn.segment(6).origin ];% 5'T8' 
traj_mvn.segment(6).joint=[traj_mvn.segment(6).origin traj_mvn.segment(7).origin ];% 6'Neck' 
traj_mvn.segment(7).joint=[traj_mvn.segment(7).origin traj_mvn.segment(7).blm(:,4:6)];% 7'Head'

traj_mvn.segment(8).joint=[traj_mvn.segment(8).origin traj_mvn.segment(9).origin];% 8'Right Shoulder' 
traj_mvn.segment(9).joint=[traj_mvn.segment(9).origin traj_mvn.segment(10).origin];% 9'Right Upper' 
traj_mvn.segment(10).joint=[traj_mvn.segment(10).origin traj_mvn.segment(11).origin];% 10'Right Forearm' 
traj_mvn.segment(11).joint=[traj_mvn.segment(11).origin traj_mvn.segment(11).blm(:,4:6)];% 11'Right Hand' 

traj_mvn.segment(12).joint=[traj_mvn.segment(12).origin traj_mvn.segment(13).origin];% 12'Left Shoulder' 
traj_mvn.segment(13).joint=[traj_mvn.segment(13).origin traj_mvn.segment(14).origin];% 13'Left Upper' 
traj_mvn.segment(14).joint=[traj_mvn.segment(14).origin traj_mvn.segment(15).origin];% 14'Left Forearm' 
traj_mvn.segment(15).joint=[traj_mvn.segment(15).origin traj_mvn.segment(15).blm(:,4:6)];% 15'Left Hand' 

traj_mvn.segment(16).joint=[traj_mvn.segment(16).origin traj_mvn.segment(17).origin];% 16'Right Upper Leg' 
traj_mvn.segment(17).joint=[traj_mvn.segment(17).origin traj_mvn.segment(18).origin];% 17'Right Lower Leg' 
traj_mvn.segment(18).joint=[traj_mvn.segment(18).origin traj_mvn.segment(19).origin];% 18'Right Foot' 
traj_mvn.segment(19).joint=[traj_mvn.segment(19).origin ];% 19'Right Toe' 

traj_mvn.segment(20).joint=[traj_mvn.segment(20).origin traj_mvn.segment(21).origin];% 20'Left Upper Leg' 
traj_mvn.segment(21).joint=[traj_mvn.segment(21).origin traj_mvn.segment(22).origin];% 21'Left Lower Leg' 
traj_mvn.segment(22).joint=[traj_mvn.segment(22).origin traj_mvn.segment(23).origin];% 22'Left Foot' 
traj_mvn.segment(23).joint=[traj_mvn.segment(23).origin ];% 23'Left Toe'

%plot_3d(traj_mvn)

RSJC=traj_mvn.segment(9).joint(:,1:3);
REJC=traj_mvn.segment(10).joint(:,1:3);
RWJC=traj_mvn.segment(11).joint(:,1:3);
RFT =traj_mvn.segment(11).joint(:,4:6);

LSJC=traj_mvn.segment(13).joint(:,1:3);
LEJC=traj_mvn.segment(14).joint(:,1:3);
LWJC=traj_mvn.segment(15).joint(:,1:3);
LFT =traj_mvn.segment(15).joint(:,4:6);

RHJC=traj_mvn.segment(16).joint(:,1:3);
RKJC=traj_mvn.segment(17).joint(:,1:3);
RAJC=traj_mvn.segment(18).joint(:,1:3);

LHJC=traj_mvn.segment(20).joint(:,1:3);
LKJC=traj_mvn.segment(21).joint(:,1:3);
LAJC=traj_mvn.segment(22).joint(:,1:3);

%trunk joints
HV=traj_mvn.segment(7).joint(:,4:6);
HeadC1=traj_mvn.segment(6).joint(:,4:6);
T1C7=traj_mvn.segment(5).joint(:,4:6);
T9T8=traj_mvn.segment(4).joint(:,4:6);
L1T12=traj_mvn.segment(3).joint(:,4:6);
L4L3=traj_mvn.segment(2).joint(:,4:6);
L5S1=traj_mvn.segment(1).joint(:,4:6);
MHJC=(RHJC+LHJC)/2;
MSJC=(RSJC+LSJC)/2;

Trunk_length=HV(1,3)-MHJC(1,3);%Zatsiorski = 85cm for subject 174 long
Trunk_length/HV(1,3)*100;%85/174*100= 49%
%figure;plot([HV(1,3) HeadC1(1,3) T1C7(1,3) T9T8(1,3) L1T12(1,3) L4L3(1,3) L5S1(1,3) MHJC(1,3)])

CP_sample= 1;
[traj_mvn,Height,CP_sample] = calc_trunk_COMs(traj_mvn,CP_sample); %height and mass

gPelvis_COP  = L5S1(CP_sample,:)+[0 0 Height.L5S1_COM_pe];
gAbdomen_COP = L5S1(CP_sample,:)+[0 0 Height.L5S1_COM_ab];
gThorax_COP  = L5S1(CP_sample,:)+[0 0 Height.L5S1_COM_th];
gHead_COP    = L5S1(CP_sample,:)+[0 0 Height.L5S1_COM_he];

% if i_sub==13
% gPelvis_COP  =MHJC(1,:)+[0 0 .06];
% gAbdomen_COP =MHJC(1,:)+[0 0 .22];
% gThorax_COP  =MHJC(1,:)+[0 0 .44];
% gHead_COP    =MHJC(1,:)+[0 0 .69];
% end


%gP=origin+prod_col_vect(gRseg, segP);

lPelvis_COP   =prod_col_vect(transpose_col(traj_mvn.segment(1).R(CP_sample,:)),(gPelvis_COP  -traj_mvn.segment(1).joint(CP_sample,4:6)));
lAbdomen_COP  =prod_col_vect(transpose_col(traj_mvn.segment(3).R(CP_sample,:)),(gAbdomen_COP -traj_mvn.segment(3).joint(CP_sample,1:3)));
lThorax_COP   =prod_col_vect(transpose_col(traj_mvn.segment(5).R(CP_sample,:)),(gThorax_COP  -traj_mvn.segment(5).joint(CP_sample,1:3)));
lHead_COP     =prod_col_vect(transpose_col(traj_mvn.segment(7).R(CP_sample,:)),(gHead_COP    -traj_mvn.segment(7).joint(CP_sample,1:3)));

traj_mvn.segment(1).com=traj_mvn.segment(1).joint(:,4:6)+prod_col_vect(traj_mvn.segment(1).R, lPelvis_COP);
traj_mvn.segment(3).com=traj_mvn.segment(3).joint(:,1:3)+prod_col_vect(traj_mvn.segment(3).R, lAbdomen_COP);
traj_mvn.segment(5).com=traj_mvn.segment(5).joint(:,1:3)+prod_col_vect(traj_mvn.segment(5).R, lThorax_COP);
traj_mvn.segment(7).com=traj_mvn.segment(7).joint(:,1:3)+prod_col_vect(traj_mvn.segment(7).R, lHead_COP);


%BLMs
C7=traj_mvn.segment(5).blm(:,13:15);
PX=traj_mvn.segment(5).blm(:,16:18);
PX(1,:)-C7(1,:);% is shorter than in de Leva
(L5S1(1,3)-MHJC(1,3)) / (MSJC(1,3)-MHJC(1,3));

%upper arm
traj_mvn.segment(9).com=RSJC+(REJC-RSJC)*0.5772;  
traj_mvn.segment(13).com=LSJC+(LEJC-LSJC)*0.5772; 
%lower arm
traj_mvn.segment(10).com=REJC+(RWJC-REJC)*0.4574;  
traj_mvn.segment(14).com=LEJC+(LWJC-LEJC)*0.4574; 
%lower arm
traj_mvn.segment(11).com=RWJC+(RFT-RWJC)*0.4574;  
traj_mvn.segment(15).com=LWJC+(LFT-LWJC)*0.4574; 

%upper leg
traj_mvn.segment(16).com= RHJC + (RKJC-RHJC)*0.3612;
traj_mvn.segment(20).com=LHJC + (LKJC-LHJC)*0.3612; 
%lower leg
traj_mvn.segment(17).com=RKJC + (RAJC-RKJC)*0.4395;  
traj_mvn.segment(21).com=LKJC + (LAJC-LKJC)*0.4395; 
%feet
traj_mvn.segment(18).com=traj_mvn.segment(18).blm(:,end-2:end);  
traj_mvn.segment(22).com=traj_mvn.segment(22).blm(:,end-2:end); 

%% copy the com to the origin
for i_seg=1:23;traj_mvn.segment(i_seg).origin=traj_mvn.segment(i_seg).com;end

% plot_3d(traj_mvn,[],[],[],120)%,[-.5 2 -2 .5 0 2 ],[])
traj_mvn.segment(1).mass   =sbj_weight/9.81* 11.18/100;%pelvis
traj_mvn.segment(3).mass   =sbj_weight/9.81* 16.33/100;%abdomen
traj_mvn.segment(5).mass   =sbj_weight/9.81* 15.96/100;%thotrax
traj_mvn.segment(7).mass   =sbj_weight/9.81* 6.94/100;%head
traj_mvn.segment(9).mass   =sbj_weight/9.81* 2.7/100;%RightUpperArm
traj_mvn.segment(10).mass  =sbj_weight/9.81* 1.63/100;%RightForeArm
traj_mvn.segment(11).mass  =sbj_weight/9.81* 0.61/100;%RightHand
traj_mvn.segment(13).mass  =sbj_weight/9.81* 2.7/100;%LeftUpperArm
traj_mvn.segment(14).mass  =sbj_weight/9.81* 1.63/100;%LeftForeArm
traj_mvn.segment(15).mass  =sbj_weight/9.81* 0.61/100;%LeftHand

traj_mvn.segment(16).mass   =sbj_weight/9.81* 14.17/100;% 'Right Upper Leg'
traj_mvn.segment(17).mass   =sbj_weight/9.81* 4.33/100;%    'Right Lower Leg'
traj_mvn.segment(18).mass   =sbj_weight/9.81* 1.37/100; %   'Right Foot'
traj_mvn.segment(20).mass   =sbj_weight/9.81* 14.17/100; %   'Left Upper Leg'
traj_mvn.segment(21).mass   =sbj_weight/9.81* 4.33/100; %   'Left Lower Leg'
traj_mvn.segment(22).mass   =sbj_weight/9.81* 1.37 /100;%  'Left Foot'
traj_mvn.segment([1 3 5 7 9 10 11 13 14 15 16:18 20:22]).mass;
