function [traj_mvn,Height,CP_sample] = calc_trunk_COMs(traj_mvn,CP_sample)
%height and mass of every COM

if ~exist('CP_sample','var'); CP_sample =47; end
maxRow = size(traj_mvn.segment(1).blm,1); %added by Ben
if CP_sample>maxRow
    CP_sample = maxRow;
end

RASIS =traj_mvn.segment(1).blm(CP_sample,7*3-2:7*3);
LASIS =traj_mvn.segment(1).blm(CP_sample,8*3-2:8*3);
pelvic_width         =norm(RASIS-LASIS);

Height.MASIS  =(RASIS(3)+LASIS(3))/2;
Height.MPSIS  =traj_mvn.segment(1).blm(CP_sample,5*3);
Height.L5S1   =(Height.MASIS+Height.MPSIS)/2+0.111790*pelvic_width;
Height.Navel  =Height.MASIS;
Height.Xiph   =traj_mvn.segment(5).blm(CP_sample,6*3);
Height.C7     =traj_mvn.segment(5).blm(CP_sample,5*3);
Height.HV     =traj_mvn.segment(7).blm(CP_sample,2*3);
Height.Hips   =Height.L5S1-0.30*pelvic_width;
Height.COM_pe =Height.Hips+(Height.Navel-Height.Hips)*0.3885;
Height.COM_ab =Height.Navel+(Height.Xiph-Height.Navel)*0.5498;
Height.COM_th =Height.Xiph +(Height.C7-Height.Xiph)*0.4934;
Height.COM_he =Height.C7   +(Height.HV-Height.C7)*0.4998;

Height.L5S1_COM_pe =Height.COM_pe-Height.L5S1;
Height.L5S1_COM_ab =Height.COM_ab-Height.L5S1;
Height.L5S1_COM_th =Height.COM_th-Height.L5S1;
Height.L5S1_COM_he =Height.COM_he-Height.L5S1;

% figure(1); hold on;
% fields = fieldnames(Height);
% 
% for i_field = 1:length(fields)
%     fdata = Height.(fields{i_field});
%     plot(0,fdata,'*')
% end
% 
% figure(2); hold on;
% for i_field = 1:length(traj_mvn.segment)
%     fdata = traj_mvn.segment(i_field).origin(1,:);
%     plot3(fdata(:,1),fdata(:,2),fdata(:,3),'*')
% end
% 
% traj_mvn2 = trim_traj(traj_mvn,1:10); 
% plot_3d(traj_mvn2)

end

