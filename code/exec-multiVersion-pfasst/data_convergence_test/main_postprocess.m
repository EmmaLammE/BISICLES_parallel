clear;close all;clc;

% filename='pf_nlevel2_nproc1_stream.L1L2.0128.petsc.000000.2d.hdf5';
% hinf=h5info(filename);
% 
% % quick disp of data
% h5disp(filename)
% 
% datasets_names=hinf.Datasets.Name;
% group_names=hinf.Groups.Name;


%% data grid size 64
dt_ref_pf=[0.1,0.02,0.01,0.002,0.001];

% ref(t=0.1)-pf(t=0.1) etc
maxnorm_ref_pf=[6.5728,6.9339,6.9775,7.0122,7.0166];

% ref(t=0.1)-ref(t=0.001),   ref(t=0.02)-ref(t=0.001) 
% ref(t=0.01)-ref(t=0.001),  ref(t=0.002)-ref(t=0.001)
dt_ref_ref=[0.1,0.02,0.01,0.002];
L2_norm_ref_ref=[9.6482e3,1.8063e3,8.5301e2,9.455e1];


% pf(t=0.1)-pf(t=0.001),   pf(t=0.02)-pf(t=0.001) 
% pf(t=0.01)-pf(t=0.001),  pf(t=0.002)-pf(t=0.001)
dt_pf_pf=dt_ref_ref;
L2_norm_pf_pf=[1.6494e-5,2.2372e-7,2.2281e-7,2.5486e-7];

%% data grid size 128
clear;clc
dt_pf_pf       =[0.14,     0.12,     0.1,     0.08,     0.06,...
                 0.05,     0.04,     0.02,     0.01,     0.002];
% L2_norm_pf_pf  =[3.6215e-4,2.1259e-4,8.7030e-5,3.0458e-5,1.0413e-5,...
%                  5.4361e-6,2.3664e-7,1.9514e-7,1.9882e-7,2.2120e-7]; % not computing velocity
L2_norm_pf_pf  =[1.7683e-1,1.0261e-1,4.0954e-2,1.4159e-2,4.8407e-3,...
                 2.5286e-3,1.0431e-3,7.2203e-5,1.5273e-5,1.9871e-5]; % computing velocity

dt_ref_ref     =[0.14,     0.12,     0.1,     0.08,     0.06,...
                 0.05,     0.04,     0.02,     0.01,     0.002,    0.001];
L2_norm_ref_ref=[1.3830e+4,1.2017e+4,9.5206e+3,7.2573e+3,5.5089e+3,...
                 4.6649e+3,3.7149e+3,1.8350e+3,9.0506e+2,1.6577e+2,7.3650e+1];

pf_slope=4;
dt_pf_pf_slope=dt_pf_pf(1:end-2);dt_pf_pf_slope=[0.16,dt_pf_pf_slope];
L2_norm_pf_pf_slope=exp(pf_slope*log(dt_pf_pf_slope)+...
    1*(log(L2_norm_pf_pf(1))-pf_slope*log(dt_pf_pf_slope(1))));
dt_pf_pf_slope=dt_pf_pf_slope+0.01;

ref_slope=1;
dt_ref_ref_slope=dt_ref_ref(1:end);
L2_norm_ref_ref_slope=exp(ref_slope*log(dt_ref_ref_slope)+...
    1.05*(log(L2_norm_ref_ref(1))-ref_slope*log(dt_ref_ref_slope(1))));

dt_pf_ref     =[0.14,     0.12,     0.1,     0.08,     0.06,...
                 0.05,     0.04,     0.02,     0.01,     0.002,    0.001];
L2_norm_pf_ref=[1.3821e+4,1.2013e+4,9.5260e+3,7.2763e+3,5.5459e+3,...
                 4.7154e+3,3.7877e+3,2.0141e+3,1.2567e+3,9.1875e+2,9.1053e+2];




%% self convergence plot
close all;

%%% plot convergence loglog grid 128
left_color=  [ 64,  125,  82 ]/255;
right_color=[128,0,32 ]/255;
middle_color=[26 ,  85, 153 ]/255;

f=figure(1);
figpos=[.1 .1 .8 .8];
set(f,'defaultAxesColorOrder',[left_color; right_color]);

yyaxis left
loglog(dt_ref_ref,L2_norm_ref_ref,'LineWidth',1.5,'LineStyle','-','Marker','o','color',left_color)
hold on
% loglog(dt_pf_ref,L2_norm_pf_ref,'LineWidth',1.5,'LineStyle','-','Marker','s','color',middle_color)
loglog(dt_ref_ref_slope,L2_norm_ref_ref_slope,'LineWidth',1.5,'LineStyle','--','color',left_color)
hold off
ylabel('Serial L2 norm, reference to serial dt=1e-3');
ylim([1e0,1e5]);

yyaxis right
loglog(dt_pf_pf,L2_norm_pf_pf,'LineWidth',1.5,'LineStyle','-','Marker','^','color', right_color)
hold on
loglog(dt_pf_pf_slope,L2_norm_pf_pf_slope,'LineWidth',1.5,'LineStyle','--','color',right_color)
hold off
ylim([1e-5,1e0]);xlim([5e-4,3e-1]);
grid on
xlabel('dt, year');
ylabel('PFASST L2 norm, reference to PFASST dt=1e-3');
hold off

legend({'serial convergence, ref soln: serial soln at dt=2e-4','',...
    'parallel convergence, ref soln: parallel soln at dt=1e-3',''},'FontSize',14);

f.Position=[100,100,600,600];

figure_name=['./figs/fig_self_temporal_converg_pfpf_and_refref_grid128.png'];
saveas(gcf,figure_name);

%% plot ref pf L2 norm on loglog grid 128
close all;

%%% plot convergence loglog grid 128
left_color=  [ 64,  125,  82 ]/255;
right_color=[128,0,32 ]/255;
middle_color=[26 ,  85, 153 ]/255;

f=figure(2);
figpos=[.1 .1 .8 .8];
loglog(dt_ref_ref,L2_norm_ref_ref,'LineWidth',1.5,'LineStyle','-','Marker','o','color',left_color)
hold on
loglog(dt_pf_ref,L2_norm_pf_ref,'LineWidth',1.5,'LineStyle','-','Marker','s','color',middle_color)
loglog(dt_ref_ref_slope,L2_norm_ref_ref_slope,'LineWidth',1.5,'LineStyle','--','color',left_color)
hold off
ylabel('L2 norm');
ylim([1e1,1e5]);xlim([5e-4,3e-1]);
grid on
xlabel('dt, year');

legend({'serial convergence, ref soln: serial soln at dt=2e-4',...
    'serial-parallel convergence, ref soln: serial soln at each dt'},'FontSize',14);

f.Position=[100,100,600,600];

figure_name=['./figs/fig_temporal_converg_pfref_grid128.png'];
saveas(gcf,figure_name);

%% plot ref pf difference semi log, grid 128
close all;

left_color=  [26 ,  85, 153 ]/255;
right_color=[128,0,32 ]/255;

f=figure(2);
figpos=[.1 .1 .8 .8];
set(f,'defaultAxesColorOrder',[left_color; right_color]);

semilogx(dt_ref_pf,ref_point_val,'LineWidth',1.5,'LineStyle','-','Marker','o','color', left_color);hold on
semilogx(dt_ref_pf,pf_point_val,'LineWidth',1.5,'LineStyle','--','Marker','^','color', left_color);hold off


xlabel('dt, year');
ylabel('Thickness at ~(79400, 116800)');


% loglog(dt_ref_pf,maxnorm_ref_pf,'LineWidth',1.5,'LineStyle','-','Marker','x')
% xlabel('dt, year');
% ylabel('max norm, between Serial and PFASST');
% figure_name=['./figs/fig_self_temporal_converg_pfref.png'];
% saveas(gcf,figure_name);

figure_name=['./figs/fig_converge_to_same_field?.png'];
saveas(gcf,figure_name);
