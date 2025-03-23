clear; clc; close all;
warning off backtrace
%% COCO compatible continuation problem of 2 DoF impact Duffing oscillator
% considering harmonic excitation and asymetric air gaps
% using nondimensionalized form

% Dependencies
addpath(genpath('_toolbox'))
addpath('pwsdde cont','plot tools');

% Initialize COCO
addpath('<COCO_dir>'); % location of COCO installation
startup;

% Default plot options
set(0, 'DefaultLineLineWidth', 1);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');


%% System definition
addpath('_system def/two dof impact oscillator');
nsnl_def;   % system definition file (returns 'struct')
sys = struct2func(struct);


%% Finding a periodic orbit

% Periodic orbit initialisation based on simulation
p0 = [1.5 1.0 1.0 0.05 0.9 0.02 0.127 0.952 0.1 0.1]; % parameter set
y0 = @(t) [0; 0; 0; 0; 1; 0]; % initial conditions
o_init.N = 2;                   % number of events
o_init.M = 100;                 % initial guess mesh dimension
sim_opts.t_end = 150*pi/p0(1);  % simulation time
sim_opts.h_act = [1 1 0 0];     % only normal events are active
[orb0,res]= pwsdde_sim(y0,p0,sys,o_init,sim_opts);
% figure(); plot_sim_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb1,'nsnl_p'); title('Simulated guess');


%% Follow the periodic orbit with COCO

% define a COCO compatible problem
p_names = {'$\omega$','$\delta_1$','$\delta_2$', '$\mu$', '$r$', ...
    '$\zeta$', '$\xi$', '$\phi$', '$\eta_1$', '$\eta_2$'};
[prob1, data1] = pwsdde_coco_prob(sys,orb0,p_names);

% set continuation options
prob1 = coco_set(prob1, 'cont', 'NAdapt', 0);
prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',5);
prob1 = coco_set(prob1,'cont','ItMX',[100 100]);

% bifurcation detection routines
prob1 = event_graze_bd(prob1,data1); % boundary grazing
prob1 = event_graze_int(prob1,data1); % interior grazing

% create a branch of solutions with COCO
bd1 = coco(prob1,'nsd_r1',[],1,{'$\omega$','$|\mu_c|$'},[1, 5]);


%% Initialize and follow a different orbit

% Initialization through simulation
p0(1) = 0.7; sim_opts.t_end = 100*pi/p0(1);
[orb1,res]= pwsdde_sim(y0,p0,sys,o_init,sim_opts);
% figure(); plot_sim_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb1,'nsnl_p'); title('Simulated guess');

% % define and run the continuation problem
p_names = {'$\omega$','$\delta_1$','$\delta_2$', '$\mu$', '$r$', ...
    '$\zeta$', '$\xi$', '$\phi$', '$\eta_1$', '$\eta_2$'};
[prob2, data2] = pwsdde_coco_prob(sys,orb1,p_names);
prob2 = coco_set(prob2,'cont','h_min',1e-3,'h_max',5,'ItMX',[100 100]);
prob2 = event_graze_bd(prob2,data2);
prob2 = event_graze_int(prob2,data2);
bd2 = coco(prob2,'nsd_r2',[],1,{'$\omega$','$|\mu_c|$'},[0.2, 2.5]);


%% Plot continuation results with COCO

% read continuation data
nsd_r1 = coco_bd_read('nsd_r1');
nsd_r2 = coco_bd_read('nsd_r2');

% bifurcation diagram
thm = struct();
thm.special = {'EP','BGR','IGR','SC'};
thm.BGR = {'dk','MarkerSize',7,'MarkerFaceColor','r'};
thm.IGR = {'ok','MarkerSize',7,'MarkerFaceColor','r'};
thm.SC = {'ob','MarkerSize',7};
thm.ustab = '$|\mu_c|$';
thm.ustabfun =  @(x) (x<1) + 1;
thm.lspec = {{'b--', 'LineWidth', 1}, {'b-', 'LineWidth', 1}};
figure(); box on; hold on
coco_plot_bd(thm,'nsd_r1','$\omega$','$|A_1|$')
coco_plot_bd(thm,'nsd_r2','$\omega$','$|A_1|$')
hold off

% display an orbit

% lab = coco_bd_labs('nsd_r1','BGR');
% orb_bgr = orb_get_data('nsd_r1',lab);
% figure(); plot_orb(orb_bgr,'nsnl_p')

% lab = coco_bd_labs('nsd_r1','IGR');
% orb_bgr = orb_get_data('nsd_r1',lab);
% figure(); plot_orb(orb_bgr,'nsnl_p')

% lab = coco_bd_labs('nsd_r2','BGR');
% orb_bgr = orb_get_data('nsd_r2',lab);
% figure(); plot_orb(orb_bgr,'nsnl_p')

lab = coco_bd_labs('nsd_r2','IGR');
orb_bgr = orb_get_data('nsd_r2',lab);
figure(); plot_orb(orb_bgr,'nsnl_p')

% lab = coco_bd_labs('nsd_r2','SC');
% orb_sn = orb_get_data('nsd_r2',lab(2));
% figure(); plot_spectrum(orb_sn.mu)

% animate orbit progression
% anim_bd_orb('nsd_r1','nsnl_p')