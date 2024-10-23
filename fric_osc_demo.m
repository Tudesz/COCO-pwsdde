clear; clc; close all;
warning off backtrace
%% COCO compatible continuation problem of a 1 DoF dry friction oscillator
% using nondimensionalized form of governing equation of motion
% artificial DDE formalism, phase shift cosidered as time delay

% Dependencies
addpath(genpath('_toolbox'))
addpath('pwsdde cont','plot tools');
addpath('C:\Users\Tudo\Documents\MATLAB\COCO');
startup;

% Default plot options
set(0, 'DefaultLineLineWidth', 1);
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextInterpreter','latex');
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter','latex');


%% System definition
addpath('_system def/dry friction oscillator');
sys.f = 'fricosc_f';        % vector field modes
sys.e = 'fricosc_e';        % event functions
sys.tau = 'fricosc_tau';    % time delays
sys.mode_no = 3;            % number of modes
sys.event_no = 9;           % number of possible events
sys.tau_no = 1;             % number of time delays


%% Finding a periodic orbit

% Periodic orbit initialisation via simulation
p0 = [0.05, 0.9, 1.0, 1.5, 1.1];   % parameter vector [zeta, om, tau, f_0, eta]
y0 = @(t) [0 5.0 p0(2)*t];  % initial history function
o_init.N = 3;   % number of events
o_init.M = 20;  % Chebyshev mesh dimension
sim_opts.h_act = [1 1 0 0 0 0 1 1 0];   % active events during simulation (no sliding motion)
sim_opts.m0 = 1;    % initial vector field mode
sim_opts.t_end = 50*2*pi/p0(2); % simulation time
[orb0,res]= pwsdde_sim(y0,p0,sys,o_init,sim_opts);
% figure(); plot_sim_res(res,2); title('Transient simulation result');
% figure(); plot_orb(orb_data,'fricosc_p'); title('Simulated guess');


%% Follow the periodic orbit with COCO

% define a COCO compatible problem
p_names = {'$\zeta$', '$\omega$', '$\tau$', '$f_0$', '$\eta$'};
[prob1, data1] = pwsdde_coco_prob(sys,orb0,p_names,5);

% set continuation options
prob1 = coco_set(prob1,'cont', 'NAdapt', 0);
prob1 = coco_set(prob1,'cont','h_min',1e-2,'h_max',2);
prob1 = coco_set(prob1,'cont','ItMX', [100 100]);
% prob1 = coco_set(prob1,'cont','LogLevel',[3,1]);

% add sliding event condition
prob1 = event_slide(prob1,data1);

% create a branch of solutions with COCO
bd1 = coco(prob1,'fro_r1',[],1,{'$\omega$','$|A_1|$','$|\mu_c|$'},...
    [0.0, 2.0]);


%% Follow sliding orbit in 2 parameters

% convert an SL endpoint to a sliding orbit
lab = coco_bd_labs('fro_r1','SL');
orb_sl = orb_get_data('fro_r1',lab(1));
% figure(); plot_orb(orb_sl,'fricosc_p')

% define a new COCO compatible problem to include a sliding condition
[prob_sl, data_sl] = pwsdde_coco_prob(sys,orb_sl,p_names);
prob_sl = cond_slide(prob_sl,data_sl,2);
prob_sl = coco_set(prob_sl,'cont','h_min',1e-2,'h_max',2);

% create a branch of solutions with COCO
bd_sl = coco(prob_sl,'fro_sl',[],1,{'$\omega$','$f_0$',...
    '$|A_1|$','$|\mu_c|$'},[0.5, 1.5]);


%% Plot continuation results with COCO

% read continuation data
fro_r1 = coco_bd_read('fro_r1');
fro_sl = coco_bd_read('fro_sl');

% bifurcation diagram
thm = struct();
thm.special = {'EP','FP','VA','SL'};
thm.VA = {'dk','MarkerSize',7,'MarkerFaceColor','y'};
thm.SL = {'ok','MarkerSize',7,'MarkerFaceColor','g'};
figure(); box on; hold on;
coco_plot_bd(thm,'fro_r1','$\omega$','$f_0$','$|A_1|$')
coco_plot_bd(thm,'fro_sl','$\omega$','$f_0$','$|A_1|$')
hold off

% display an orbit

% lab = coco_bd_labs('fro_r1','SL');
% orb1 = orb_get_data('fro_r1',lab(1));
% figure(); plot_orb(orb1,'fricosc_p')

lab = coco_bd_labs('fro_sl','EP');
orb_sl = orb_get_data('fro_sl',lab(1));
figure(); plot_orb(orb_sl,'fricosc_p')

% animate orbit progression
% anim_bd_orb('fro_r1','fricosc_p')