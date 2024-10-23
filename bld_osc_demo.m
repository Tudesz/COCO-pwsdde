clear; clc; close all;
warning off backtrace
%% COCO compatible continuation problem of 1DoF bilinear delayed oscillator
% considering harmonic excitation and delayed velocity feedback control
% using nondimensionalized form of governing equation of motion
% further details on the model at: https://doi.org/10.1016/j.cnsns.2019.105095

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
addpath('_system def/bilinear oscillator');
sys.f = 'bldosc_f';     % vector field modes
sys.e = 'bldosc_e';     % event functions
sys.tau = 'bldosc_tau'; % time delays
sys.tau_no = 1;         % number of time delays
sys.mode_no = 2;        % number of modes
sys.event_no = 3;       % number of possible events


%% Finding a periodic orbit

% Periodic orbit initialisation via simulation
p0 = [3.8 0.12 0.01 29 0.9*0.8^2/1.26 0.8]; % parameter vector [tau, k, zeta, beta, f, om]
y0 = @(t) [0 0 cos(p0(6)*t) sin(p0(6)*t)];  % initial history function
o_init.N = 2;                   % number of events
o_init.M = 100;                 % initial guess mesh dimension
sim_opts.h_act = [1 1 0];       % active events during simulation
sim_opts.t_end = 20*2*pi/p0(6); % simulation time
[orb0,res]= pwsdde_sim(y0,p0,sys,o_init,sim_opts);
% figure(); plot_sim_res(res,1); title('Transient simulation result');
% figure(); plot_orb(orb_data,'bldosc_p'); title('Simulated guess');


%% Follow the periodic orbit with COCO

% define a COCO compatible problem
p_names = {'$\tau$', '$k$', '$\zeta$', '$\beta$', '$f$', '$\omega$'};
[prob1, data1] = pwsdde_coco_prob(sys,orb0,p_names);

% set continuation options
prob1 = coco_set(prob1,'cont', 'NAdapt', 0);
prob1 = coco_set(prob1,'cont','h_min',1e-2,'h_max',2);
prob1 = coco_set(prob1,'cont','ItMX', [0 100]);
% prob1 = coco_set(prob1,'cont','LogLevel',[3,1]);

% add vanishing event condition
prob1 = event_vanish(prob1,data1,1e-2);

% create a branch of solutions with COCO
bd1 = coco(prob1,'bld_r1',[],1,{'$\tau$','T_min','$|\mu_c|$'},[3.5, 5]);

% create a second branch of solutions with COCO
bd2 = coco(prob1,'bld_r2',[],1,{'$k$','T_min','$|\mu_c|$'},[0.1, 0.20]);


%% Follow grazing orbit in 2 parameters

% convert a VA endpoint to a grazing orbit
lab = coco_bd_labs('bld_r1','VA');
orb_va = orb_get_data('bld_r1',lab(1));
orb_gr = orb_convert(orb_va,data1,0);
% figure(); plot_orb(orb_gr,'bldosc_p')

% define a COCO compatible problem including a grazing condition
[prob_gr, data_gr] = pwsdde_coco_prob(sys,orb_gr,p_names);
prob_gr = cond_graze(prob_gr,data_gr,1);
prob_gr = coco_set(prob_gr,'cont','h_min',1e-2,'h_max',2);

% create a branch of solutions with COCO
bd_gr = coco(prob_gr,'bld_gr',[],1,{'$\tau$', '$k$', 'T_min',...
    '$|\mu_c|$'},[1.5, 5]);


%% Plot continuation results with COCO

% read continuation data
bld_r1 = coco_bd_read('bld_r1');
bld_r2 = coco_bd_read('bld_r2');
bld_gr = coco_bd_read('bld_gr');

% bifurcation diagram
thm = struct();
thm.special = {'FP','VA'};
thm.VA = {'dk','MarkerSize',7,'MarkerFaceColor','y'};
figure(); box on; hold on
coco_plot_bd(thm,'bld_r1','$\tau$','$k$')
coco_plot_bd(thm,'bld_r2','$\tau$','$k$')
coco_plot_bd(thm,'bld_gr','$\tau$','$k$')
hold off

% display an orbit
lab = coco_bd_labs('bld_r1','VA');
orbp = orb_get_data('bld_r1',lab(1));
% lab = coco_bd_labs('bld_gr','EP');
% orbp = get_orb_data('bld_gr',lab(1));
figure(); subplot(1,2,1); plot_orb(orbp,'bldosc_p')
subplot(1,2,2); plot_spectrum(orbp.mu);

% animate orbit progression
% anim_bd_orb('bld_r1','bldosc_p')
% anim_bd_spectrum('bld_r1')
