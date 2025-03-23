clear; clc; close all;
warning off backtrace
%% Example continuation problem of turning with a robotic arm
% governed by a neutral delay differential equation
% more data on the model at https://doi.org/10.1016/j.ijnonlinmec.2022.104239

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
addpath('_system def/nonlinear robot arm');
sys.f = 'nlne_robot_f';     % vector field modes
sys.e = 'nlne_robot_e';     % event functions
sys.tau = 'nlne_robot_tau'; % time delays
sys.tau_no = 1;             % number of time delays
sys.mode_no = 1;            % number of modes
sys.event_no = 1;           % number of possible events


%% Finding a periodic orbit

% Periodic orbit initialisation via simulation
p0 = [0.05 1 1 -1 0.48 0 0.1]; % parameter vector [chi, r, gamma, mu, k, knl, tau]
y0 = @(t) [p0(7); 0; 0; -0.1]; % initial history function
o_init.N = 2;               % number of events
o_init.M = 20;              % Chebyshev mesh dimension
sim_opts.m0 = 1;            % first event
sim_opts.t_end = 100;       % simulation time
[orb0,res]= pwsdde_sim(y0,p0,sys,o_init,sim_opts);
figure(); plot_sim_res(res,2); title('Transient simulation result');
% figure(); plot_orb(orb0); title('Simulation guess');


%% Follow the periodic orbit with COCO

% define a COCO compatible problem
p_names = {'$\xi$', '$r$', '$\gamma$', '$\mu$', '$k$', '$k_{nl}$', '$\tau$'};
[prob1, data1] = pwsdde_coco_prob(sys,orb0,p_names);

% COCO numeric options
prob1 = coco_set(prob1,'cont','h_min',1e-3,'h_max',0.5);
prob1 = coco_set(prob1,'cont','ItMX',[100 100]);

% create a branch of solutions with COCO
bd1 = coco(prob1,'nlro_r1',[],1,{'$k$','$|\mu_c|$'},[0.4, 0.5]);


%% Plot continuation results

% read continuation data
nlro_r1 = coco_bd_read('nlro_r1');

% bifurcation diagram
thm = struct();
thm.special = {'EP','SC'};
thm.SC = {'ob','MarkerSize',7};
thm.ustab = '$|\mu_c|$';
thm.ustabfun =  @(x) (abs(x)<1) + 1;
thm.lspec = {{'b--', 'LineWidth', 1}, {'b-', 'LineWidth', 1}};
figure(); box on; hold on
coco_plot_bd(thm,'nlro_r1','$k$','$|A_1|$')
hold off

% display an orbit
lab = coco_bd_labs('nlro_r1','SC');
orbp = orb_get_data('nlro_r1',lab(1));
figure();
plot_orb(orbp,'nlne_robot_p');

% animate orbit progression
anim_bd_orb('nlro_r1','nlne_robot_p')