function plot_sim_tau(res,tau_ind)
%PLOT_SIM_TAU Plot the state dependent delays obtained via pwsdde_sim.m
% Input:
%   res: DDE solver output structure for troubleshooting
%   tau_ind: index of delayed term to plot (default all)

colors = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840];

figure(); hold on;
for i=1:length(res)
    color = colors(mod(res(i).mode,size(colors,1)),:);
    if nargin<2
        plot(res(i).x,res(i).tau.','Color',color);
    else
        plot(res(i).x,res(i).tau(tau_ind,:).','Color',color);
    end
    if i==1
        hold on
    end
end
xlabel('$t$');  ylabel('$\tau(x(t))$'); 
box on; hold off

end