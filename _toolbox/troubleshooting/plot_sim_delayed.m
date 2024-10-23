function plot_sim_delayed(res,x_ind,tau_ind,tau)
%PLOT_SIM_DELAYED Plot the delayed terms obtained via pwsdde_sim.m for
%troubleshootingS
% Input:
%   res: DDE solver output structure for troubleshooting
%   x_ind: index of state dimension to plot (default all)
%   tau_ind: index of delayed term to plot
%   tau: value of time delay, if it is given shift the results backward in
%       time for validation

if nargin<4
    tau = 0; % no shift by default
end
 
for i=1:length(res)
    if isempty(x_ind)
        plot(res(i).x-tau,squeeze(res(i).Z(:,tau_ind,:)),':k');
    else
        plot(res(i).x-tau,squeeze(res(i).Z(x_ind,tau_ind,:)),':k');
    end
    if i==1
        hold on
    end
end
xlabel('$t$');  ylabel('$x(t-\tau)$'); 
box on; hold off

end