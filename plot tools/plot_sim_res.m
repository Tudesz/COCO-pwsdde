function plot_sim_res(res,ind,p)
%PLOT_SIM_RES Plot the transient simulation results obtained via
% sim_pwsddde.m or sim_pws_sd_dde.m
% Input:
%   res: DDE solver output structure for troubleshooting
%   ind: index of state dimension to plot (default all), if a name of a
%       funciton is given plot according to that
%   p: system parameter vector neccessary for the plotting function

colors = [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330;
    0.6350 0.0780 0.1840];
 
for i=1:length(res)
    color = colors(mod(res(i).mode,size(colors,1)),:);
    if nargin<2
        plot(res(i).x,res(i).y(:,:),'Color',color);
    elseif (ischar(ind) || isa(ind,'function_handle')) && nargin>2
        [px,py] = feval(ind,res(i).x,res(i).y,p);
        plot(px,py,'Color',color);
    else
        plot(res(i).x,res(i).y(ind,:),'Color',color);
    end
    if i==1
        hold on
    end
end
if nargin<2 || ~ischar(ind)
   xlabel('$t$');  ylabel('$x$'); 
else
    xlabel('$f_1(t,x)$'); ylabel('$f_2(t,x)$');
end
box on; hold off

end