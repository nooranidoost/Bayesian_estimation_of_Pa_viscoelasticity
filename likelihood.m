function [ss] = likelihood(theta,data)

strain_obs     = data.ydata; 
xdata          = data.xdata;

[strain_model] = forward_burger(theta,xdata);
ss = sum((strain_model - strain_obs).^2);
