function [strain_model] = forward_burger(theta,xdata)

E1    = theta(1);
eta1  = theta(2);
E2    = theta(3);
eta2  = theta(4);

% creep part for Kelvin, Maxwell, and Burger
for i= 1:60
compliance_c_k(i) = (1/E1)*(1-exp(-xdata(i)*(E1/eta1)));
compliance_c_m(i) = (1/E2) + (1/eta2)*(xdata(i));
compliance_c_b(i) = compliance_c_m(i) + compliance_c_k(i);
end

% recovery part for Kelvin, Maxwell, and Burger
for i= 61:180
compliance_r_k(i) = (1/E1)*(exp(((E1/eta1)*(60)))-1)*(exp(-xdata(i)*(E1/eta1)));
compliance_r_m(i) = (1/eta2)*(60);
compliance_r_b(i) = compliance_r_m(i) + compliance_r_k(i);
end

for i= 1:60
strain_c(i) = 0.5*compliance_c_b(i);
end
for i= 1:120
strain_r(i) = 0.5*compliance_r_b(i+60);
end

% Sum of creep and recovery
strain_model    = [strain_c strain_r]';
end




