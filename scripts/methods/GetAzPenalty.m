% DESCRIPTION: Computes the azimuthal energy contribution:
function EPair = GetAzPenalty
%Set the parameters for the angular distribution: 
t_w = pi/9; %Width of the strip 
t_d = (2*pi)/7;%(5*pi)/18; %Separation between strip centerlines 
t_flat = 0.06; %0.17; %0.06; %'Flat' width of the strip around the peak 
q = 100; % q^(-1) is the value of EAz at the strip edge
p = 1.0004; %p^(-1) is the value of EAz at the edges of the flat strip part
t_p = 0.5; %Constant shift factor for EAz. !IF YOU CHANGE THIS MAKE ADJUSTMENTS TO ClusterEnergy.m!

%Define the angular component of the energy: 
f = @(t) sin((pi*t)/(2*t_d));
b = (log(q-1)*f(t_flat) - log(p-1)*f(t_w))/(f(t_w)-f(t_flat));
a = (log(q-1)+b)/f(t_w);
EAz = @(az) 1./ (1 + exp(a.*abs( sin(pi.*(az-t_p)./t_d) ) - b)); 

%Define the pairwise interaction energy: 
EPair = @(t_i,t_j) EAz(t_i).*EAz(t_j);
end
