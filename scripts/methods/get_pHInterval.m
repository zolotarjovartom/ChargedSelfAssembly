% DESCRIPTIONL Compute the ionisation interval on the pH scale for a set of
% residues 

% INPUTS: 

% alphaRange - 2 element vector denoting the starting and ending ionisation
% fractions. 

% ResType - 'A' or 'B' denoting whether the residue is an acid or a base

% NOTE: 

% All the inputs need to be row vectors! 

function [pH_window,f_act_out] = get_pHInterval(alphaRange,ResType,z,eps_r,T,I,pKa)
%Set constant(s):
A = (eps_r*T)^(-1.5)*1.82e6;

%Compute the activity coefficient:
I = I*1e-3; %Convert the units of I into M from mM 
f_act = -A.*z.^(2).*(sqrt(I)/(1+sqrt(I)) - 0.2*I);

%Identify the acidic residues 
TF_acid = matches(ResType, 'A');

%Define the output: 
pH_window = zeros(numel(ResType),2);

%Compute the interval half-width along the pH scale: 
delta_pH = 0.5*log10(alphaRange(:,2).*alphaRange(:,1).^(-1).*(1-alphaRange(:,1))...
    .*(1-alphaRange(:,2)).^(-1));

%Compute the pH window: 
pH_window(TF_acid,:) = [pKa(TF_acid) + f_act(TF_acid) - delta_pH(TF_acid),...
    pKa(TF_acid) + f_act(TF_acid) + delta_pH(TF_acid)]; 
pH_window(~TF_acid,:) = [pKa(~TF_acid) - f_act(~TF_acid) - delta_pH(~TF_acid),...
    pKa(~TF_acid) - f_act(~TF_acid) + delta_pH(~TF_acid)]; 

%Output the activity coefficient: 
f_act_out = unique(f_act);
if numel(f_act_out) > 1 
    error('Activity coefficients differ between residues!')
end
end