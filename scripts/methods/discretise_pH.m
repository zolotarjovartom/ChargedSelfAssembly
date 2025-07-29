% DESCRIPTION: Discretises the pH input, according to the desired step size
% in the ionisation fraction

% INPUTS: 

% pH_int - 2 element vector describing the pH interval that needs to be
% discretised 

function disc_pH = discretise_pH(pH_int,z,I,eps_r,T,pKa,ResType,...
    delta_alpha,alphaBounds)
%Initialise Output vector: 
disc_pH = zeros(2*round(diff(alphaBounds)/delta_alpha,0),1);
disc_pH(1) = pH_int(1);

%Set the constant(s):
A = (eps_r*T)^(-1.5)*1.82e6;

%Compute the activity coefficient:
I = I*1e-3; %Convert the units of I into M from mM 
f_act = -A.*z.^(2).*(sqrt(I)/(1+sqrt(I)) - 0.2*I);

%Identify the acidic residues 
TF_acid = matches(ResType, 'A');

%Define some useful quantities: 
delta = pKa' + 1*f_act.*TF_acid - 1*f_act.*~TF_acid;
delta_alpha_vec = TF_acid.*delta_alpha + ~TF_acid.*(-delta_alpha);

%Obtain the pH discretisation recursively: 
alpha = zeros(numel(pKa),1);

%Initialise the loop variable: 
idx = 1; 

%Check if we have reached the upper end of the interval
while disc_pH(idx) < pH_int(2)
    num = 10.^(delta - disc_pH(idx)) + 1;
    denom = TF_acid.*delta_alpha_vec.*(10.^(delta-disc_pH(idx)) + 1) ...
        - ~TF_acid.*delta_alpha_vec.*(10.^(delta-disc_pH(idx)) + 1) +1;
    
    %Compute the acid dissociation at the current pH step: 
    alpha(TF_acid) = (1 + 10.^(pKa(TF_acid)'-disc_pH(idx)+1*f_act(TF_acid))).^(-1);  
    alpha(~TF_acid) = 1 - (1 + 10.^(pKa(~TF_acid)'-disc_pH(idx)-1*f_act(~TF_acid))).^(-1);
    inc_alpha = alpha + delta_alpha_vec;
    TF = or(inc_alpha > 1,inc_alpha < 0);
    
    %Identify curves for which we automatically satisfy the step-size
    %condition: 
    pH_candidates = delta - log10(num./denom - 1);

    %Compute the next discrete point: 
    disc_pH(idx + 1) = min(pH_candidates(~TF));

    %Increment the loop variable: 
    idx = idx + 1; 
end

%Remove unnecessary zeros: 
disc_pH(disc_pH == 0) = [];
end