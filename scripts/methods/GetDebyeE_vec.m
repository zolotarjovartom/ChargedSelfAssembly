%DESCRIPTION: At a given pH and I value, computes the Debye-Huckel
%interaction energy, accepting pKa's, charges and distances as vector
%inputs

%NOTE: a0 MUST BE IN METERS 

function E = GetDebyeE_vec(r,I,zi,zj,pH,pKa_i,pKa_j,a0,resType_i,resType_j)
%Define/compute the relevant constants: 
e_charge = 1.602177e-19; %elementary charge
eps_0 = 8.854188e-12; %vaccuum permittivity
eps_r = 76.5; %relative permittivity, take value for water at 303.15K
k_b = 1.380649e-23; %Boltzmann constant
T = 303.15; %Temperature in K 
NA = 6.022141e23; %Avogadro's number
l_bjer = e_charge^(2)/(4*pi*eps_0*eps_r*k_b*T); %Bjerrum length in meters
k = sqrt(8*pi*l_bjer*NA*I); %Inverse Debye Length in meters^(-1)
fac_deb = exp(k*a0)/(1+k*a0); %Debye Factor 

%Compute the effective charges: 
zi = GetEffCharge(pH,pKa_i,I,zi,eps_r,T,resType_i);
zj = GetEffCharge(pH,pKa_j,I,zj,eps_r,T,resType_j);

%Compute the Debye-Huckel interaction energy (in units of kbT):
E = zeros(size(r));

%Identify distances greater or equal than a0:
TF_DH = r >= a0*10^(10); %r is an ANGSTROMS 

%Compute the DH interaction: 
E(TF_DH) = (zi(TF_DH).*zj(TF_DH).*l_bjer.*fac_deb.*exp(-k.*r(TF_DH).*10^(-10)))./(r(TF_DH).*10^(-10));

%For the remaining distances closer than a0, assign them the same value as at r = a0: 
E(~TF_DH) = (zi(~TF_DH).*zj(~TF_DH).*l_bjer.*fac_deb.*exp(-k.*a0))./a0;
end


%Compute the effective charge on the residues: 

%INPUTS: 

%SCALARS: 

%pH, I (ionic strength) 

%eps_r, T - relative permittivity and temperature 

%VECTORS: 

%pKa, z - pKa value and charge of the residue 
function z_eff = GetEffCharge(pH,pKa,I,z,eps_r,T,resType) 
%Set constant(s):
A = (eps_r*T)^(-1.5)*1.82e6;

%Compute the activity coefficient:
I = I*1e-3; %Convert the units of I into M from mM 
f_act = -A.*z.^(2).*(sqrt(I)/(1+sqrt(I)) - 0.2*I);

%Allocate a container for fraction of charged residues: 
a_d = zeros(size(f_act));

%Identify acidic residues: 
%TF_acid = pKa < 7.0; 
TF_acid = matches(resType, 'A');

%Compute the fraction of charged ACIDIC residues: 
a_d(TF_acid) = 1./(1 + 10.^(pKa(TF_acid) - pH + f_act(TF_acid)));

%Compute the fraction of charged BASIC residues: 
a_d(~TF_acid) = 1 - 1./(1 + 10.^(pKa(~TF_acid) - pH - f_act(~TF_acid)));

%Compute the effective charge: 
z_eff = a_d.*z;
end 