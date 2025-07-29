%DESCRIPTION: receives vectorised inputs and computes the electrostatic interaction
%energy curve. Adds up all interactions calculated at the same stagger
%values. 

function [E,SVec] =  GetAllDH_E(dist_vec,z_charge,pKa_vec,RowStatic,RowDyn,V_lc,pH,I,a0,resType)
%Compute the Debye-Huckel energies: 
E = GetDebyeE_vec(dist_vec,I,z_charge{1}(RowStatic),z_charge{2}(RowDyn),...
    pH,pKa_vec{1}(RowStatic),pKa_vec{2}(RowDyn),a0,resType{1}(RowStatic),resType{2}(RowDyn));
%Add up the energies at the same lattice points: 
[SVec,~,ix] = unique(reshape(V_lc.',[],1));
E = accumarray(ix,E);

%Transpose the outputs, so that they are row vectors: 
SVec = SVec'; E = E';
end
