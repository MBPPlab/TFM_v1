% function [lestime, input_deflt, dfin]
%   = detect_et_estime_part_1vue_deflt(input, wn, r0, pfa, n_deflt)
%
% EN/
% ldetect = [num, i, j, alpha, sig^2]
% lestime = [num, i, j, alpha, sig^2, rayon, ok]
% 
% radius at half max: rayon_mih = sqrt(2*ln(2))*rayon
%
% alpha: amplitude of the Gaussian of power 1
% power of the signal P = alpha^2
% power of the noise: sig^2
%
% amplitude max of the signal of the Gaussian
% alpha_max = alpha / (sqrt(pi)*rayon))
%
% ok = 0 if during the estimation one gets a solution
% on one side of the research area (i,j,r)
%
% macro which uses the deflation
%
% input_deflt: result after deflation
% if perfect, there must be only noise remaining!!!
% dfin: binary map of detection
% n_deflt: number of iteration of deflation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/
% ldetect = [num, i, j, alpha, sig^2]
% lestime = [num, i, j, alpha, sig^2, rayon, ok]
% 
% rayon a mi hauteur : rayon_mih = sqrt(2*ln(2))*rayon
%
% alpha : amplitude de la gaussienne de puissance 1
% puissance du signal P = alpha^2
% puissance du bruit : sig^2
%
% amplitude max du signal de gaussienne
% alpha_max = alpha / (sqrt(pi)*rayon))
%
% ok = 0 si lors de l_estimation on a une solution
% sur l_un des bords des domaines de recherche (i,j,r)
%
% macro qui utilise la deflation
%
% input_deflt : resultat apres deflation
% si parfait, il ne doit rester que du bruit !!!
% dfin : carte binaire de detection
% n_deflt : nombre d_iteration de deflation


function [lestime, input_deflt, dfin] = detect_et_estime_part_1vue_deflt(input, wn, r0, pfa, n_deflt)

[lest,ldec,dfin] = detect_et_estime_part_1vue (input, wn, r0, pfa) ;
input_deflt = deflat_part_est(input, lest, wn);
lestime = lest ;

for n=1:n_deflt
   [l,ld,d,N] = detect_et_estime_part_1vue (input_deflt, wn, r0, pfa) ;
   if (N==0), return, end
   lestime = [lestime ; l] ;
   %%dfin += d 
   dfin = dfin | d ;
   input_deflt = deflat_part_est(input_deflt, l, wn);
end%for

end%function

