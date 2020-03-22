% function [lestime, ldetect] 
%    = detect_et_estime_part_1vue(input, wn, r0, pfa, pas_ijr)
%
% EN/ INPUTS
%
% wn size of thewindow for the detection and the estimation
% by default: wn = 13
%
% r0 radius for the detection
% by default: r0 = 1.0
%
% pfa, by default 28 (cf. carte_H0H1_1vue.m)
%
% pas_ijr vector 1x3 setting the research steps
% by default: pas_ijr = [0.125 0.125 0.08]
%
% 
% OUTPUTS
%
% ldetect = [num, i, j, alpha, sig^2]
% lestime = [num, i, j, alpha, sig^2, rayon, ok]
% 
% radius at half max: rayon_mih = sqrt(2*ln(2))*rayon
%
% alpha: amplitude of the Gaussian of power 1
% power of the signal P = alpha^2
%
% power of the noise: sig^2
%
% amplitude max of the signal of the Gaussian
% alpha_max = alpha / (sqrt(pi)*rayon))
%
% ok = 0 if during the estimation one gets a solution
% on one side of the research area (i,j,r)
%
% Nestime nb of detected particles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/ INPUTS
%
% wn taille de la fenetre pour la detection et l_estimation
% par defaut : wn = 13
%
% r0 rayon pour la detection
% par defaut : r0 = 1.0
%
% pfa, par defaut 28 (cf. carte_H0H1_1vue.m)
%
% pas_ijr vecteur 1x3 fixant les pas de recherche
% par defaut : pas_ijr = [0.125 0.125 0.08]
%
% 
% OUTPUTS
%
% ldetect = [num, i, j, alpha, sig^2]
% lestime = [num, i, j, alpha, sig^2, rayon, ok]
% 
% rayon a mi hauteur : rayon_mih = sqrt(2*ln(2))*rayon
%
% alpha : amplitude de la gaussienne de puissance 1
% puissance du signal P = alpha^2
%
% puissance du bruit : sig^2
%
% amplitude max du signal de gaussienne
% alpha_max = alpha / (sqrt(pi)*rayon))
%
% ok = 0 si lors de l_estimation on a une solution
% sur l_un des bords des domaines de recherche (i,j,r)
%
% Nestime nb particules detectees


function [lestime, ldetect, d, Nestime] = detect_et_estime_part_1vue(input, wn, r0, pfa)%%, pas_ijr)

if (nargin < 2)
     wn = 9 ;
end%if
if (nargin < 3 )
     r0 = 1.3 ;
end%if
if (nargin < 4)
     pfa = 28 ;
end%if
% if (nargin < 5)
%      pas_ijr = [0.125 0.125 0.08] ;
% end%if


[Ni, Nj] = size(input) ;

%% positions des parametres
Nparam = 7 ;
detect_i = 2 ;
detect_j = 3 ;
alpha = 4 ;
% sig2 = 5 ;

%% Pour Matlab
stderr = 1 ;

%% detection pour un rayon moyen r0
[c,ldetect,d] = carte_H0H1_1vue(input, r0, wn, wn, pfa);

%% estimation pour des pas de recherche
%% interval [-1,1] pour ij
%% interval [0.6,1.4] pour le rayon

Ndetect = size(ldetect, 1) ;
if (Ndetect==0)
     lestime = zeros(1,Nparam) ;
     warning('no particle detected => exit') ; % #ok
     Nestime = 0 ;
     return ;
end%if

Nestime = 0 ;
bord = ceil(wn/2) ;

fprintf(stderr,'nb T+cted : %d\n', Ndetect);
for n=1:Ndetect
     test_bord = (ldetect(n,detect_i) < bord) || (ldetect(n,detect_i) > Ni-bord) || ...
                 (ldetect(n,detect_j) < bord) || (ldetect(n,detect_j) > Nj-bord) ;
     if ((ldetect(n,alpha) > 0.0) && (~test_bord) )
       %%fprintf(stderr, "%d/%d\n", Nestime, n) ;
        Nestime = Nestime + 1 ;
        lestime(Nestime, :) = estim_param_part_GN(input, wn, ldetect(n,:),r0) ;
     end%if
end%for

% a la bonne taille
if (Nestime==0)
  lestime = zeros(1,Nparam) ;
else
  lestime = lestime(1:Nestime,:) ;
end%if
fprintf(stderr, 'nb part detected : %d, estimated : %d, correct : %d\n', Ndetect, Nestime, sum(lestime(:,7))) ;

end%function

