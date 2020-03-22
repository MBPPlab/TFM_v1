% function output = deflat_part_est(input, liste_est, wn)
%
% EN/ function deflation of particles 
% already estimated & validated (result_ok)
%
%
% FR/ fonction deflation des particules 
% deja estimees & validees (result_ok)


function output = deflat_part_est(input, liste_est, wn)

%%
if (nargin < 3)
     wn = 9 ;
else
     wn = wn - 2 ; %% si probleme sur les bords
end%if


% [idim, jdim] = size(input) ;
nb_part = size(liste_est, 1) ;

output = input ;

%% parametre dans liste_est :
%% liste_param = [num, i, j, alpha, sig^2, rayon, ok]

for part=1:nb_part
if (liste_est(part, 7) == 1)
     i0 = liste_est(part, 2) ;
     j0 = liste_est(part, 3) ;
     alpha = liste_est(part, 4) ;
     r0 = liste_est(part, 6) ;

     pos_i = round(i0) ;
     dep_i = i0 - pos_i ;
     pos_j = round(j0) ;
     dep_j = j0 - pos_j ;

     alpha_g = alpha * gausswin2(r0, wn, wn, dep_i, dep_j) ; 
 
     dd = (1:wn) - floor(wn/2) ;
     di = dd + pos_i ;
     dj = dd + pos_j ;
     output(di, dj) = output(di, dj) - alpha_g ;
end%if
end%for


end%function

