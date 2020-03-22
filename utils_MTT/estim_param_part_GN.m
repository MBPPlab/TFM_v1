%function liste_param
%   = estim_param_part_GN(im, wn, liste_info_param, r0, bornes_ijr)
%
% EN/ sub-pixel estimate of the peak position by Gauss Newton regression
% liste_info_param is a line of the matrix liste_detect
% wn odd
%
%
% FR/ estimation sub-pixel de la position du pic par Gauss Newton
% liste_info_param est une ligne de la matrice liste_detect
% wn impaire


%%% version 11 07 05

function liste_param = estim_param_part_GN(im, wn, liste_info_param, r0, bornes_ijr)


if (nargin < 5)
     bornes_ijr(1) = -1.5 ;
     bornes_ijr(2) = 1.5 ;
     bornes_ijr(3) = -1.5 ;
     bornes_ijr(4) = 1.5 ;
     bornes_ijr(5) = 0.3 ;
     bornes_ijr(6) = 3.0 ;
end%if

Pi = liste_info_param(2) ;
Pj = liste_info_param(3) ;
di = (1:wn)+Pi-floor(wn/2) ;
dj = (1:wn)+Pj-floor(wn/2) ;
im_part = im(di, dj) ;


r = r0 ;
i = 0.0 ;
j = 0.0 ; 
% dr = 1 ;
di = 1 ;
dj = 1 ;
fin = 0.01 ;
sig2 = inf ;
cpt = 0 ;
test = 1 ;
ITER_MAX = 50 ;
while (test)
   %%[r, i, j, dr, di, dj, alpha, sig2] = deplt_GN_estimation (r, i, j, im_part) ;
   [r, i, j, dr, di, dj, alpha, sig2] = deplt_GN_estimation (r, i, j, im_part, sig2, di, dj) ;
   cpt = cpt + 1 ;
   test = max([abs(di), abs(dj), abs(dr)]) > fin ;
   if (cpt > ITER_MAX) 
     test = 0 ;
   end%if

   %% on stop si l_on sort des bornes  
   result_ok = ~((i < bornes_ijr(1)) || (i > bornes_ijr(2)) || ...
		 (j < bornes_ijr(3)) || (j > bornes_ijr(4)) || ...
		 (r < bornes_ijr(5)) || (r > bornes_ijr(6)) ) ;
   test = test & result_ok ;

end%while


% liste_info_param = [num, i, j, alpha, sig^2]
% liste_param = [num, i, j, alpha, sig^2, rayon, ok]

liste_param = [liste_info_param(1), ...
               Pi+i , ...
               Pj+j , ...
               alpha , ...
               sig2 , ...
	           r , ...		      
               result_ok ];

end%fonction

