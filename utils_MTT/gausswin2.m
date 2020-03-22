% function g = gausswin2(sig, wn_i, wn_j, offset_i, offset_j)
%
% EN/ generation of a 2D Gaussian
% of width sig, of power 1
% in a squared window wn
% ex: g = gausswin2(1.3, 8)
%
% if two supplementary arguments 
% then offset_x & offset_y (in pixel)
% ex: g = gausswin2(1.3, 8, 0.3, -0.8)
%
% if an odd number of arguments, different size 
% of the window in i & j:
% ex: g = gausswin2(1.3, 8, 9, 0.3, -0.8)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FR/ generation d_une gaussienne 2D
% de largeur sig, de puissance 1
% sur une fenetre carree wn
% ex : g = gausswin2(1.3, 8)
%
% si deux arguments supplementaires
% alors offset_x & offset_y (en pixel)
% ex : g = gausswin2(1.3, 8, 0.3, -0.8)
%
% si un nombre d_arguments impaire la taille 
% de la fenetre differente en i & j :
% ex : g = gausswin2(1.3, 8, 9, 0.3, -0.8)


function g = gausswin2(sig, wn_i, wn_j, offset_i, offset_j)

if (nargin < 5)
  offset_i = 0.0 ;
  offset_j = 0.0 ;
end%if

if (nargin < 3)
     wn_j = wn_i ;
end%if

refi = 0.5 + (0:(wn_i-1)) - wn_i/2 ;
i = refi - offset_i ;
refj = 0.5 + (0:(wn_j-1)) - wn_j/2 ;
j = refj - offset_j ;
ii = i' * ones(1,wn_j) ; %'
jj = ones(wn_i,1) * j ;   

%%% puissance unitaire
g = (1/(sqrt(pi)*sig)) * exp(-(1/(2*sig^2))*(ii.*ii + jj.*jj)) ;

end %function

