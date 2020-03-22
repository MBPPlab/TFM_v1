%  function [carte_MV, liste_detect, detect_pfa] 
%     = carte_H0H1_1vue(im, rayon, wn_x, wn_y, s_pfa)
%
% EN/ hypothesis map: no particle /
% presence of a particle in the centre of the
% research window, under the hypotheses of
% Gaussian iid background noise with
% the signal of the particle being a Gaussian.
% The amplitude of the signal (particle)
% is unknown, as well as the power of
% the noise. By hypothesis, the background is 
% supposed constant, it is thus preferable
% to limit the size of the window (<12)
%
% im: input image
% rayon = width of the Gaussian
% wn_x: width of the detection window
% wn_y: width in y, squared if absent
% s_pfa: threshold for the detection pfa
%
% exemple: cmv = carte_H0H1_1vue(im, 1.3, 8) ;
% exemple: cmv = carte_H0H1_1vue(im, 1.3, 8, 7) ;
% exemple: [cmv, liste_detect] = carte_H0H1_1vue(im, 1.3, 8, 8, 33) ;
% by default (absence of parameter) pfa = 1E-7
%
% <liste_detect> returns the matrix of coordinates i,j of
% the detected particles as well as the amplitude of the peaks and
% the variance of the noise (local)
% liste_detect = [num, i, j, ampli, var_bruit]
%
% remark: better choose the size of the window even
% if the size of the image is odd in the considered direction
% and reciprocally if even (avoid shifts...)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FR/ Carte d_hypothese : pas de particule /
%  presence de particule au centre de la 
%  fenetre de balayage, sous l_hypothese que
%  de bruit de fond gaussien iid et que
%  le signal de la particule est une gaussienne.
%  L_amplitude du signal (particule)
%  est inconnue, ainsi que la puissance
%  du bruit. Par hypothese, le fond est 
%  suppose constant, il est donc preferable
%  de limiter la taille de la fenetre (<12)
%
%  im : image d_entree
%  rayon = largeur de la gaussienne
%  wn_x : taille de la fenetre de detection
%  wn_y : taille en y, carre si absent
%  s_pfa : seuil pour la pfa de detection
%
% exemple : cmv = carte_H0H1_1vue(im, 1.3, 8) ;
% exemple : cmv = carte_H0H1_1vue(im, 1.3, 8, 7) ;
% exemple : [cmv, liste_detect] = carte_H0H1_1vue(im, 1.3, 8, 8, 33) ;
% par defaut (absence de parametre) la pfa = 1E-7
%
% <liste_detect> renvoie la matrice des coordonnées i,j des
% particules detectees ainsi que l_amplitude des pics et
% la variance du bruit (local)
% liste_detect = [num, i, j, ampli, var_bruit]
%
% remarque : prendre plutot une taille de fenetre paire
% si la taille de l_image est paire dans la direction consideree
% et reciproquement si impaire (evite des decalages...)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EN/ detection at given false alarm rate
% chi2 law with degree of liberty = 1 (3-2)
% 
%
% FR/ detection a taux de fausse alarme donnee
% loi du chi2 de degree de liberte = 1 (3-2)
%
%     s_pfa      1-pfa
%    3.77060   0.94784
%    6.63106   0.98998
%   10.79172   0.99898
%   15.00000   0.999892488823270
%   20.00000   0.999992255783569
%   24.00000   0.999999036642991  (1E-6)
%   25.00000   0.999999426696856
%   28.00000   0.999999878684549  (1E-7)
%   30.00000   0.999999956795369
%   33.00000   0.999999990784113  (1E-8)
%   37.50000   0.999999999085870  (1E-9)
%   40.00000   0.999999999746037
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [carte_MV, liste_detect, detect_pfa] = carte_H0H1_1vue(im, rayon, wn_x, wn_y, s_pfa)

if (nargin < 5) 
  s_pfa = 28.0 ;
end%if ;
if (nargin < 4) 
  wn_y = wn_x ;
end%if ;



[N,M] = size(im) ; 
T = wn_x*wn_y ; % nombre de pixel dans la fenetre 


%% Hypothese H0
%% pas de particule dans la fenetre
m = ones(wn_x,wn_y) ;
hm = expand_w(m, N, M) ;
tfhm = fft2(hm) ;
tfim = fft2(im) ;
m0 = real(fftshift(ifft2(tfhm .* tfim))) /T ;

im2 = im .* im ;
tfim2 = fft2(im2) ;
Sim2 = real(fftshift(ifft2(tfhm .* tfim2)));

%% H0 = T/2*log(2*pi*sig0^2)-T/2 ;
T_sig0_2 = Sim2 - T*m0.^2 ;

%% Hypothèse H1
%% une particule est au centre de la fenetre
%% amplitude inconnue, rayon fixe

%% generation masque gaussien de largeur (sigma)
%% egal a rayon

%%g = gausswin2(rayon, wn_x, wn_y) ;
g = gausswin2(rayon, wn_x, wn_y, 0, 0) ;
gc = g - sum(g(:))/T ;
Sgc2 = sum(gc(:).^2) ;
hgc = expand_w(gc, N, M) ;
tfhgc = fft2(hgc) ;

alpha = real(fftshift(ifft2(tfhgc .* tfim))) / Sgc2 ;

%% H1 = T/2*log(2*pi*sig1^2)-T/2 ;
%%sig1_2 = sig0_2 - alpha.^2 * Sgc2 / T ;

%% pour test
%sig1_2 = T_sig0_2/T - alpha.^2 * Sgc2 / T ;
%%imagesc(T_sig0_2/T);
%imagesc(sig1_2);
%%imagesc(sig1_2 ./ (T_sig0_2/T));

%%  carte_MV = -0.5*(H0 - H1) ;
% carte_MV = - T * log(1 - (Sgc2 * alpha.^2) ./ T_sig0_2) ; 
test = 1 - (Sgc2 * alpha.^2) ./ T_sig0_2 ;
test = (test > 0) .* test + (test <= 0) ;
carte_MV = - T * log(test) ; 

%% detection et recherche des maximas
%% s_pfa = 28 ;
detect_masque = carte_MV > s_pfa ;
if (sum(detect_masque(:))==0)
    warning('No target detected !') ; %#ok
    liste_detect = zeros(1,6) ;
    detect_pfa = zeros(size(detect_masque)) ; % ajout AS 4/12/7
else
    detect_pfa = all_max_2d(carte_MV) & detect_masque ;

    [di, dj] = find(detect_pfa) ;
    n_detect = size(di, 1) ;
    vind = N*(dj-1)+di ;
    valpha = alpha(:) ;
    alpha_detect = valpha(vind) ;

    sig1_2 = ( T_sig0_2 - alpha.^2 * Sgc2 ) / T ;
    vsig1_2 = sig1_2(:) ;
    sig2_detect = vsig1_2(vind) ;

    %% g de puissance unitaire
    %%RSBdB_detect = 10*log10(alpha_detect.^2  ./ sig2_detect) ;

    %%liste_detect = [(1:n_detect)', di, dj, alpha_detect, sqrt(sig2_detect), RSBdB_detect] ;
    liste_detect = [(1:n_detect)', di, dj, alpha_detect, sig2_detect, rayon*ones(n_detect,1),ones(n_detect,1)] ;

end%if
end %function




