% function [r,v,b] 
%   = affiche_trajectoire(im_t, t, maxi, mini, liste_part, num_t, num_traj)
%
% EN/ displays all particles
% with a color corresponding
% to the particle number
% the diameter of the motif corresponds to the std
% blink is displayed lighter
%
%
% FR/ affiche toute les particules
% avec une couleur correspondant
% au num de la particule
% le diametre du motif correspond a la std
% le blink est affiche en plus clair


function [r,v,b] = affiche_trajectoire(im_t, t, maxi, mini, liste_part, num_t, num_traj)

global tab_param ;
% global tab_var ;
global T_off ;
global Boule_free ;


[M,N] = size(im_t) ; 
dim_max = N*M ;

%% def du cercle
pas_th = pi/18;
theta = pas_th : pas_th : 2*pi ;
cosinus = cos(theta) ;
sinus = sin(theta) ;

%% recalage 0:255
r = 255*(im_t-mini)/(maxi-mini) ;
v = r ;
b = r ;
rr = r ;
vr = v ;
br = b ;

if (nargin<6)
  num_t = t ;
end%if

if (nargin<7)
  num_t = t ;
  num_traj = 1 ; %% affichage numero traj
end%if


%% afffichage du temps
liste_pix = mat2mat3x5(M+1, N, M, num_t) ;
r(liste_pix) = 255 ;
v(liste_pix) = 255 ;
b(liste_pix) = 255 ;

nb_traj_t = size(tab_param, 1) ;

%% particules a afficher (>0)
%% ou a ne pas afficher (<0) (toutes)
if (nargin < 5) 
    liste_part = 1:nb_traj_t ;
end%if
if (liste_part == 0)
  liste_part = 1:nb_traj_t ;
end %if

if (liste_part(1) < 0)
  liste_tmp = 1:nb_traj_t ;
  liste_tmp(-liste_part) = liste_part ;
  liste_part = liste_tmp ;
else

end%if



for traj = liste_part
if ((traj > 0) && (traj <= nb_traj_t))
  if (tab_param(traj, 7*(t-2)+8)>T_off) 
    ii = round(tab_param(traj, 7*(t-1)+3))  ;
    jj = round(tab_param(traj, 7*(t-1)+4))  ;
    stdij = sigij_blink(traj, t-1) ; %% confine
    stdijf = Boule_free * sigij_free_blink(traj, t-1) ; %% libre

    %% couleur du motif
    rand ('seed', traj) ;
    cr = 255*(0.6*rand+0.4) ;
    cv = 255*(0.6*rand+0.4) ;
    cb = 255*(0.6*rand+0.4) ;

    %% intensite blink/non_blink
    if (tab_param(traj, 7*(t-1)+8) > 0)
      coef = 1 ;
    else
      coef = 0.5 ;
    end %if 
    
    %% confine
    pas_theta = 2 ;
    ci = ii + ceil_(stdij*cosinus(1:pas_theta:end)) ;
    sj = jj + ceil_(stdij*sinus(1:pas_theta:end)) ;
    %% limitation zone image
    out_vert = ci < M ;
    ci = ci(out_vert) ;
    sj = sj(out_vert) ;
    mi = ci + sj*M ;
    %% limitation zone image
    mi = 1 + mod(mi, dim_max) ;

    %% on affiche le motif
    r(mi) = coef*cr ;
    v(mi) = coef*cv ;
    b(mi) = coef*cb ;
 
    %% free (domaine de recherche)
    pas_theta = 1 ;
    ci = ii + ceil_(stdijf*cosinus(1:pas_theta:end)) ;
    sj = jj + ceil_(stdijf*sinus(1:pas_theta:end)) ;
    %% limitation zone image
    out_vert = ci < M ;
    ci = ci(out_vert) ;
    sj = sj(out_vert) ;

    mi = ci + sj*M ;
    %% limitation zone image hori
    mi = 1 + mod(mi, dim_max) ;

    %% on affiche le motif
    r(mi) = coef*cr ;
    v(mi) = coef*cv ;
    b(mi) = coef*cb ;
    

    if (num_traj)
    %% affichage numero traj
    liste_pix = mat2mat3x5(max(mi)+1, N, M, traj) ;
    %% limitation zone image hori
    %% liste_pix = 1+mod(liste_pix, dim_max) ;
    liste_pix = 1 + liste_pix(liste_pix<dim_max) ;
    %% on affiche le num
    r(liste_pix) = coef*cr ;
    v(liste_pix) = coef*cv ;
    b(liste_pix) = coef*cb ;
    end %if
  end %if
end %if
end %for

r = [rr, 200*ones(M,1), r] ; %% modif bug 12-11-07 (N->M)
v = [vr, 200*ones(M,1), v] ; %% modif bug 12-11-07
b = [br, 200*ones(M,1), b] ; %% modif bug 12-11-07 

%% on remet aleatoire le generateur aleatoire !!!
rand ('seed', cputime) ;


end %function


function f = ceil_(x)
  sng = sign(x);
  f = sng .* ceil(sng.*x) ;
  f = f .* (abs(x) > 1e-4) ;
end %function 


%% function num2mat3x5
function liste_pix = mat2mat3x5(pos, N, M, num) %#ok

n0 = [...
      1 1 1;...
      1 0 1;...
      1 0 1;...
      1 0 1;...
      1 1 1] ;%#ok
n1 = [...
      0 1 0;...
      0 1 0;...
      0 1 0;...
      0 1 0;...
      0 1 0] ;%#ok
n2= [...
      1 1 1;...
      0 0 1 ;...
      1 1 1;...
      1 0 0;...
      1 1 1] ;%#ok
n3= [...
      1 1 1;...
      0 0 1 ;...
      0 1 1;...
      0 0 1;...
      1 1 1] ;%#ok
n4= [...
      1 0 1;...
      1 0 1 ;...
      1 1 1;...
      0 0 1;...
      0 0 1] ;%#ok
n5 = [...
      1 1 1;...
      1 0 0;...
      1 1 1;...
      0 0 1;...
      1 1 1] ;%#ok
n6= [...
      1 1 1;...
      1 0 0 ;...
      1 1 1;...
      1 0 1;...
      1 1 1] ;%#ok
n7 = [...
      1 1 1;...
      0 0 1;...
      0 1 0;...
      0 1 0;...
      0 1 0] ;%#ok
n8 = [...
      1 1 1;...
      1 0 1;...
      1 1 1;...
      1 0 1;...
      1 1 1] ;%#ok
n9 = [...
      1 1 1;...
      1 0 1;...
      1 1 1;...
      0 0 1;...
      1 1 1] ;%#ok

strnum = num2str(num) ;

cmd = sprintf('[ii, jj] = find(n%s) ;', strnum(1) );
eval(cmd) ;
liste_pix = pos + (ii') + jj'*M ;

for c=2:size(strnum, 2)

  cmd = sprintf('[ii, jj] = find(n%s) ;', strnum(c) );
  eval(cmd) ;
  %% limitation zone image verticale

  liste_pix = [liste_pix, pos + ii' + (jj'+(c-1)*4)*M ];

end%for
end%function
