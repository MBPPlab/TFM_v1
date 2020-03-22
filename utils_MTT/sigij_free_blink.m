% function sig_free_blk = sigij_free_blink(traj, t)
%
% EN/ if blink, increase of the std
% traj scalar or vector of trajectories
%
%
% FR/ si blink, augmentation de la std
% traj scalaire ou vecteur de trajectoires


function sig_free_blk = sigij_free_blink(traj, t)

global tab_param ;
global sig_free ;

traj = traj(:) ;
nb_traj = size(traj, 1) ;

%% offset
offset = tab_param(traj, 7*t+8) ;

%% offset de zone de blink nb_blink==(-offset)
%% offset = 0 si non_blink (>0)
%% idem sinon ;
nb_blink = - (offset .* (offset < 0))  ;

%% sigi = sigj
sig_free_blk = sig_free*ones(nb_traj, 1) ; 

%% prise en compte du blink pour sig_i/j
sig_free_blk = sig_free_blk .* sqrt(1+nb_blink) ;

end%function
