% function sig_blk = sigij_blink(traj, t)
%
% EN/ if blink, increase of the std
%
%
% FR/ si blink, augmentation de la std


function sig_blk = sigij_blink(traj, t)

global tab_param ;
global tab_var ;

%% offset de zone de blink nb_blink==(-offset)
if (tab_param(traj, 7*t+8) < 0)
  offset = tab_param(traj, 7*t+8)  ;
else
  offset = 0 ;
end %if

%% sigi = sigj
sig_blk = tab_var(traj, 7*t+3)  ;  %% sur i 

%% prise en compte du blink pour sig_i/j
if (offset < 0)
    sig_blk = sig_blk * sqrt(1-offset) ;
end %if

end%function
