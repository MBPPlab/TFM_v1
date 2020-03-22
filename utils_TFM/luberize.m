% MIT License
% 
% Copyright (c) 2018 HoffmanLab-Public
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

function [newtracks] = luberize(tracks)

% reassigns the unique ID# to 0,1,2,3...
% /presort will sort on ID# first, then reassign
% start will begin with that ID#

% function returns a new track array

ndat=length(tracks(1,:));
% if (keyword_set(presort)) then begin
%     newtracks=tracks(*,sort(tracks(ndat,*)))
% endif else begin
%     newtracks=tracks
% endelse

newtracks=tracks;

u=unq((newtracks(:,ndat))',[]);
ntracks=length(u);
u=[0,u];
for i=1:ntracks,  newtracks(u(i)+1:u(i+1),ndat) = i; end

% if (keyword_set(start)) then newtracks(ndat,*)=newtracks(ndat,*)+start

