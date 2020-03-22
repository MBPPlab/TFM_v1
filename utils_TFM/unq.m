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

%THIS CODE IS PLANNED FOR ARRAY BEING A ROW VECTOR

function [ret] = unq(array,idx)

s = size(array);
if s, else, warning('array must be an array'), end    %warning if s is empty
if idx
    q = array(idx);
    qshift = circshift(q,[0,-1]);
    indices = find(q~=qshift);
    if indices, ret = idx(indices);, else, ret = length(q);, end
else
    array=array;
    arrayshift = circshift(array,[0,-1]);
    indices = find(array~=arrayshift);
    if indices, ret = indices;, else, ret = length(array);, end
end