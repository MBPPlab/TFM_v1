%%  ****************************TFM_V1*************************************
%
%  Copyright 2014 by Martial Balland and Irene Wang (Univ. Grenoble Alpes,
%  LIPhy, F-38000 Grenoble, France. 2 CNRS, LIPhy, F-38000 Grenoble) 
%  Email: <martial.balland@ujf-grenoble.fr> 
%  
%  Changes: 
%  Copyright 2019 by Paolo Pierobon (Institut Curie, PSL Research
%  University, INSERM U932, 26 rue d’Ulm, 75248 Paris, Cedex 05, France.)
%  Email: <pierobon@curie.fr>
%  
%  Licensed under GNU General Public License 3.0 or later. 
%  Some rights reserved. See COPYING, AUTHORS.
%  @license GPL-3.0+ <http://spdx.org/licenses/GPL-3.0+>
% 
% The functions luberize, trackmem and unc are Copyright (c) 2018 HoffmanLab-Public
% under MIT Licence
%
% 
% If this software is use to generate analysis for publications please cite
% the following articles: 
%
% [1] Mandal K., Wang I., Vitiello E., Orellana L. A. C. and Balland
% M. (2014). Cell dipole behaviour revealed by ECM sub-cellular geometry.
% Nat. Commun. 5, 5749.
%
% [2] Kumari A., Pineau J., Sáez P. J., Maurin M., Lankar D., Roman
% M. S., Voituriez R., Hennig K., Boura V. F., Karlsson M., Balland M.,
% Lennon-Duménil A.-M. and Pierobon P.  (2019). Actomyosin-driven force
% patterning controls endocytosis at the immune synapse. Nat. Commun. 10,
% 2870.
%
% [3] Kumari A., Pineau J., Lennon-Duménil A.-M., Balland M. and Pierobon
% P. (2020), Traction force microscopy to study B lymphocyte activation,
% Journal of Visualized Experiments
%
% [4] Bertaux N., Marguet D. and Serge A., Dynamic multiple-target tracing
% to probe spatiotemporal cartography of cell membranes, Nat. Methods,
% vol. 5, no. 8, pp. 687–694, 2008.
%
%% *************************************************************************
% function map_max = all_max_2d(img)
%
% EN/ provides a binary map of all local max


function map_max = all_max_2d(img)
img(isnan(img))=0;
[N,M] = size(img) ; 
img0 = img(2:N-1, 2:M-1) ;
map_max = zeros(N,M) ;
map_max(2:N-1, 2:M-1)= img(1:N-2, 2:M-1) < img0 & ...
            img(3:N  , 2:M-1) < img0 & ...
            img(2:N-1, 1:M-2) < img0 & ...
            img(2:N-1, 3:M  ) < img0 & ...
            img(1:N-2, 1:M-2) < img0 & ...
            img(3:N  , 3:M  ) < img0 & ...
            img(3:N  , 1:M-2) < img0 & ...
            img(1:N-2, 3:M  ) < img0;
end