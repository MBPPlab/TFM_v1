%%  ****************************TFM_V1*************************************
%
%  Copyright 2014 by Martial Balland and Irene Wang (Univ. Grenoble Alpes,
%  LIPhy, F-38000 Grenoble, France. 2 CNRS, LIPhy, F-38000 Grenoble) 
%  Email: <martial.balland@ujf-grenoble.fr> 
%  
%  Changes: 
%  Copyright 2019 by Paolo Pierobon (Institut Curie, PSL Research
%  University, INSERM U932, 26 rue d'Ulm, 75248 Paris, Cedex 05, France.)
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

%% funtion to launch dialog for analysis
function param=TFM_dialog

prompt = {  'Pixel size (m):' ...
            'Time interval (s):'...
            'Young modulus (Pa):'...
            'Poisson ratio:'...
            'Regularization:'...
            'Window size (pixel):'...
            'Reference frame (0 to select file):'};
dlgtitle = 'Input parameters';
dims = [1 40];
definput  = {'0.106e-6', '5', '500', '0.5', '5e-19', '4', '1'};
answer = inputdlg(prompt,dlgtitle,dims,definput)

param.pix=str2num(answer{1});
param.dt=str2num(answer{2});
param.E=str2num(answer{3});
param.nu=str2num(answer{4});
param.alphadef=str2num(answer{5});
param.interval=double(uint16(str2num(answer{6})));
param.reference=double(uint16(str2num(answer{7})));

answer = questdlg('Do you need movies as output (slow)?', ...
    'Output options',...
	'Yes', 'No','No');
% Handle response
if strcmp(answer, 'Yes')
    param.movie=1;
else
    param.movie=0;
end

end