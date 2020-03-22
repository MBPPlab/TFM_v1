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
% The functions luberize and trackmem are Copyright (c) 2018 HoffmanLab-Public
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

%% Function to generate movies of stress and displacement maps

function [frame1, frame2, frame3, frame4]=movie_quiverTFM(pathname, fileIn)
if nargin<1
[fileIn, pathname]=uigetfile('*movie.tif');
end
fileDispl=fullfile(pathname,[fileIn(1:end-10) '_displ.mat']);
fileForce=fullfile(pathname,[fileIn(1:end-10) '_force.mat']);

m=importTIF(fullfile(pathname,fileIn));
load(fileDispl);
load(fileForce);
%traj=TM_convert_log_mlab(d);
%frame=TM_PTV_dataconv(traj);

% axis ij;
% axis off;
% % axis tight
% set(gca,'nextplot','replacechildren');
% f=figure('Renderer','zbuffer');

% imagesc(imadjust(m(:,:,1)));
% colormap gray
% hold on
% plot(displ(1).positionX, displ(1).positionY,'ro');
% axis ij tight manual;
% axis off;
%set(gca,'nextplot','replacechildren');
%f=figure('Renderer','zbuffer');

% generate the limits for the full frame
bw=bwconncomp(displ(end).BWmask);
rp=regionprops(bw,'all');
s=12;
xl1=[rp.BoundingBox(1)-s,rp.BoundingBox(1)+rp.BoundingBox(3)+s];
yl1=[rp.BoundingBox(2)-s,rp.BoundingBox(2)+rp.BoundingBox(4)+s];

% generate the limits for the cropped frame
mask=logical(displ(end).BWmask(displ(end).gridY(:,1),displ(end).gridX(1,:)));
bw=bwconncomp(mask);
rp=regionprops(bw,'all');
s=3;
xl2=[rp.BoundingBox(1)-s,rp.BoundingBox(1)+rp.BoundingBox(3)+s];
yl2=[rp.BoundingBox(2)-s,rp.BoundingBox(2)+rp.BoundingBox(4)+s];

% compute scale for right rescaling
% displacement
sd=max(max((sqrt([displ.deplacementX].^2+[displ.deplacementY].^2))));
% force
sf=max(max((sqrt([force.TractionX].^2+[force.TractionY].^2))));
% energy
for i=1:length(force)
    u(:,:,i)=force(i).TractionX.*displ(i).deplacementX+force(i).TractionY.*displ(i).deplacementY;
end
su=max(u(:));
% number frame
n=length(force);


%% beads quiver
f=figure('Renderer','zbuffer');
axis ij; axis off; set(gca,'nextplot','replacechildren');
for i=1:n
    imagesc(imadjust(m(:,:,i)));
    colormap gray
    hold on
    plot(displ(i).positionX, displ(i).positionY,'ro');
    quiver(displ(i).positionX, displ(i).positionY, displ(i).deplaceX, displ(i).deplaceY,0,'k')
    B = bwboundaries(displ(i).BWmask);
    plot(B{1}(:,2), B{1}(:,1),'w','LineWidth',1);
    xlim(xl1); ylim(yl1);
    axis off;
    drawnow;
    pause(0.1);
    frame1(i)=getframe(gcf);
end
close all

%% displacement+quiver

f=figure('Renderer','zbuffer');
axis ij; axis off; set(gca,'nextplot','replacechildren');
for i=1:n
    colormap jet
    hold on
    imagesc(sqrt(displ(i).deplacementX.^2+displ(i).deplacementY.^2),[0 sd]);
    quiver(1/sd*displ(i).deplacementX, 1/sd*displ(i).deplacementY,0,'k');
    mask=logical(displ(i).BWmask(displ(i).gridY(:,1),displ(i).gridX(1,:)));
    B = bwboundaries(mask);
    plot(B{1}(:,2), B{1}(:,1),'w','LineWidth',1);
    xlim(xl2); ylim(yl2);
    axis off;
    drawnow;
    pause(0.1);
    frame2(i)=getframe(gcf);
end
close all;

%% force+quiver
f=figure('Renderer','zbuffer');
axis ij; axis off; set(gca,'nextplot','replacechildren');
for i=1:n
    colormap jet
    hold on
    imagesc(sqrt(force(i).TractionX.^2+force(i).TractionY.^2),[0 sf]);
    quiver(1/sf*force(i).TractionX, 1/sf*force(i).TractionY,0,'k');
    mask=logical(displ(i).BWmask(displ(i).gridY(:,1),displ(i).gridX(1,:)));
    B = bwboundaries(mask);
    plot(B{1}(:,2), B{1}(:,1),'w','LineWidth',1);
    xlim(xl2); ylim(yl2);
    axis off;
    drawnow;
    set(gca,'nextplot','replacechildren');
    pause(0.1);
    frame3(i)=getframe(gcf);
end
close all;

%% energy
f=figure('Renderer','zbuffer');
axis ij; axis off; set(gca,'nextplot','replacechildren');
for i=1:n
    colormap jet
    hold on
    imagesc(u(:,:,i),[0 su]);
    mask=logical(displ(i).BWmask(displ(i).gridY(:,1),displ(i).gridX(1,:)));
    B = bwboundaries(mask);
    plot(B{1}(:,2), B{1}(:,1),'w','LineWidth',1);
    xlim(xl2); ylim(yl2);
    axis off;
    drawnow;
    set(gca,'nextplot','replacechildren');
    pause(0.1);
    frame4(i)=getframe(gcf);
end
close all;
%
%% save
try
    transfermovie(frame1, fullfile(pathname,[fileIn(1:end-7) '_movie_bead.avi']))
    transfermovie(frame2, fullfile(pathname,[fileIn(1:end-7) '_movie_displ.avi']))
    transfermovie(frame3, fullfile(pathname,[fileIn(1:end-7) '_movie_force.avi']))
    transfermovie(frame4, fullfile(pathname,[fileIn(1:end-7) '_movie_energy.avi']))
catch
    disp('could not generate the movie');
end

end



function transfermovie(frame, fileOut)
writerObj = VideoWriter(fileOut, 'Uncompressed AVI');
open(writerObj);
for i=1:length(frame)
    sz1(i)=size(frame(i).cdata,1);
    sz2(i)=size(frame(i).cdata,2);
end
msz1=min(sz1);
msz2=min(sz2);

for i=1:length(frame)
    frame(i).cdata=frame(i).cdata(1:msz1,1:msz2,:);
    writeVideo(writerObj,frame(i));
end
close(writerObj);
end

function  mask=importTIF(filename);
    nFrames=length(imfinfo(filename));
    for i=1:nFrames
        mask(:,:,i)=imread(filename,i);
    end
end
%
%
