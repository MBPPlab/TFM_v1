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

%% Main TFM function; to be launched from TFM_dialog
function TFM_fun(param,myDir,myFile)    

    pathname=myDir;
    basename=myFile(1:end-10);
    %    disp(basename);
    BKname=myFile;
    
    % AKname=[basename '_AK.tif'];
    % BFname=[basename '_BF.tif'];
    Maskname=[basename '_Mask.tif'];
    
    
    
    % beads
    tempfilename=BKname;
    if tempfilename
        filename=tempfilename;
        Name=filename;
        info = imfinfo(fullfile(pathname,filename),'tif');
        Nb_frames=numel(info);
        I=imread(fullfile(pathname,filename),'tif');
        cl=class(I);
        stressed=zeros([size(I) Nb_frames], cl);
        stressed(:,:,1)=I;
        for frame=2:Nb_frames
            stressed(:,:,frame)=imread(fullfile(pathname,filename),frame);
        end
    end
    
    if param.reference==0
        [refName, refPath]=uigetfile;
        nonstressed=imread(fullfile(refPath,refName));
    else
        nonstressed= stressed(:,:,param.reference);
    end
    cell=stressed;
    
    % mask
    tempfilename=Maskname;
    if tempfilename
        I=imread(fullfile(pathname,tempfilename),'tif');
        mask=zeros([size(I) Nb_frames], class(I));
        mask(:,:,1)=I;
        for frame=2:Nb_frames
            mask(:,:,frame)=imread(fullfile(pathname,tempfilename),frame);
        end
    end
    mask(mask~=0)=1;
    
    disp([basename ' loaded']);
    
    %%Parameters
    pix = param.pix;% taille du pixel (commonly we used SD4 at 60x px=0.106e-6)
    E=param.E;
    nu=param.nu; %Poisson ratio
    alphadef=param.alphadef;
    interval=param.interval;
    window=floor(min(size(I(:,:,1)))/interval-1)*interval;
    
    Energie=zeros(1,Nb_frames);%�nergie contractile
    Mu=zeros(1,Nb_frames);% moment contractile calul� dans l'espce de Fourier
    Moment_bis=zeros(1,Nb_frames);% moment contractile dans espace R�el
    Pmoy=zeros(1,Nb_frames);%contrainte moyenne
    Tmax=zeros(1,Nb_frames);%contrainte max
    Ftot=zeros(1,Nb_frames);
    Fvect=zeros(1,Nb_frames);
    surface=zeros(1,Nb_frames);
    nbeads=zeros(1,Nb_frames);
    noiseFtot=zeros(1,Nb_frames);
    noiseFvect=zeros(1,Nb_frames);
    noiseU=zeros(1,Nb_frames);
    
    
    %% Begin analysis
    pas=interval*pix;
    %Tracking parameters
    featsize=3; % partical radius in pixel, it is rather robust between 2 and 4 
    barrg=7; % giration ratio (useless in new version)
    barcc=0.25; % excentricity (useless in new version)
    IdivRg=0; %intensit� int�gr�e minimale pour detecter une particule
    masscut=270; %intensit� de coupure pour la d�tection des particules , important mais attention si trop grand on vire les particules sous le pattern
    maxd=4;%deplacement max d'une bille entre les 2 images
    threshMTT=15; 
    ndefl=1;
    
        
    
    [n,s,p_corr,m_corr,regx,regy]=registration_ssinterp_pattern(nonstressed,stressed(:,:,1),cell(:,:,1),mask(:,:,1));
    [dimy,dimx]=size(n);
    nbx=floor(dimx/window);
    nby=floor(dimy/window);
    [posx,posy]=meshgrid((0:nbx-1)*window+floor((dimx-(nbx-1)*window)/2)+1,(0:nby-1)*window+floor((dimy-(nby-1)*window)/2)+1);
    nodeX=floor(window*nbx/interval);
    nodeY=floor(window*nby/interval);
    [gridX,gridY]=meshgrid((0:nodeX-1)*interval+floor((dimx-nbx*window)/2)+1+interval/2,(0:nodeY-1)*interval+floor((dimy-nby*window)/2)+1+interval/2);
    
    
    param.window=window;
    param.featsize=featsize;% rayon de la particule en pixel, detect�e par un masque de taille featsize typiquement 2 ou 3
    param.barrg=barrg; % rayon de giration (reli� � la taille de la particule en pixel)
    param.barcc=barcc; % excentricit� (0 = rond 1 = pas rod du tout) typiquement entre 0.1 OU 0.2
    param.IdivRg=IdivRg;%intensit� int�gr�e minimale pour detecter une particule
    param.masscut=masscut; %intensit� de coupure pour la d�tection des particules , important mais attention si trop grand on vire les particules sous le pattern
    param.maxd=maxd;%deplacement max d'une bille entre les 2 images
    param.threshMTT=threshMTT;
    param.ndefl=ndefl;
    
    % replace "parfor" with "for" if you don't have the Parallel Computing
    % Toolbox
    parfor ff=1:Nb_frames
        
        
        % disp(['Analysis image ',num2str(ff)]);
        % [n,s,p_corr,regx,regy]=registration_ssinterp(nonstressed,stressed(:,:,ff),cell(:,:,ff));
        
        [n,s,p_corr,m_corr,regx,regy]=registration_ssinterp_pattern(nonstressed,stressed(:,:,ff),cell(:,:,ff),mask(:,:,ff));
        [dimy,dimx]=size(n);
        nbx=floor(dimx/window);
        nby=floor(dimy/window);
        [posx,posy]=meshgrid((0:nbx-1)*window+floor((dimx-(nbx-1)*window)/2)+1,(0:nby-1)*window+floor((dimy-(nby-1)*window)/2)+1);
        nodeX=floor(window*nbx/interval);
        nodeY=floor(window*nby/interval);
        [gridX,gridY]=meshgrid((0:nodeX-1)*interval+floor((dimx-nbx*window)/2)+1+interval/2,(0:nodeY-1)*interval+floor((dimy-nby*window)/2)+1+interval/2);
        BWmask=logical(m_corr);
        restex=regx-round(regx);
        restey=regy-round(regy);
        deplaceX=zeros(10000,1);
        deplaceY=zeros(10000,1);
        positionX=zeros(10000,1);
        positionY=zeros(10000,1);
        compteur=0;
        mauvaises_fenetres=0;
        mauvais_track=false;
        
        
        
        %     disp(['Nombre de fenètres: ',num2str(nby*nbx)]);
        for i=1:nbx
            for j=1:nby
                % crop windows
                n_win=n(posy(j,1)-window/2:posy(j,1)+window/2-1,posx(1,i)-window/2:posx(1,i)+window/2-1);
                s_win=s(posy(j,1)-window/2:posy(j,1)+window/2-1,posx(1,i)-window/2:posx(1,i)+window/2-1);
                
                % compute registration
                [offsetx,offsety,Cmax]=decalage(s_win,n_win);
                
                %faire une interpolation pour d�terminer n d�cal�
                xi=(1:window)-offsetx+posx(1,i)-window/2-1;
                yi=(1:window)-offsety+posy(j,1)-window/2-1;
                n_interp = interp2(double(n),xi',yi);
                %nouvelle image s_corr avec un petit d�placement
                %       s(posy(j,1)-window/2:posy(j,1)+window/2-1,posx(1,i)-window/2:posx(1,i)+window/2-1)=s_interp;
                if strcmp(cl,'uint16')
                    %n_win=uint16(n_interp);
                    n_win= uint16(65535/(max(n_interp(:))-min(n_interp(:)))*(n_interp-min(n_interp(:))));
                    s_win= uint16(65535/double(max(s_win(:))-min(s_win(:)))*double(s_win-min(s_win(:))));
                elseif strcmp(cl,'uint8')
                    n_win= uint8(255/(max(n_interp(:))-min(n_interp(:)))*(n_interp-min(n_interp(:))));
                    s_win= uint8(255/double(max(s_win(:))-min(s_win(:)))*double(s_win-min(s_win(:))));%n_win=uint8(n_interp);
                end
                
                
                %Tracking sur les petites fen�tres
                % nM = feature2D(n_win,1,featsize,masscut);
                % % replace the feature2D by the detection
                % function of MTT and adapting the parameters
                % r=detect2D(input, win, radius, thresh, ndefl);
                nM=detect2D(double(n_win), barrg, featsize, threshMTT, ndefl);
                
                if numel(nM)==1||isempty(nM)
                    figure(1)
                    imshow(n_win,[])
                    mauvais_track=true;
                else
                end
                
                % sM = feature2D(s_win,1,featsize,masscut);
                % %  replace the feature2D by the detection
                % function of MTT and adapting the parameters
                % r=detect2D(input, win, radius, thresh, ndefl);
                sM=detect2D(double(s_win), barrg, featsize, threshMTT, ndefl);
                
                if numel(sM)==1||isempty(sM)
                    figure(1)
                    imshow(s_win,[])
                    mauvais_track=true;
                else
                    
                end
                if not(mauvais_track)
                    Mtrack=[nM(:,1:2),ones(size(nM,1),1),ones(size(nM,1),1); sM(:,1:2),2*ones(size(sM,1),1),2*ones(size(sM,1),1)];
                    [lub] = trackmem(Mtrack,maxd,2,2,0);
                end
                
                %num�ro des images: null=1, stressed=2
                if lub==-1
                    figure(1)
                    subplot(1,2,1), imshow(s_win,[])
                    subplot(1,2,2), imshow(n_win,[])
                    mauvais_track=true;
                end
                if mauvais_track
                    mauvaises_fenetres=mauvaises_fenetres+1;
                    deplaceX(compteur+1:compteur+2)=offsetx+restex;
                    deplaceY(compteur+1:compteur+2)=offsety+restey;
                    positionX(compteur+1:compteur+2)=posx(1,i);
                    positionY(compteur+1:compteur+2)=posy(j,1);
                    compteur=compteur+1;
                else
                    nb_part=max(lub(:,5));
                    
                    deplaceX(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),1)-lub((lub(:,4)==1),1)+offsetx+restex;
                    deplaceY(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),2)-lub((lub(:,4)==1),2)+offsety+restey;
                    positionX(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),1)+(posx(1,i)-window/2)-1;
                    positionY(compteur+1:compteur+nb_part)=lub((lub(:,4)==2),2)+(posy(j,1)-window/2)-1;
                    compteur=compteur+nb_part;
                end
            end
        end
        %     disp([num2str(compteur-mauvaises_fenetres) ' features tracked']);
        %     disp([num2str(mauvaises_fenetres) ' mauvaises fenetres']);
        deplaceX(compteur+1:10000)=[];
        deplaceY(compteur+1:10000)=[];
        positionX(compteur+1:10000)=[];
        positionY(compteur+1:10000)=[];
        
        %% computation
        %% Interpolation code
        % nodeX=floor(window*nbx/interval);
        % nodeY=floor(window*nby/interval);
        % [gridX,gridY]=meshgrid((0:nodeX-1)*interval+floor((dimx-nbx*window)/2)+1+interval/2,(0:nodeY-1)*interval+floor((dimy-nby*window)/2)+1+interval/2);
        
        deplacementX=pix*griddata(positionX,positionY,deplaceX,gridX,gridY);
        deplacementY=pix*griddata(positionX,positionY,deplaceY,gridX,gridY);
        deplacementX(isnan(deplacementX))=0;
        deplacementY(isnan(deplacementY))=0;
        
        %%  Inversion of the problem
        
        [TractionX,TractionY,mu1,mu2,theta]=regularized(deplacementX,deplacementY,E,nu,pas,alphadef);
        Tmagn=sqrt(TractionX.^2+TractionY.^2);
        mu=mu1+mu2;
        
        %% Computation on mask
        
        BWmask=logical(mask(:,:,ff));
        cc=bwconncomp(BWmask);
        if cc.NumObjects>0
            rp=regionprops(cc,'MajorAxisLength');
            SE = strel('disk',double(int16(0.1*rp(1).MajorAxisLength))); %dilatation du masque
        else
            SE = strel('disk',20);
        end
        BWmask2=imdilate(BWmask,SE);
        BWgrid=logical(BWmask2(gridY(:,1),gridX(1,:)));
        Txcrop=TractionX.*BWgrid;
        Tycrop=TractionY.*BWgrid;
        Tmagncrop=sqrt(Txcrop.^2+Tycrop.^2);
        nbeads(ff)=inROI(positionX,positionY,BWmask2);
        numbeads= nbeads(ff);
        
        %% Memorize results
        displ(ff)=struct('deplacementX',{deplacementX},'deplacementY',{deplacementY}, 'BWmask',{BWmask},'gridX',{gridX},'gridY',{gridY},'positionX',{positionX},'positionY',{positionY},'deplaceX',{deplaceX},'deplaceY',{deplaceY},'interval',{interval},'pix',{pix},'regx',{regx},'regy',{regy},'numbeads',{numbeads});
        force(ff)=struct('TractionX',{TractionX},'TractionY',{TractionY},'gridX',{gridX},'gridY',{gridY},'alphadef',{alphadef});
        
        %% Global computation
        %calcul de l'énergie contractile
        u=Txcrop.*deplacementX+Tycrop.*deplacementY;  %produit scalaire traction par deplacement
        U=pas^2/2*sum(u(:));
        
        %Pression max sur la cellule
        TmagnMAX=max(Tmagncrop(:));
        
        %surface de la cellule
        Scell=bwarea(BWmask)*pix^2;%en m au carr�
        Sgrid=bwarea(BWgrid);
        
        %Somme des forces sous la cellule
        Ttot=sum(Tmagncrop(:));
        Ftot(ff)=Ttot*(pix*interval)^2;
        
        %Traction moyenne
        TmagnMOY=Ttot/Sgrid;% pression moyenne exerc�e par la cellule
        
        Fvect(ff)=sqrt((sum(Txcrop(:)))^2+(sum(Tycrop(:)))^2)*(pix*interval)^2;%sum magnitude forces
        Energie(ff)=U;%energie contractile
        Mu(ff)=mu*(pix*interval)^2;%moment contractile
        Mu1(ff)=mu1*(pix*interval)^2; %moment contractile first eig
        Mu2(ff)=mu2*(pix*interval)^2; %moment contractile second eig
        Theta(ff)=theta;%angle du moment
        Pmoy(ff)=TmagnMOY;%contrainte moyenne
        Tmax(ff)=TmagnMAX;%contrainte max
        surface(ff)=Scell;%surface cellule
        
        
        % if (nargin>2 && verb_opt) displayres(U, mu, TmagnMAX, Scell, Ftot, TmagnMOY, noiseU, pix, interval, ff); end;
        %                    if noiseanalysis
        %                         [noiseFtot(ff), noiseFvect(ff), noiseU(ff)]=computeNoise(noise(:,:,ff), deplacementX, deplacementY, TractionX, TractionY, gridX, gridY, pas);
        %                    end
        
        % outputFig(stressed(:,:,ff),BWmask(:,:,ff),displ,'dfd',ff);
        %                     catch
        %                         disp(['problem in ' num2str(ff)]);
        %                     end
    end
    
    %% outputs
    %   nbeads=[displ.numbeads];
    FinalTable=[(1:Nb_frames)',Energie',Mu', Ftot',Pmoy',Tmax', Fvect',surface', Theta', Mu1', Mu2', nbeads'];  % in abscence of noise analysis the noise associated vectors will be set to zero
    % pathname='/home/paolo/Desktop/test_TFMLSM';
    % pathname='/Volumes/LACIE SHARE/BC_TFM/FTM_Analysis_September2015/Analysed_2015/20150609_analysis/Concentrations';
    % pathname=myDir;
    save(fullfile(pathname,[basename '_results.mat']), 'FinalTable');
    save(fullfile(pathname,[basename '_force.mat']), 'force');
    save(fullfile(pathname,[basename '_displ.mat']), 'displ');
    save(fullfile(pathname,[basename '_param.mat']), 'param');
    disp('DONE!!');
    outputFigNonoise(nonstressed,mask(:,:,end),displ,fullfile(pathname,[basename '.fig']),1);
    % dlmwrite(fullfile(PathName,filename),[(1:Nb_frames)',Energie',Mu', Ftot',Pmoy',Tmax', Fvect',surface', DTheta'],'-append','delimiter','\t')
    figure
    plot([0:Nb_frames-1]*param.dt, Energie);
    xlabel('Time');
    ylabel('Energy (J)');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% additional functions
% this function is active when verbose option verb_opt is set to 1
function displayres(U, mu, TmagnMAX, Scell, Ftot, TmagnMOY, pix, interval, ff);
disp(['Contractile energy inside cell: ',num2str(U), ' J']);
disp(['Contractile moments trace: ',num2str(mu*(pix*interval)^2), ' N.m']);
disp(['Max traction: ',  num2str(TmagnMAX),' Pa']);
disp(['Cell surface: ',  num2str(Scell*1e12),' µm^2']);
disp(['Total forces: ',  num2str(Ftot(ff)),' N']);
disp(['Mean traction: ',  num2str(TmagnMOY),' Pa']);
end

% this function output the figure with the beads and save it
function outputFig(data,mask,noise,displ,fn,i);
figure;
% imagesc(double(data)+double(mask*256))
data=double(data);
imagesc(data,[min(data(:)) (min(data(:))+max(data(:)))/2]);
colormap gray
hold on
B = bwboundaries(mask);
plot(B{1}(:,2), B{1}(:,1),'r','LineWidth',1);
B = bwboundaries(noise);
plot(B{1}(:,2), B{1}(:,1),'b','LineWidth',1);
plot(displ(i).positionX,displ(i).positionY,'or');
saveas(gca, fn);
close;
end

% same as above WITHOUT NOISE region
function outputFigNonoise(data,mask,displ,fn,i);
figure;
% imagesc(double(data)+double(mask*256))
data=double(data);
imagesc(data,[min(data(:)) (min(data(:))+max(data(:)))/2]);
colormap gray
hold on
B = bwboundaries(mask);
plot(B{1}(:,2), B{1}(:,1),'r','LineWidth',1);
plot(displ(i).positionX,displ(i).positionY,'or');
saveas(gca, fn);
close;
end

% this function computes the properties of the noise IF required
function [noiseFtot, noiseFvect, noiseU]=computeNoise(noiseMask, deplacementX, deplacementY, TractionX, TractionY, gridX, gridY, pas);
BWmask=logical(noiseMask);
cc=bwconncomp(BWmask);
if cc.NumObjects>0
    rp=regionprops(cc,'MajorAxisLength');
    SE = strel('disk',double(int16(0.1*rp(1).MajorAxisLength))); %dilatation du masque
else
    SE = strel('disk',20);
end
BWmask2=imdilate(BWmask,SE);
BWgrid=logical(BWmask2(gridY(:,1),gridX(1,:)));
Txcrop=TractionX.*BWgrid;
Tycrop=TractionY.*BWgrid;
Tmagncrop=sqrt(Txcrop.^2+Tycrop.^2);
u=Txcrop.*deplacementX+Tycrop.*deplacementY;  %produit scalaire traction par deplacement
noiseU=pas^2/2*sum(u(:));
noiseFtot=sum(Tmagncrop(:))*(pas)^2;
noiseFvect=sqrt((sum(Txcrop(:)))^2+(sum(Tycrop(:)))^2)*pas^2;%sum magnitude forces
end


function r=detect2D(input, win, radius, thresh, ndefl);
s = detect_et_estime_part_1vue_deflt(input, win, radius, thresh, ndefl);
r=s; % recast in matrix structure given by feature2D
r(:,1)=s(:,3)+1; % x
r(:,2)=s(:,2)+1; % y
r(:,3)=s(:,3);
r(:,4)=s(:,6);
r(:,5)=0;
end

% this function counts the number of beads in the mask
function n=inROI(x,y,mask)
in=zeros(size(x));
for i=1:length(x)
    if (round(x(i))>0 & round(x(i))<size(mask,2) & round(y(i))>0 & round(y(i))<size(mask,1))
        in(i)=mask(round(y(i)),round(x(i)));
    else
        in(i)=0;
    end
end
n=sum(in);
end



function [Tractionx,Tractiony,mu1,mu2,theta]=regularized(deplacementx,deplacementy,E,nu,ech,alpha)
%ech=interval*pix;
[M,N]=size(deplacementx);
%padding with zeros
m=2;
while (2^m < M)||(2^m < N )
    m=m+1;
end
M2=2^m;
N2=M2;

uxt=fft2(deplacementx,M2,N2);
uyt=fft2(deplacementy,M2,N2);
%suppression de la translation
uxt(1,1)=0;
uyt(1,1)=0;

Kx=[(2*pi/ech)/N2*[0:N2/2],-(2*pi/ech)/N2*(N2-[N2/2+1:N2-1])];
Ky=[(2*pi/ech)/M2*[0:M2/2],-(2*pi/ech)/M2*(M2-[M2/2+1:M2-1])];
[kx, ky]=meshgrid(Kx,Ky);
k=sqrt(kx.^2+ky.^2);
Txt=zeros(M2,N2);
Tyt=zeros(M2,N2);

for i=1:M2
    for j=1:N2
        if (i==M2/2+1)||(j==N2/2+1)  %Fr�quences de Nyquist
            Gt=2*(1+nu)/(E*k(i,j)^3)*[(1-nu)*k(i,j)^2+nu*ky(i,j)^2  0; 0 (1-nu)*k(i,j)^2+nu*kx(i,j)^2] ;
            Tt=(Gt'*Gt+alpha*eye(2))^-1*Gt'*[uxt(i,j); uyt(i,j)];
            Txt(i,j)=Tt(1);
            Tyt(i,j)=Tt(2);
        elseif ~((i==1) && (j==1))
            Gt=2*(1+nu)/(E*k(i,j)^3)*[(1-nu)*k(i,j)^2+nu*ky(i,j)^2  -nu*kx(i,j)*ky(i,j); -nu*kx(i,j)*ky(i,j) (1-nu)*k(i,j)^2+nu*kx(i,j)^2] ;
            Tt=(Gt'*Gt+alpha*eye(2))^-1*Gt'*[uxt(i,j) ; uyt(i,j)];
            Txt(i,j)=Tt(1);
            Tyt(i,j)=Tt(2);
            
        end
    end
end


Tx=ifft2(Txt);
Ty=ifft2(Tyt);
Tractionx=real(Tx([1:M],[1:N]));
Tractiony=real(Ty([1:M],[1:N]));

%calcul des moments contractile(selon Butler)

% careful with i!!!!!
if ((~isempty(Tractionx)) && (~isempty(Tractiony)))
    Mxx=-1i*(Txt(1,2)+Txt(1,N2))/(2*kx(1,2));
    Myy=-1i*(Tyt(2,1)+Tyt(M2,1))/(2*ky(2,1));
    Mxy=-1i/2*((Tyt(1,2)+Tyt(1,N))/(2*kx(1,2))+(Txt(2,1)+Txt(M,1))/(2*ky(2,1)));
    
    %net contractile moment
    m=real([Mxx,Mxy; Mxy,Myy]);
    eigenv=eig(m);
    mu1=eigenv(1);
    mu2=eigenv(2);
    
    %angle de l'axe principal en deg
    theta=real(180/pi*atan(2*Mxy/(Mxx-Myy))/2);
    
else
    Mxx=nan;
    Myy=nan;
    Mxy=nan;
    mu1=nan;
    mu2=nan;
    theta=nan;
end
end

function [n_corr,s_corr,p_corr,t_corr,regx,regy]=registration_ssinterp_pattern(nonstressed,stressed,widefield,pattern)
%corrige le d�calge au pixel pr�s pour les 4 images: billes stress�es, non
%stress�es, phase et fluo pattern
%soustraction de background
background = imopen(stressed,strel('disk',12));
stressed2=imsubtract(stressed,background);
background = imopen(nonstressed,strel('disk',12));
nonstressed2=imsubtract(nonstressed,background);
%augmentation du contraste
n=contrast(nonstressed2,0.1,0.1);
s=contrast(stressed2,0.1,0.1);
[regx,regy]=decalage(n, s);

deltax = round(regx);
deltay = round(regy);

%correction du d�calage
if (deltax>0) && (deltay>0)
        n_corr=n(deltay+1:size(nonstressed,1),deltax+1:size(nonstressed,2));
        s_corr=s(1:size(stressed,1)-deltay,1:size(stressed,2)-deltax);
        p_corr=widefield(1:size(stressed,1)-deltay,1:size(stressed,2)-deltax);
        t_corr=pattern(1:size(stressed,1)-deltay,1:size(stressed,2)-deltax);
    elseif (deltax>0) && (deltay<=0)
        n_corr=n(1:size(nonstressed,1)-abs(deltay),deltax+1:size(nonstressed,2));
        s_corr=s(abs(deltay)+1:size(stressed,1),1:size(stressed,2)-deltax);
        p_corr=widefield(abs(deltay)+1:size(stressed,1),1:size(stressed,2)-deltax);
        t_corr=pattern(abs(deltay)+1:size(stressed,1),1:size(stressed,2)-deltax);
    elseif (deltax<=0) && (deltay>0)
        n_corr=n(deltay+1:size(nonstressed,1),1:size(nonstressed,2)-abs(deltax));
        s_corr=s(1:size(stressed,1)-deltay,abs(deltax)+1:size(stressed,2));
        p_corr=widefield(1:size(stressed,1)-deltay,abs(deltax)+1:size(stressed,2));
        t_corr=pattern(1:size(stressed,1)-deltay,abs(deltax)+1:size(stressed,2));
    else
        s_corr=s(abs(deltay)+1:size(stressed,1),abs(deltax)+1:size(stressed,2));
        n_corr=n(1:size(nonstressed,1)-abs(deltay),1:size(nonstressed,2)-abs(deltax));
         p_corr=widefield(abs(deltay)+1:size(stressed,1),abs(deltax)+1:size(stressed,2));
         t_corr=pattern(abs(deltay)+1:size(stressed,1),abs(deltax)+1:size(stressed,2));
end
end

function [deltax,deltay,Cmax]=decalage(im1,im2)
%Attention!! la taille des 2 images doivent �tre identiques
Nx=size(im1,2);
Ny=size(im1,1);
if rem(Nx,2)
    Nx=Nx-1;
    im1=im1(1:Ny,1:Nx);
    im2=im2(1:Ny,1:Nx);   
end
if rem(Ny,2)
     Ny=Ny-1;
     im1=im1(1:Ny,1:Nx);
    im2=im2(1:Ny,1:Nx);
end

A1=real(ifft2(fft2(im1) .* fft2(rot90(im1,2))));
A2=real(ifft2(fft2(im2) .* fft2(rot90(im2,2))));
M1=max(A1(:));
M2=max(A2(:));
C = real(ifft2(fft2(im1) .* fft2(rot90(im2,2))));
C2=fftshift(C);
%figure, imshow(C2,[])
[Vmax,imax]=max(C2(:));
Cmax=Vmax/sqrt(M1*M2);
[I,J]=ind2sub(size(C),imax);
if (I==size(C,1)) || (I==1)
    Isub=I;
else
%mesure d�calage pr�cision subpixel
    Isub=I+(log(C2(I-1,J))-log(C2(I+1,J)))/(2*(log(C2(I-1,J))+log(C2(I+1,J))-2*log(C2(I,J))));
end
if (J==size(C,2)) || (J==1)
    Jsub=J;
else
Jsub=J+(log(C2(I,J-1))-log(C2(I,J+1)))/(2*(log(C2(I,J-1))+log(C2(I,J+1))-2*log(C2(I,J))));
end
deltax=Jsub-Nx/2; 
deltay=Isub-Ny/2;
% si deltay>0 im2 d�cal�e vers le haut par rapport � im1
% si deltax>0 im2 d�cal�e vers la gauche par rapport � im1
end

function [Iout,Imax]=contrast(Iin,Fsat,Fzero)
%Iin is the input image
%Fsat in the percentage of saturated pixels
%Fzero is the percentage of pixels at zero intensity
cl=class(Iin);
Nbpix=numel(Iin);
Iin=double(intmax(cl))/double(max(Iin(:))-min(Iin(:)))*double(Iin-min(Iin(:)));
Imax=intmax(cl);
Imin=0;
if strcmp(cl,'uint16')
   X=Iin>Imax;
   while sum(X(:))<0.01*Fsat*Nbpix
       Imax=Imax-100;
       X=Iin>Imax;
   end
   Y=Iin<Imin;
   while sum(Y(:))<0.01*Fzero*Nbpix
       Imin=Imin+100;
       Y=Iin<Imin;
   end
Iout=uint16(double(intmax(cl))/double(Imax-Imin)*double(Iin-Imin));
elseif strcmp(cl,'uint8')
   X=Iin>Imax;
   while sum(X(:))<0.01*Fsat*Nbpix
       Imax=Imax-5;
       X=Iin>Imax;
   end
   Y=Iin<Imin;
   while sum(Y(:))<0.01*Fzero*Nbpix
       Imin=Imin+5;
       Y=Iin<Imin;
   end
Iout=uint8(double(intmax(cl))/double(Imax-Imin)*double(Iin-Imin));
end   
end

function  mask=importTIF(filename);
    nFrames=length(imfinfo(filename));
    for i=1:nFrames
        mask(:,:,i)=imread(filename,i);
    end
end
