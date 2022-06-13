clear all;
close all;
clc;

%% Load data and pre-processing
p = 5;
q = 3;
venc_integer = [q*(p-q);    p*(p-q);     p*q];                              % enforce integer to avoid possible by applying lcm to noninteger values
scaling = 21/venc_integer(1);
venc = scaling*venc_integer;

load('VIVO_FS_Data_Venc_21_35_52.5.mat')                                    % still has 25% partial fourier
% Coil compression or coil combining is optional, less coils will make PRoM
% faster. Here only coil compression is illustrated.

[Nx,Ny,Nt,Ne,Nc] = size(K);

K = reshape(K,[],Nc);
[U,S,V] = svd(K,'econ');
Nc = 8;                                                                     
K = reshape(U(:,1:Nc)*S(1:Nc,1:Nc),[Nx,Ny,Nt,Ne,Nc]);

I = K2I(K);                                                                 % k-space to image 

v1 = angle( sum(I(:,:,:,2,:).*conj(I(:,:,:,1,:)),5) )./pi*venc(1);
v2 = angle( sum(I(:,:,:,3,:).*conj(I(:,:,:,1,:)),5) )./pi*venc(2);
v3 = angle( sum(I(:,:,:,2,:).*conj(I(:,:,:,3,:)),5) )./pi*venc(3);

v123 = [v1(:),v2(:),v3(:)]';                                                % concatenation noisy possible wrapped velocity measurements

%% SDV
tic
disp('SDV Processing')
V_SDV = round((v3-v1)./(venc(1)*2))*venc(1)*2+v1;
toc

%% PRoM
Omega = 2*lcm(venc_integer(1),venc_integer(2));
KnownRange_integer = [0, 65/scaling, Omega-65/scaling, Omega];

Numofcand = 2;                                                              % number of most possible velocity candidates to use
V_PRoM_Plus = zeros(size(v1));
V_PRoM_Cand = zeros([numel(v1),Numofcand]);   
L_PRoM_Cand = zeros([numel(v1),Numofcand]);   

s1 = reshape(abs(I(:,:,:,1,:))/2+abs(I(:,:,:,3,:))/2,[],size(I,5));         % assuming first and third encoding have same noiseless 
s2 = reshape(abs(I(:,:,:,2,:)),[],size(I,5));
s3 = s1;

clear PRoM3
tic
disp('PRoM3 Processing')
parfor i = 1:numel(v1)
    COV = ...
        [   (sum(s2(i,:).^2)+sum(s1(i,:).^2))/(2*sum(s2(i,:).*s1(i,:))^2)           size(I,ndims(I))/(2*sum(s1(i,:).^2))                                size(I,ndims(I))/(2*sum(s2(i,:).^2));
            size(I,ndims(I))/(2*sum(s1(i,:).^2))                                    (sum(s3(i,:).^2)+sum(s1(i,:).^2))/(2*sum(s3(i,:).*s1(i,:))^2)       -size(I,ndims(I))/(2*sum(s3(i,:).^2));
            size(I,ndims(I))/(2*sum(s2(i,:).^2))     	                            -size(I,ndims(I))/(2*sum(s3(i,:).^2))   	                        (sum(s2(i,:).^2)+sum(s3(i,:).^2))/(2*sum(s2(i,:).*s3(i,:))^2)];
    
    % Project estimated covariance matrix to be positive semidefinite matrix
    [V,D] = eig(COV);
    COV = V*max(D,0)*inv(V);
    COV = diag(venc_integer/pi)*COV*diag(venc_integer/pi);
    COV_Inv = pinv(COV);

    [~,v_k,L] = PRoM3 (v123(:,i), venc,COV_Inv);
    V_PRoM_Cand(i,:) = v_k(1:Numofcand);
    L_PRoM_Cand(i,:) = L(1:Numofcand);
end
toc

V_PRoM_Cand = reshape(V_PRoM_Cand,size(v1,1),size(v1,2),size(v1,3),[]);
L_PRoM_Cand = reshape(L_PRoM_Cand,size(v1,1),size(v1,2),size(v1,3),[]);

I_SSOS = SSOS(I(:,:,:,1,:));

%% PRoM +
tic
disp('PRoM + postprocessing')
[x,y] = meshgrid(1:size(V_PRoM_Cand,2),1:size(V_PRoM_Cand,1));              % space position
parfor frame = 1: size(I_SSOS,3)
    Mask = logical((I_SSOS(:,:,frame) >= 0.3*max(I_SSOS(:,:,frame),[],'all')));

    % Erosion and dialiation to smooth the mask
    Mask = imdilate(imerode(Mask,strel("disk",1)),strel("disk",4));
        
    z = V_PRoM_Cand(:,:,frame,1);                                           % current choice of velocity estimate
    z_Cand = squeeze(V_PRoM_Cand(:,:,frame,:));                             % all velocity estiamtes for current frame                         
    z_old = zeros(size(z_Cand));
    
    Lambda = 1;                                                             % regularization weight
    
    counts = 0;
    while sum(abs(z-z_old)) ~=0
        counts = counts + 1;
        z_old = z;
        sf = fit([x(Mask),y(Mask)],z(Mask),'loess','Span',0.03);            % other curve fitting choices can be chosen
        [~,indx] = min(Lambda*(z_Cand - sf(x,y)).^2+z_Cand,[],3,'linear');
        Tmp = z_Cand(indx);
        z(Mask) = Tmp(Mask);
        disp(['     Frame ', num2str(frame),', number of iterations for post precessing is ', num2str(counts)])
        
    end
    
    V_PRoM_Plus(:,:,frame) = z;
end
toc

V_PRoM = V_PRoM_Cand(:,:,:,1);

%% ODV
V_ODV = zeros(size(v1));
tic
parfor i =1:numel(v1)
    V_ODV(i) = ODV(v1(i),v2(i),venc(1)*2/q,q,p);
end
toc


%% Analysis
Comparison = cat(2,V_SDV, V_ODV, V_PRoM,V_PRoM_Plus);
Encoding_Mag = cat(2,SSOS(I(:,:,:,1,:)),SSOS(I(:,:,:,3,:)),SSOS(I(:,:,:,2,:)));

figure;
imshow(Encoding_Mag(:,:,3),[]);
title('Magnitude images of all encodings for frame 3')

figure;
imshow(cat(2, V_SDV(80:125,25:65, 3),V_ODV(80:125,25:65, 3),V_PRoM(80:125,25:65, 3),V_PRoM_Plus(80:125,25:65, 3)),[]);
title('Velocity estimate of SDV, ODV, PRoM and PRoM+ for frame 3')

figure;
imshow(cat(2, V_SDV(:,:,3),V_ODV(:,:,3), V_PRoM(:,:,3), V_PRoM_Plus(:,:,3)),[]);
rectangle('Position',[25, 80, 40, 45],'LineWidth',4,'EdgeColor','#D95319','LineStyle','-.')

% If image processing toolbox is installed, you can view a the 3D volume
figure; orthosliceViewer(Comparison)
figure; orthosliceViewer(Encoding_Mag)


%% Functions
function I = K2I(Kdata)                                                                                                 % k-space to image domain
I = sqrt(size(Kdata,1)*size(Kdata,2))*fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(Kdata,1),2),[],1),[],2),1),2);
end

function Kdata = I2K(I)                                                                                                 % image domain to k-sapce
Kdata = 1/sqrt(size(I,1)*size(I,2))*fftshift(fftshift(fft(fft(ifftshift(ifftshift(I,2),1),[],2),[],1),2),1);
end


function I_SSOS = SSOS(I)                                                                                               % image to square root of sum of square
I_SSOS = sum(abs(I).^2,ndims(I)).^0.5;
end

function snr = SNR(X,Ref)                                                                                               % calculate the SNR
snr = -20*log10(norm(X(:)-Ref(:))/norm(Ref(:)));
end