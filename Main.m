addpath('/Users/pinaud/WORK/NUMERICS/Randy/RobustICA');
addpath('/tempest/a/accounts/pinaud/WORK/NUMERICS/Randy/RobustICA');
addpath('/tempest/a/accounts/pinaud/WORK/NUMERICS/SpeckleICA/DATA');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Delta_x=window/N; %stepsize
Delta_x_padded=window/Npad; %padded stepsize

x=Delta_x*((-N/2+1):(N/2)); %main grid
[Xx,Xy]=meshgrid(x,x); %main grid

x_pad=Delta_x_padded*((-Npad/2+1):(Npad/2)); %padded main grid
[Xx_pad,Xy_pad]=meshgrid(x_pad,x_pad); %padded main grid

Ux=Xx/(Delta_x*window); %spatial frequencies on the CCD
Uy=Xy/(Delta_x*window); %spatial frequencies on the CCD

Ux_pad=(Xx_pad)/(Delta_x_padded*window); %spatial frequencies padded grid on the CCD
Uy_pad=(Xy_pad)/(Delta_x_padded*window); %spatial frequencies padded grid on the CCD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Generate Illuminations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GeneIllum=='yes'
    disp('Generate Illuminations')
    GenerateIlluminations
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Generate Reflection Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if GeneReflec=='yes'
    disp('Generate Reflection Matrix')
    GenerateReflectionMatrix
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SVD+source separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Sep=='yes'
    disp('Separation')
    Separation
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Deconvolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Deconv=='yes'
    disp('Deconvolution')

    load SepData
    NscatSVD=SepData{3};

    for k=NSVDbeg:(NscatSVD+NSVDbeg-1) % with NSVDbeg we can choose a range of singular values
        EBest{k-NSVDbeg+1}=reshape(SepData{2}(:,k),Nc,Nc); %ICA separated field
        EBestU{k-NSVDbeg+1}=reshape(SepData{1}(:,k),Nc,Nc); %SVD separated field
        EBest{k-NSVDbeg+1}=EBest{k-NSVDbeg+1}/max(max(abs(EBest{k-NSVDbeg+1}))); %set maximum to 1 
        EBestU{k-NSVDbeg+1}=EBestU{k-NSVDbeg+1}/max(max(abs(EBestU{k-NSVDbeg+1}))); %set maximum to 1
    end

    [Global_image2]=Deconvolution(NscatSVD,Nc,EBest,mu,opts); %image with ICA separated field, deconvolution with UCSD code
    %[Global_image2]=DeconvolutionGelma(NscatSVD,Nc,EBest,mu,opts); %deconvolution with Gelma
    %[Global_image2]=DeconvolutionFista(NscatSVD,Nc,EBest,mu,opts); %deconvolution with Fista
    save Global_image2 Global_image2

    %[Global_image2U]=Deconvolution(NscatSVD,Nc,EBestU,mu,opts); %image with SVD separated field
    %save Global_image2U Global_image2U

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Plotting image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if PlotFig=='yes'
    PlotImage
end
    