clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    PARAMETERS FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% What the main file does

GeneIllum="yes"; %set to 'yes' for generating illuminations
GeneReflec="ye"; %set to 'yes' for generating reflection matrix
Sep="ye"; %set to 'yes' for separation
Deconv="ye"; %set to 'yes' for deconvolution
PlotFig="ye"; %set to 'yes' for plotting image

%%%% Geometry

lambda=1; %wavelength
f=500*lambda; %focal
k0=2*pi/lambda; %wavenumber
NA=0.75; %numerical aperture
L0=0.8*f; %distance of the 1st random screen from the sample plane
L1=f/10; %distance of the 2nd random screen from the first screen
L=L0+L1; %distance from sample to second random screen
window=500*lambda; %size of the domain
R=NA/lambda; %radius of input pupil at the SLM, resolution is then the usual lambda/2NA
width=100; %1/2 width of the cropping window in the CCD plane 
Nc=2*width+2; %width of the cropping window in the CCD plane
widthSample=100; %1/2 width of the cropping window in the sample plane 

%%%% Scatterers

Nillum=1000; %number of illuminations to generate the reflection matrix
NillumSVD=Nillum; %number of illuminations used in the SVD (should be <= Nillum)
Startillum=1; %index of 1st illumination used in SVD

ScatType="RandScat";%"RandScat";%"2scat";"cross"; "circle", "2groups", "smiley" scatterer type

Nscat=10; %total number of scatterers
NscatSVD=10; %number of singular values kept after SVD 
%Nscat=284; %number of scatterers for smiley
%NscatSVD=66; %number of singular values kept after SVD for smiley
%NL=24;
%Nscat=4*NL+1; %number of scatterers for cross
%NscatSVD=44; %number of scatterers for cross

shift=15; %shift btw scatterers for cross, 2scat, 2 groups, or radius for circle
windowScat=25; %1/2 window size for placing random scatterers

rho_amp=0.; %amplitude susceptiblity fluctuations
rho_phase=0.; %phase susceptiblity fluctuations

NSVDbeg=1;% index of first singular value to be accounted in deconvolution


%%%% Random Fields

losc1=2*pi*lambda*2; %minimal wavelength length of the random screen 1
losc2=2*pi*lambda*2; %minimal wavelength length of the random screen 2
loscSLM=2*pi*lambda*2; %minimal wavelength length of the random screen SLM
kmax1=2*pi/losc1; %maximal wavenumber for the random field 1;
kmax2=2*pi/losc2; %maximal wavenumber for the random field 2;
kmaxSLM=2*pi/loscSLM; %maximal wavenumber for the random field SLM;
sig1=1/pi; %std random potential 1 
sig2=1/pi; %std random potential 2
sigSLM=1/pi; %std random potential SLM 

%%%% Numerics stuff

N=2^ceil(log2(window*5*NA/lambda)); % # of discretization points of the window
padding_f=0; %padding forward
padding_b=0; %padding backward
Npad=N*2^padding_b; %# of points when padding backward
Ncor=20; %number of points for computing the correlation

%% ICA parameters

%arguments = {'prewhi', false,'tol',1e-19, 'maxiter',15000}; %parameters for the ICA
arguments = {'prewhi', false,'tol',1e-19, 'maxiter',15000,'deftype', 'regression'}; %parameters for the ICA

%% Deconvolution parameters

mu = 2;
opts.rho_r   = 1; 
opts.rho_o   = 1;
opts.beta    = [1.0 1.0 0];
opts.gamma  = 2;
opts.print   = false;
opts.alpha   = 0.1;
opts.tol   = 1e-6;
opts.method  = 'l2';
opts.max_itr  = 1000;
opts.thres=0.0; %denoising parameter

save param

Main