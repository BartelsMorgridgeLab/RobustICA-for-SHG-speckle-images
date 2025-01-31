function V=RandPot(N,kmax,dx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generates a mean-zero stationary random field whose
% power spectrum is the indicator function of the disk of radius kmax.
% The minimal wavelength length is 2pi/kmax and the variance 1 when kmax <=pi/dx. 
% The correlation length is 1/kmax.
%
%N: # of points along one direction
%kmax: maximal wavenumber
%dx: stepsize for spatial grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=((-N/2+1):(N/2))*2*pi/(dx*N); %wavenumber grid
[Kx,Ky]=meshgrid(k,k); %wavenumber grid

C=(-0.5+rand(N)).*exp(1i*2*pi*rand(N)); %random Fourier coefficients with 0 mean

V=C.*((Kx.^2+Ky.^2)<=kmax^2)*sqrt(24*N^2*(2*pi)^2/(pi*kmax^2*dx^2)); %Filter and normalize
V=real(fftshift(ifft2(ifftshift(V)))); %inverse FT and take real part



