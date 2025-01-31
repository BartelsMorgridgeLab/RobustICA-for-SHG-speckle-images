function [outputArg2] = PropagField2F(Phase,padding,k0,L,R,E,Xx,Xy,Ux,Uy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the backward field at the CCD knowing the field at a plane at a distance L 
%from the sample plane assuming infinite aperture in the microscope lens
%
%Phase: transfer function of the screen at L 
%padding: # of points multiplied by 2^padding for the Fourier transform from L to CCD
%k0:wavenumber
%L: distance from sample plane
%R: circular aperture size at plane L
%E: Input field on the plane at L
%Xx, Xy spatial grid
%Ux, Uy spatial wavenumbers grid on the CCD, i.e. u=y/(lambda f) for y a pixel on the CCD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,N]=size(Phase);
lambda=2*pi/k0; %wavelength

outputArg2=fft2(ifftshift(Phase.*E.*((Xx.^2+Xy.^2)<=R^2)),N*2^padding,N*2^padding);
outputArg2=exp(1i*pi*lambda*L*(Ux.^2+Uy.^2)).*fftshift(outputArg2);

%ifftshift puts the origin (0,0) in the first quadrant since fft2 wants (0,0) there
%fftshift moves the origin (0,0) from the first quadrant to the center

end