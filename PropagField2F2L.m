function [outputArg1] = PropagField2F2L(k0,L,R,E,Ux,Uy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the forward field from the SLM through the
%microscope (with infinite aperture) to the plane at L from the sample
%plane
%
%k0: wavenumber
%L: distance from sample plane
%R: input pupil radius at SLM
%E: Input field at the SLM
%Ux, Uy spatial wavenumber grid on the SLM, i.e. u=y/(lambda f) for y a pixel on the SLM. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda=2*pi/k0;
fact=exp(1i*pi*lambda*L*(Ux.^2+Uy.^2));

outputArg1=fft2(ifftshift(fact.*E.*((Ux.^2+Uy.^2)<=R^2)));
outputArg1=fftshift(outputArg1);

end