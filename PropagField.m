function [outputArg1] = PropagField(Phase,padding,k0,L,R,E,Xx,Xy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function propagates the input field E multiplied by a transfer
%function Phase over a distance L using paraxial (i.e. Schrodinger) propagator
%
%Phase: transfer function
%padding: # of points multiplied by 2^padding for the inverse Fourier transform 
%E: Input field
%k0: wavenumber
%L: total distance traveled
%R: circular aperture size
%Xx,Xy spatial grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[N,N]=size(Phase); %N should be even
dx=Xx(1,2)-Xx(1,1); %spatial grid stepsize
k=((-N/2+1):(N/2))*2*pi/(dx*N); %grid in Fourier space associated with the spatial grid Xx,Xy.
[Kx,Ky]=meshgrid(k,k);

fact=1i*exp(-1i*0.5*L*(Kx.^2+Ky.^2)/k0);

outputArg1=(ifft2(ifftshift(fact).*fft2(ifftshift(Phase.*E.*((Xx.^2+Xy.^2)<=R^2))),N*2^padding,N*2^padding));
outputArg1=fftshift(outputArg1);

%ifftshift puts the origin (0,0) in the first quadrant since fft2 wants (0,0) there
%fftshift moves the origin (0,0) from the first quadrant to the center

end