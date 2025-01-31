function EB = BackwardField(Phase1,Phase2,padding,k0,L0,L1,R,Xx,Xy,Ux,Uy,E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the backward field from the sample plane to the CCD
%with random screens at distances L0 and L=L0+L1 from the sample plane
%
%Phase1: transfer function screen 1
%Phase2: transfer function screen 2
%padding: # of points multiplied by 2^padding for the Fourier transform from L to CCD
%k0: wavenumber
%L0: distance from the sample plane to the 1st screen
%L1: distance from the 1st screen to the 2nd screen
%R: circular aperture size on the screens and CCD if needed
%Xx,Xy spatial grid
%Ux, Uy spatial wavenumber grid on the CCD
%E: Input field on the sample plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=L0+L1; %distance btw the sample plane and the 2nd screen
RM=4000*max(max(Xx)); %aperture set large to model infinite aperture
[N,N]=size(Xx);

EB=PropagField(ones(N),0,k0,L0,RM,E,Xx,Xy); %field from sample to screen at L0, infinite aperture, transfer function equal to 1
EB=PropagField(Phase1,0,k0,L1,RM,EB,Xx,Xy); %field from random screen at L0 with transfer Phase1 
% and infinite aperture propagated over distance L1

EB=PropagField2F(Phase2,padding,k0,L,RM,EB,Xx,Xy,Ux,Uy); %field from plane L with infinite aperture multiplied 
%by transfer function Phase2 to CCD

end