function EL = ForwardField(Phase1,Phase2,PhaseSLM,padding,k0,L0,L1,R,Xx,Xy,Ux,Uy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function computes the forward field from the CCD to the sample plane
%with random screens at distances L0 and L=L0+L1 from the sample plane
%
%Phase1: transfer function screen 1
%Phase2: transfer function screen 2
%PhaseSLM: field on the SLM 
%padding: # of points multiplied by 2^padding for the Fourier transform from SLM to L
%k0:wavenumber
%L0: distance from the sample plane to the 1st screen
%L1: distance from the 1st screen to the 2nd screen
%R: radius of the input pupil at the SLM
%Xx, Xy spatial grid
%Ux, Uy spatial wavenumbers grid on the SLM, i.e. u=y/(lambda f) for y a pixel on the SLM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=L0+L1;
RM=400*max(max(Xx)); %set large to model infinite aperture at screens 1 and 2`

EL=PropagField2F2L(k0,L,R,PhaseSLM,Ux,Uy); %field from SLM through microscope to plane at L
EL=PropagField(Phase2,0,k0,L1,RM,EL,Xx,Xy); %field from L to L0
EL=PropagField(Phase1,padding,k0,L0,RM,EL,Xx,Xy); %field from L0 to sample plane

end