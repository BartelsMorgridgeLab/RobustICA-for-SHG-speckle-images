%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Random Fields for the screens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Random potentials
V1=sig1*randpot(N,kmax1,Delta_x); %random potential for screen 1
Phase1=exp(2*pi*1i*V1); %this is the random screen 1, this is FIXED

V2=sig2*randpot(N,kmax2,Delta_x); %random potential for screen 2
Phase2=exp(2*pi*1i*V2); %this is the random screen 2, this is FIXED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Illuminations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Illum=cell(Nillum,1);
IllumData=cell(3,1);

parfor j=1:Nillum %loop over realizations of the phase on the SLM

    jp=j;
    fprintf('Forward Field calculation illumination %d out of %d \n',jp,Nillum)

    V=sigSLM*randpot(N,kmaxSLM,Delta_x); %random potential for SLM
    PhaseSLM=exp(2*pi*1i*V); %field on SLM

    EF=ForwardField(Phase1,Phase2,PhaseSLM,padding_f,k0,L0,L1,R,Xx,Xy,Ux,Uy); %forward field from SLM to sample plane
    Illum{j}=EF((-widthSample+N/2):(widthSample+N/2),(-widthSample+N/2):(widthSample+N/2)); %store the cropped field for each illumination
    
end

IllumData{1}=Illum; %store field
IllumData{2}=Phase1; %store phase screen 1
IllumData{3}=Phase2; %store phase screen 2

save IllumData IllumData -v7.3
