%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Reflection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%load the data
load ReflData; %load the reflection matrix data

G=ReflData{1}; %load backward field
rho=ReflData{2}; %load susceptibilities, not always needed, see below
H=ReflData{3}; %load forward field

%%%%%Susceptility of the scatterers (here set a random)
rho=(1+rho_amp*(rand(Nscat,1)-0.5)).*exp(rho_phase.*2*pi*1i*rand(Nscat,1));
rho=diag(rho);

%%%%%Build the reflection matrix + SVD
Refl=G*rho*transpose(H(Startillum:(NillumSVD+Startillum-1),:)); %reflection matrix
[U,S,V] = svds(Refl,NscatSVD); %SVD of the reflection matrix keeping NscatSVD singular values

%%%%% Store the singular value before filtering
SS=diag(S);
save SS SS 

%%%%% Removing small singular values (useful when keeping all scatterers)
if min(SS)<0.01
    disp('Small singular value detected and removed');   
    I=find(SS<=0.01);
    NscatSVD=NscatSVD-length(I);
    U=U(:,1:NscatSVD);
    S=S(1:NscatSVD,1:NscatSVD);
    V=V(:,1:NscatSVD);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Source separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Using V
[X AA iter W]= robustica(V',arguments); %separation of the V
Gest=U*S*AA; %separation of the G from the H, gives the estimated separated backward fields

%%%% Using U
%[X HH iter2 WW]= robustica(U(:,1:NscatSVD)',arguments); %separation of the U
%Gest=X'; %the estimated separated backward fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Test the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NscatSVD==Nscat
    [EEH]=errorICA(H(Startillum:(NillumSVD+Startillum-1),:),transpose(X)) %test on H
    [EEG]=errorICA(G,Gest) %test on G
    [EEU]=errorICA(G,U) %test the SVD separation on G
    [EEV]=errorICA(conj(H(Startillum:(NillumSVD+Startillum-1),:)),V) %test the SVD separation on H
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Save the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SepData{1}=U; %store U to make an image with the SVD
SepData{2}=Gest; %store Gest to make an image with the ICA
SepData{3}=NscatSVD; %store the new number of singular values after filtering

save SepData SepData


