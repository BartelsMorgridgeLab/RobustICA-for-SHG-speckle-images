
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Some stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Susceptility of the scatterers (here just equal to one)
rho=(1+0.*(rand(Nscat,1)-0.5)).*exp(0.*2*pi*1i*rand(Nscat,1));
rho=diag(rho);

%%%Initialization of the forward and backward fields
G=zeros(Nc^2,Nscat);
H=zeros(Nillum,Nscat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Backward fields G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%loading illuminations
%%%% we load smaller files when fewer illuminations are needed

if Nillum<=1000
    load IllumData1000;
    IllumData=IllumData1000;
    clear IllumData1000

elseif (Nillum<=2500)&&(Nillum>1000)
    load IllumData2500;
    IllumData=IllumData2500;
    clear IllumData2500

else
    load IllumData
end

%%%% Initialization
E=cell(Nscat,1);
EB=cell(Nscat,1);

for i=1:Nscat
    E{i}=zeros(N); 
end

Eref=zeros(N); 

%%%% Definition of the scatterers (choose which type in parameters file)
%%%% E is equal to one where there is a scatterer on the sample plane
%%%% Pos stores the exact XY positions of the scatterers

if ScatType=='2scat'

    %%% 2 scatterers
    E{1}(N/2,N/2)=1; E{2}(N/2,N/2+shift)=1;
    Ry(1)=0;Rx(1)=0;
    Rx(2)=shift; Ry(2)=0;
    Pos= [(Rx)' (Ry)']*Delta_x;

elseif ScatType=='RandScat'
    
    %%%Random scatterers
    Rx = randi([-windowScat windowScat],Nscat,1);
    Ry = randi([-windowScat windowScat],Nscat,1);
    Pos= [(Rx) (Ry)]*Delta_x;

    for i=1:Nscat
     E{i}(N/2+Ry(i),N/2+Rx(i))=1;
    end

elseif ScatType=='cross'

    %cross
    for i=1:NL
        E{i}(N/2,N/2+shift*i)=1;
        E{i+NL}(N/2+shift*i,N/2)=1;
        Rx(i)=shift*i; Ry(i)=0;
        Rx(i+NL)=0; Ry(i+NL)=shift*i;
    end

    for i=1:NL
     E{i+2*NL}(N/2,N/2-shift*i)=1;
     E{i+3*NL}(N/2-shift*i,N/2)=1;
     Rx(i+2*NL)=-shift*i; Ry(i+2*NL)=0;
     Rx(i+3*NL)=0; Ry(i+3*NL)=-shift*i;
    end

    E{4*NL+1}(N/2,N/2)=1;
    Rx(4*NL+1)=0;Ry(4*NL+1)=0;

    Pos= [(Rx)' (Ry)']*Delta_x;

elseif ScatType=='circle'

    %%% circle
    for i=1:Nscat    
      E{i}(N/2+ceil(shift*cos(2*pi*i/Nscat)),N/2+ceil(shift*sin(2*pi*i/Nscat)))=1;
      Ry(i)=shift*cos(2*pi*i/Nscat);Rx(i)=shift*sin(2*pi*i/Nscat);
    end

    Pos= [(Rx)' (Ry)']*Delta_x;

elseif ScatType=='2groups'

    %scatterers in 2 groups
    E{1}(N/2,N/2+3*shift)=1; E{2}(N/2,N/2-shift)=1; E{3}(N/2,N/2)=1;
    E{4}(N/2-shift,N/2+shift)=1; E{5}(N/2,N/2+shift)=1; E{6}(N/2+shift,N/2+shift)=1;
    E{7}(N/2-shift,N/2-shift)=1;E{8}(N/2+shift,N/2-shift)=1;E{9}(N/2-shift,N/2)=1;
    E{10}(N/2+shift,N/2)=1;
 
    Rx = [3*shift ; -shift ; -shift; -shift ; 0 ; 0; 0 ; shift ; shift ; shift];
    Ry = [0 ; -shift ; 0 ; shift ; -shift ; 0 ; shift ;-shift ; 0 ; shift];
    Pos= [(Rx) (Ry)]*Delta_x;

elseif ScatType=='smiley'

    %% Smiley
    load Smiley
    %Smiley=Smiley2;
    [u v]=find(Smiley~=0);
    for i=1:Nscat
     E{i}(u(i),v(i))=1;
     Rx(i)=u(i)-N/2; Ry(i)=v(i)-N/2;
    end
    Pos= [(Rx) (Ry)]*Delta_x;

end

%%%%Reference
for i=1:Nscat
    Eref=Eref+E{i};
end

save Pos Pos
save Eref Eref

Phase1=IllumData{2};
Phase2=IllumData{3};

%%%% Computing the backward fields for each scatterers
parfor i=1:Nscat
    ip=i;
    fprintf('Backward Field calculation scatterer %d out of %d \n',ip,Nscat)
    temp=BackwardField(Phase1.^2,Phase2.^2,padding_b,2*k0,L0,L1,R,Xx,Xy,Ux_pad,Uy_pad,E{i});
    EB{i}=temp((N/2-width-1):(N/2+width),(N/2-width-1):(N/2+width)); %crop the field     
    G(:,i)=(EB{i}(:)); %Flattening each field
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Illuminations values at the scatterers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

center=widthSample; %needed to shift the origin to the center

for j=1:Nillum

    EF=IllumData{1}{j};

    %%Values of the forward field at the location of the scatterers

    if ScatType=='2scat'

        %%%2 scatterers
        H(j,1)=(EF(center,center))^2;H(j,2)=(EF(center,center+shift))^2;
    
    elseif ScatType=='RandScat'
    
        %%Random scatterers
        for i=1:Nscat
          H(j,i)=(EF(center+Ry(i),center+Rx(i)))^2;
        end

    elseif ScatType=='cross'

        %%%cross
        for i=1:NL
          H(j,i)=(EF(center,center+shift*i))^2;
          H(j,i+NL)=(EF(center+shift*i,center))^2;
          H(j,i+2*NL)=(EF(center,center-shift*i))^2;
          H(j,i+3*NL)=(EF(center-shift*i,center))^2;
        end

        H(j,4*NL+1)=(EF(center,center))^2;

    elseif ScatType=='circle'

        %%circle
        for i=1:Nscat
          H(j,i)=(EF(center+ceil(shift*cos(2*pi*i/Nscat)),center+ceil(shift*sin(2*pi*i/Nscat))))^2;
        end
    

    elseif ScatType=='2groups'

        %scatterers in 2 groups
        H(j,1)=(EF(center,center+3*shift))^2;H(j,2)=(EF(center,center-shift))^2;H(j,3)=(EF(center,center))^2;     
        H(j,4)=(EF(center-shift,center+shift))^2;H(j,5)=(EF(center,center+shift))^2;H(j,6)=(EF(center+shift,center+shift))^2;     
        H(j,7)=(EF(center-shift,center-shift))^2;H(j,8)=(EF(center+shift,center-shift))^2;H(j,9)=(EF(center-shift,center))^2;     
        H(j,10)=(EF(center+shift,center))^2;

    elseif ScatType=='smiley'

        %%smiley
        for i=1:Nscat
            H(j,i)=(EF(u(i)-N/2+widthSample,v(i)-N/2+widthSample))^2;
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Reflection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ReflData=cell(5,1);

ReflData{1}=G;
ReflData{2}=rho;
ReflData{3}=H;
ReflData{4}=Eref;

save ReflData ReflData