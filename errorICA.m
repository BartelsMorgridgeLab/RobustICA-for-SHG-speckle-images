function [E Er C] = Error(H,X)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Ni Ns]=size(H);

for i=1:Ns
    H(:,i)=H(:,i)/sqrt(sum(abs(H(:,i)).^2));
    X(:,i)=X(:,i)/sqrt(sum(abs(X(:,i)).^2));
end

C=X'*H;

for i=1:Ns
    for j=1:Ns

        Error(i,j)=sum(abs(H(:,i)-C(j,i)*X(:,j)).^2);

    end
end

dim=Ns;
E=0;
%[M,I1]=min(Error)
%[MM,I2]=min(M)
%Error

Er=zeros(Ns,1);

while dim>0
    [M,I1]=min(Error);
    [MM,I2]=min(M);
    E=E+Error(I1(I2),I2);
    Er(dim)=Error(I1(I2),I2);
    Error(I1(I2),:)=[];
    Error(:,I2)=[];
    %Error
    dim=dim-1;    
end

E=E/Ns;
end