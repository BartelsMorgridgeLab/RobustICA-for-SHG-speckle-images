function x=deconvtvGelma(Psi,yy,tau,toler)

y=yy(:);
N=length(y);
[Nc Nc]=size(Psi);

%dt=1/sqrt(sum(sum(abs(Psi.^2))))*Nc*Nc;


%Initialization
x=zeros(N,1);
%Ax=convo(Psi,reshape(x,Nc,Nc),Nc^4);
%z=Ax;
ww=x;
ak=1;

%tau=mean(convostar(Psi,yy,Nc^4))


LastRelativeError=0.1;


%convo = @(K,xx) ifft2(fft2(K).*fft2(xx));
%convostar = @(K,xx) ifft2(conj(fft2(K)).*fft2(xx));

while (LastRelativeError>toler)
    
    x_old=x;
    Ax=convo(Psi,reshape(x,Nc,Nc),Nc^4);
    Astar=convostar(Psi,reshape(Ax-y,Nc,Nc),Nc^4);
    ww=x-2*ak*Astar;
    x=eta(ww,tau*ak);
    
    ak_old=ak;
    ak=(1+sqrt(1+4*ak^2))/2;
    xi=x+(ak_old-1)/ak*(x-x_old);

    % figure(1);
    % mesh(abs(Psi))
    % figure(2);
    % plot(abs(Ax))
    % figure(3);
    % plot(abs(x))
    % pause
    %znew=(y-Ax)*dt+z;
    %Astar=convostar(Psi,reshape(y-Ax+z,Nc,Nc),Nc^4);
    %norm(Astar)
    %ww=Astar*dt+x;
    %z=znew;
    %xold=x;
   %  norm(ww)
    
    %norm(x)
    %pause

    LastRelativeError=mean(abs(x-x_old))/mean(abs(x))
end

x=fftshift(reshape(x,Nc,Nc));
end

function c=convo(K,xx,NN)
    
    c=ifft2(fft2(K).*fft2(xx));
    c=c(:)/NN;
end

function c=convostar(K,xx,NN)

    c=ifft2(conj(fft2(K)).*fft2(xx));
    c=c(:)/NN;
end