function x=deconvtvGelma(Psi,yy,tau,toler)

y=yy(:);
N=length(y);
[Nc Nc]=size(Psi);

dt=1/sqrt(sum(sum(abs(Psi.^2))))*Nc*Nc;


%Initialization
x=zeros(N,1);
Ax=convo(Psi,reshape(x,Nc,Nc),Nc^4);
z=Ax;
ww=x;
%tau=mean(convostar(Psi,yy,Nc^4));


LastRelativeError=0.1;


%convo = @(K,xx) ifft2(fft2(K).*fft2(xx));
%convostar = @(K,xx) ifft2(conj(fft2(K)).*fft2(xx));

while (LastRelativeError>toler)

    Ax=convo(Psi,reshape(x,Nc,Nc),Nc^4);
    % figure(1);
    % mesh(abs(Psi))
    % figure(2);
    % plot(abs(Ax))
    % figure(3);
    % plot(abs(x))
    % pause
    znew=(y-Ax)*dt+z;
    Astar=convostar(Psi,reshape(y-Ax+z,Nc,Nc),Nc^4);
    %norm(Astar)
    ww=Astar*dt+x;
    z=znew;
    xold=x;
   %  norm(ww)
    x=eta(ww,tau*dt);
    %norm(x)
    %pause

    LastRelativeError=mean(abs(x-xold))/mean(abs(x))
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