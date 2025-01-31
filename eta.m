function z=eta(x,tau)
y=round((sign(x.^2-tau^2)+1)/2).*sign(x);
zz=abs(y);
z=zz.*(x-tau*y);