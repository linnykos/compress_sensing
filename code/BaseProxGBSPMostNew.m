function [z,omega,domega]=BaseProxGBSPMostNew(xi,dim,scale,p);
domega=zeros(dim.z,1);
z=zeros(dim.z,1);
if scale.type(1)=='e',
    domega(1:dim.x)=-xi(1:dim.x)/scale.x;
    mx=max(domega(1:dim.x));
    z(1:dim.x)=exp(domega(1:dim.x)-mx);
    sx=sum(z(1:dim.x));
    z(1:dim.x)=z(1:dim.x)/sx;
    lnsx=log(sx);
    domega(1:dim.x)=scale.x*(domega(1:dim.x)-mx-lnsx);
    omega=z(1:dim.x)'*domega(1:dim.x);
    domega(1:dim.x)=domega(1:dim.x)+scale.x;
else
    domega(1:dim.x)=-xi(1:dim.x)/scale.x;
    z(1:dim.x)=sign(domega(1:dim.x)).*NewProj(1/(scale.dgx-1),dim.x,abs(domega(1:dim.x)));
    domega(1:dim.x)=scale.x*sign(domega(1:dim.x)).*abs(z(1:dim.x)).^(scale.dgx-1);
    omega=domega(1:dim.x)'*z(1:dim.x)/scale.dgx;
end;
%%%
if p==2,
    z(dim.x+1:end)=-xi(dim.x+1:end)/scale.y;
    nz=norm(z(dim.x+1:end));
    if nz>1,
        z(dim.x+1:end)=z(dim.x+1:end)/nz;
    end;
    domega(dim.x+1:end)=scale.y*z(dim.x+1:end);
    omega=omega+0.5*z(dim.x+1:end)'*domega(dim.x+1:end);
elseif scale.type(1)=='e',
    domega(dim.x+1:end)=-xi(dim.x+1:end)/scale.y;
    mx=max(domega(dim.x+1:end));
    z(dim.x+1:end)=exp(domega(dim.x+1:end)-mx);
    sx=sum(z(dim.x+1:end));
    z(dim.x+1:end)=z(dim.x+1:end)/sx;
    lnsx=log(sx);
    domega(dim.x+1:end)=scale.y*(domega(dim.x+1:end)-mx-lnsx);
    omega=omega+z(dim.x+1:end)'*domega(dim.x+1:end);
    domega(dim.x+1:end)=domega(dim.x+1:end)+scale.y;
else
    domega(dim.x+1:end)=-xi(dim.x+1:end)/scale.y;
    z(dim.x+1:end)=sign(domega(dim.x+1:end)).*...
        NewProj(1/(scale.dgy-1),dim.y,abs(domega(dim.x+1:end)));
    domega(dim.x+1:end)=scale.y*sign(domega(dim.x+1:end)).*abs(z(dim.x+1:end)).^(scale.dgy-1);
    omega=omega+domega(dim.x+1:end)'*z(dim.x+1:end)/scale.dgy;
end;

