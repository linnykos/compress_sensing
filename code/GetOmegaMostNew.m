function [omega,domega]=GetOmegaMostNew(z,dim,scale,p);
domega=zeros(dim.z,1);
if scale.type(1)=='e',
    domega(1:dim.x)=scale.x*log(z(1:dim.x));
    omega=z(1:dim.x)'*domega(1:dim.x);
    domega(1:dim.x)=domega(1:dim.x)+scale.x;
else
    domega(1:dim.x)=scale.x*sign(z(1:dim.x)).*abs(z(1:dim.x)).^(scale.dgx-1);
    omega=domega(1:dim.x)'*z(1:dim.x)/scale.dgx;
end;
if p==2,
    domega(dim.x+1:end)=scale.y*z(dim.x+1:end);
    omega=omega+0.5*z(dim.x+1:end)'*domega(dim.x+1:end);
elseif scale.type(1)=='e',
    domega(dim.x+1:end)=scale.y*log(z(dim.x+1:end));
    omega=omega+z(dim.x+1:end)'*domega(dim.x+1:end);
    domega(dim.x+1:end)=domega(dim.x+1:end)+scale.y;
else
    domega(dim.x+1:end)=scale.y*sign(z(dim.x+1:end)).*abs(z(dim.x+1:end)).^(scale.dgy-1);
    omega=omega+domega(dim.x+1:end)'*z(dim.x+1:end)/scale.dgy;
end;
