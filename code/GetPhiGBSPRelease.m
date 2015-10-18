function [Phi,xyf]=GetPhiGBSPRelease(z,dim,data,rho,dgf);
Phi=[MultGBSPRelease(z(dim.x+1:end),1,data,dgf);...
    -MultGBSPRelease(z(1:dim.x),0,data,dgf)];
if dgf(1)=='e',
    xyf.x=z(1:dim.xf)-z(dim.xf+1:dim.x);
else
    xyf.x=z(1:dim.x);
end;
xyf.xx=xyf.x;
xyf.w=1;
xyf.Ax=-Phi(dim.x+1:dim.x+dim.yf);
xyf.Axx=xyf.Ax;
Phi(dim.x+1:dim.x+dim.yf)=Phi(dim.x+1:dim.x+dim.yf)+rho*data.b;
if (data.fit==inf)&(dgf(1)=='e'),
    Phi(dim.x+dim.yf+1:dim.x+dim.y)=-Phi(dim.x+1:dim.x+dim.yf);
end;
%%%
if (data.fit==inf)&(dgf(1)=='e'),
    xyf.y=z(dim.x+1:dim.x+dim.yf)-z(dim.x+dim.yf+1:dim.x+dim.y);
else
    xyf.y=z(dim.x+1:dim.x+dim.y);
end;
xyf.yy=xyf.y;
xyf.ATy=Phi(1:dim.xf);
xyf.ATyy=xyf.ATy;
if dgf(1)=='e',
    Phi(dim.xf+1:dim.x)=-Phi(1:dim.xf);
end;
