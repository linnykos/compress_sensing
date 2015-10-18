function [ubn,lbn,bestn,rhoubn,Best]=UpdateBoundsGBSPRelease(dim,ub,lb,xyf,best,data,rho,rhoub,dgf,delta);
bestn=best;
rhoubn=rhoub;
cfrho=-xyf.yy'*data.b-delta;
cf=-max(abs(xyf.ATyy));
lbn=cf+cfrho*rho;
if cfrho>0,
    rhoubn=min(rhoubn,-cf/cfrho);
end;    
if lbn>lb,
    bestn.y=xyf.yy;
    bestn.ATy=xyf.ATyy;
else
    lbn=lb;
end;
if data.fit==inf,
    ubn=max(abs(xyf.Axx-rho*data.b))-rho*delta;
else
    ubn=norm(xyf.Axx-rho*data.b)-rho*delta;
end;
if ubn<ub,
    bestn.x=xyf.xx;
    bestn.Ax=xyf.Axx;
else
    ubn=ub;
end;
%%%
if xyf.w>1,
    cfrho=-xyf.y'*data.b/xyf.w-delta;
    cf=-max(abs(xyf.ATy))/xyf.w;
    lbnn=cf+cfrho*rho;
    if cfrho>0,
        rhoubn=min(rhoubn,-cf/cfrho);
    end;
    if lbnn>lbn,
        bestn.y=xyf.y/xyf.w;
        bestn.ATy=xyf.ATy/xyf.w;
    else
        lbnn=lbn;
    end;
    if data.fit==inf,
        ubnn=max(abs(xyf.Ax/xyf.w-rho*data.b))-rho*delta;
    else
        ubnn=norm(xyf.Ax/xyf.w-rho*data.b)-rho*delta;
    end;
    if ubnn<ubn,
        bestn.x=xyf.x/xyf.w;
        bestn.Ax=xyf.Ax/xyf.w;
    else
        ubnn=ubn;
    end;
    ubn=ubnn;
    lbn=lbnn;
end;
if dgf(1)=='e',
    Best=zeros(dim.z,1);
    Best(1:dim.x)=[max(bestn.x,0);max(-bestn.x,0)];
    tmp=max((1-sum(Best(1:dim.x)))/dim.x,0);
    Best(1:dim.x)=Best(1:dim.x)+tmp;
    if data.fit==inf,
        Best(dim.x+1:dim.x+dim.y)=[max(bestn.y,0);max(-bestn.y,0)];
        tmp=max((1-sum(Best(dim.x+1:dim.x+dim.y)))/dim.y,0);
        Best(dim.x+1:dim.x+dim.y)=Best(dim.x+1:dim.x+dim.y)+tmp;
    else
        Best(dim.x+1:dim.x+dim.y)=bestn.y;
    end;
else
    Best=[bestn.x;bestn.y];
end;


