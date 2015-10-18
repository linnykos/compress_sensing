function [w,omn,domn,lambda,nBsteps]=ProxGBSPMostNew(xi,linf,dim,scale,p);
%%%%
small=1.e-5;
nBsteps=0;
%%%%
[w,omn,domn]=BaseProxGBSPMostNew(xi,dim,scale,p);
lambda=0;
if linf.c==-inf,
    return;
end;
dv=linf.cf'*w+linf.c;
if dv<=0,
    return;
end;
%%%
lmn=0;
lmx=1;
flag=1;
best=-inf;
while(lmx-lmn>small*max(lmx,1))
    nBsteps=nBsteps+1;
    if flag==1,
        lm=lmx;
    else
        lm=0.5*(lmn+lmx);
    end;
    [ww,omnn,domnn]=BaseProxGBSPMostNew(xi+lm*linf.cf,dim,scale,p);
    dv=linf.cf'*ww+linf.c;
    v=omnn+xi'*ww+lm*dv;
    if v>best,
        w=ww;
        omn=omnn;
        domn=domnn;
        lambda=lm;
        best=v;
    end;
    if dv>=0,
        if flag==1
            lmn=lmx;
            lmx=2*lmx;
        else
            lmn=lm;
        end;
    else
        flag=0;
        lmx=lm+(best-v)/dv;
    end;
end;
            
    
