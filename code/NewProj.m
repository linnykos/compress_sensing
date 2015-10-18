function s=NewProj(alpha,n,xi);
%%% Corrected as of March 3, 2012
SMALL=1.e-8;
s=xi.^alpha;
if sum(s)<=1,
    return;
end;
[s,ind]=sort(xi);
g=sum((s(2:n)-s(1)).^alpha);
if g<=1,
    lmin=s(1)-(1/n)^(1/alpha);
    lmax=s(1);
    kmin=0;
else
    kmin=1;
    kmax=n;
    while(kmax-kmin>1);
        k=floor(0.5*(kmin+kmax));
        g=sum((s(k+1:n)-s(k)).^alpha);
        if g>=1,
            kmin=k;
        else
            kmax=k;
        end;
    end;
    lmin=s(kmin);
    lmax=s(kmin+1);
end;
lm=lmin;
nnwt=0;
bt=alpha-1;
while(1==1)
    tmp=(s(kmin+1:end)-lm).^bt;
    f=tmp'*(s(kmin+1:end)-lm);
    df=-alpha*sum(tmp);
    %disp(sprintf('%.3e',f-1));
    if abs(f-1)<SMALL*1.e-3*f,
        tmp=tmp.*(s(kmin+1:end)-lm);
        s(ind(1:kmin))=0;
        s(ind(kmin+1:end))=tmp;
        break;
    end;
    nnwt=nnwt+1;
    lm=lm-(f-1)/df;
    lm=max(lm,lmin);
    lm=min(lm,lmax);
end;
    

            
