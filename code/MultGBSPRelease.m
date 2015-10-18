function res=MultGBSPRelease(arg,Tflag,data,dgf);
if Tflag==0,
    if (data.fit==inf)&(dgf(1)=='e'),
        res=zeros(2*data.ydim,1);
        %str=sprintf('res(1:data.ydim)=%s((arg(1:data.xdim)-arg(data.xdim+1:end)),Tflag,data);',...
            %data.Mult);
        %eval(str);
        res(1:data.ydim)=feval(data.Mult,arg(1:data.xdim)-arg(data.xdim+1:end),...
            Tflag,data);
        %res(1:data.ydim)=data.A*(arg(1:data.xdim)-arg(data.xdim+1:end));
        res(data.ydim+1:end)=-res(1:data.ydim);
    elseif dgf(1)=='e',
        %str=sprintf('res=%s(arg(1:data.xdim)-arg(data.xdim+1:end),Tflag,data);',...
            %data.Mult);
        %eval(str);
        res=feval(data.Mult,arg(1:data.xdim)-arg(data.xdim+1:end),Tflag,data);
    else
        %str=sprintf('res=%s(arg,Tflag,data);',...
            %data.Mult);
        %eval(str);
        res=feval(data.Mult,arg,Tflag,data);
    end;
else
    if dgf(1)=='e',
        res=zeros(2*data.xdim,1);
        if data.fit==inf,
            %str=sprintf('res(1:data.xdim)=%s(arg(1:data.ydim)-arg(data.ydim+1:end),Tflag,data);',...
                %data.Mult);
            %eval(str);
            res(1:data.xdim)=feval(data.Mult,arg(1:data.ydim)-arg(data.ydim+1:end),...
                Tflag,data);
        else
            %str=sprintf('res(1:data.xdim)=%s(arg,Tflag,data);',...
                %data.Mult);
            %eval(str);
            res(1:data.xdim)=feval(data.Mult,arg,Tflag,data);
        end;
        res(data.xdim+1:end)=-res(1:data.xdim);
    else
        %str=sprintf('res=%s(arg,Tflag,data);',...
            %data.Mult);
        %eval(str);
        res=feval(data.Mult,arg,Tflag,data);
    end;
end;
res=full(res);