function res=MultPlain(arg,Tflag,data);
if Tflag==0,
    res=data.A*arg;
else
    res=(arg'*data.A)';
end;