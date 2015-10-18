function res=FFTMult(x,Tflag,data);
myi=sqrt(-1);
k=data.ydim/2;
n=data.xdim;
if Tflag,
    tmp=zeros(n,1);
    tmp(data.freq)=x(1:k)+myi*x(k+1:end);
    res=data.FFTscale*real(n*ifft(tmp));
else
    tmp=fft(x);
    res=data.FFTscale*[real(tmp(data.freq));imag(tmp(data.freq))];
end;