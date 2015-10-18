function sol=MirrorProx(data,control)
%% version of 01/09/2012
%MIRRORPROX implements the deterministic Mirror-Prox algorithm for 
% l_1 minimization: solving for with p=2 or p=inf the problem 
%       Opt = min{||x||_1: ||A*x - b||_p< delta},                      (1)
% or
%       Opt(rho) = min{||A*x - b||_p: ||x||_1 < 1/rho}                 (2)
%
% On successful run, MirrorProx returns a vector x such that
%       ||x||_1< (1 + omega)*Opt and ||A*x - b||_p< delta + eps,       (3)
%                                                           
% USAGE:
%   res = MirrorProx(data) 
%   res = MirrorProx(data, control)
%
% INPUTS: 
% data: data structures with fields 
%   fit:    (double, either 2 or inf) [2] type of the fit
%   xdim:   (positive int) problem X-dimension 
%   ydim:   (positive int) problem Y-dimension
%   b:      (DATA.ydim-double) problem Y
%   L:      (double) 1,p-norm off A: L=max {||Column(A)_j||_p, 1<j<n}
%   delta:  (double) [0.0005], delta=delta in the problem definition.
%           Note: data.delta is irrelevant when control.rhofixed=’on’
%   Mult:   (function_handle) handle to the 'multiplication oracle' -- 
%           a user-supplied function computing matrix-vector products 
%           involving A and A' (see MANUAL.PDF for details)
%
% control: algorithm control stucture with fields 
%   rhofixed:   (’on’/’off’) [’off’] 
%               ’on’: solving problem (2)
%               ’off’: solving problem (1)
%   rho:        (>0) [1] when control.rhofixed='on': rho in (3)
%               when control.rhofixed='off': control.rho is irrelevant
%   accuracy:   (’abs’/’rel’) [’rel’] 
%               accuracy=’abs’ means that a solution of absolute accuracy 
%               control.eps is sought;
%               accuracy=’rel’ means that a solution of relative accuracy 
%               control.eps is sought
%   eps:        (> 0) [0.0005] 
%               required absolute (when control.accuracy=’abs’)
%               or relative (when control.accuracy=’rel’) accuracy
%   omega:      (> 0) [0]
%               when control.rhofixed='off': control.omega is omega in (3),
%               when control.rhofixed='on': control.omega is irrelevant
%               (but should be well defined and positive)
%   ItrMax:     (1,2,...) [50000] maximal # of steps
%   cpumax:     (0,inf) [1800] maximum running time, sec
%   kappa:      (0.01,0.99) [0.75] (see MANUAL.PDF Section 3 for details)
%   printlevel: (0/1) [0] printlevel (0 – no display output)
%   dgf:        (’e’/’p’) [’p’] switch from ’p’ to ’e’ is expected to increase
%               somehow the iteration count while reducing the
%               computational effort per iteration 
%               For more details see MANUAL.PDF
%
% OUTPUT: 
% res: structure containing the results of the run
%   ok:         =0 for successful run (desired accuracy is reached),
%               =1 when maximum iteration count/running time is reached
%   rhofixed:   ’on’ when control.rhofixed is on, and ’off’ otherwise
%   x:          array with resulting primal solution x
%   A*x:         array A*x
%   rho:        when control.rhofixed='off': the final value of rho, 
%               (see section 3 of MANUAL.PDF)
%               when control.rhofixed='on': control.rho
%   ell1norm:   l_1 norm of the solution
%   y:          array with the resulting dual approximate solution y, 
%               (see section 3 of MANUAL.PDF)
%   ATy:        array A'*y
%   err.R:      when control.rhofixed='off': upper bound on [||x||_1- Opt]
%               when control.rhofixed='on: 0
%   err.abs:    upper bound on |A*x - b|_p- delta
%   err.rel:    when control.rhofixed='off': upper bound on 
%               [||A*x - b||_p- delta]/||A||_{1,p}/Opt
%               when control.rhofixed='on': inf
%   opt_lwb:    when control.rhofixed='off': a lower bound on Opt, see (1)
%               when control.rhofixed='on': 0
%   res:        ||A*x - b||_p
%   res_rel:    ||A*x - b||_p/||b||_p
%   res_lwb:    a lower bound on Opt(rho) for rho = res.rho, see (2)
%   cpu:        running time, sec
%   Steps:      total number of iterations
%   Calls:      total number of multiplications A'*y and A*x
%   Stages:     total number of Stages, (see section 3 of MANUAL.PDF)
%
%
nin=nargin;
error(nargchk(1, 2, nin, 'struct'))
if nin<2, control=[]; end;
%
sol.ok=1;
sol.p=[];
sol.x=[];
sol.Ax=[];
sol.y=[];
sol.ATy=[];
sol.ell1norm=[];
sol.rho=[];
sol.opt_lwb=[];
sol.res_lwb=[];
sol.res=[];
sol.res_rel=[];
sol.err.R=[];
sol.err.abs=[];
sol.err.rel=[];
sol.Steps=[];
sol.Calls=[];
sol.Stages=[];
sol.cpu=[];
sol.rhofixed=[];
sol.control=[];
%
if isfield(control,'omega'),
    cntr.omega=control.omega;
    if cntr.omega<0,
        cntr.omega=0;
    end;
else
    cntr.omega=0;
end;
if isfield(control,'kappa'),
    cntr.kappa=control.kappa;
    if cntr.kappa<0.05,
        cntr.kappa=0.05;
    end;
    if cntr.kappa>0.95,
        cntr.kappa=0.95;
    end;
else
    cntr.kappa=0.75;
end;
if isfield(control,'ItrMax'),
    cntr.StepMax=control.ItrMax;
    if (cntr.StepMax~=floor(cntr.StepMax))||(cntr.StepMax<1)
        cntr.StepMax=10000;
    end;
else
    cntr.StepMax=50000;
end;
%%% ai
if isfield(data,'delta'),
    cntr.delta=data.delta;
    if cntr.delta<0,
        cntr.delta=0;
    end;
else
    cntr.delta=0.0005;
end;
if ~isfield(data,'fit'),
    data.fit=2;
end
if isfield(control,'accuracy'),
    cntr.accuracy=control.accuracy;
    if cntr.accuracy(1)~='a',
        cntr.accuracy='rel';
    end;
else
    cntr.accuracy='rel';
end;
if isfield(control,'eps'),
    cntr.eps=control.eps;
    if cntr.eps<0,
        cntr.eps=0;
    end;
else
    cntr.eps=0.0005;
end;
cntr.gamma_up=1.2;
cntr.gamma_down=0.8;
if isfield(control,'dgf'),
    cntr.dgf=control.dgf;
    if cntr.dgf(1)~='e',
        cntr.dgf='p';
    end;
else
    cntr.dgf='p';
end;
if isfield(control,'rhofixed'),
    cntr.rhofixed=control.rhofixed;
    if cntr.rhofixed(2)=='n'
        cntr.rhofixed='on';
        if isfield(control,'rho'),
            cntr.rho=control.rho;
        else
            cntr.rho=1;
        end;
    else
        cntr.rhofixed='off';
        cntr.rho=1;
    end;
else
    cntr.rhofixed='off';
    cntr.rho=1;
end;
if isfield(control,'cpumax'),
    cntr.cpumax=control.cpumax;
    if cntr.cpumax<100,
        cntr.cpumax=100;
    end;
else
    cntr.cpumax=1800;
end;
%%%
cntr.mode='a';
cntr.printlevel=100;
cntr.nAverages=12;
cntr.restart=0.5000;
cntr.restarttype='b';
cntr.postolerance=0.2000;
if isfield(control,'printlevel'),
    cntr.printlevel=100*control.printlevel;
    if (cntr.printlevel~=floor(cntr.printlevel))||(cntr.printlevel<0)
        cntr.printlevel=0;
    end;
else
    cntr.printlevel=0;
end;
if cntr.dgf(1)=='e',
    cntr.restart=0.5;
    cntr.restarttype='b';
else
    if data.fit==2,
        cntr.restart=1;
        cntr.restarttype='b';5; %%% ai ???
    else
        cntr.restart=0.5;
        cntr.restarttype='b';
    end;
end;
%%%
Starts=[2.^(0:1:cntr.nAverages),ones(1,cntr.nAverages)];
Scores=zeros(2*cntr.nAverages,1);
%%%
ACCMAX=100000;
%%%
dim.xf=data.xdim;
dim.yf=data.ydim;
if cntr.dgf(1)=='p'
    dim.x=dim.xf;
    dim.y=dim.yf;
elseif cntr.dgf(1)=='e'
    dim.x=2*dim.xf;
    if data.fit==2,
        dim.y=dim.yf;
    else
        dim.y=2*dim.yf;
    end;
end;
dim.z=dim.x+dim.y;
%%%
if data.fit==inf,
    errorbase=max(abs(data.b))-cntr.delta;
else
    errorbase=norm(data.b)-cntr.delta;
end;
%%%
ireport=0;
decrement=sqrt(sqrt(0.5));
errortarget=errorbase*decrement;
%%% scaling
if cntr.dgf(1)=='e',
    scale.type='e';
    Omega.x=log(dim.x);
    if data.fit==inf,
        Omega.y=log(dim.y);
    else
        Omega.y=1/2;
    end;
else
    scale.type='p';
    scale.dgx=1+1/log(dim.x);
    scale.cfx=(dim.x^(scale.dgx-1))/(scale.dgx-1);
    Omega.x=scale.cfx;
    if data.fit==inf,
        scale.dgy=1+1/log(dim.y);
        scale.cfy=(dim.y^(scale.dgy-1))/(scale.dgy-1);
        Omega.y=scale.cfy;
    else
        Omega.y=1/2;
    end;
end;
L.xy=data.L;
cL=2*L.xy*sqrt(Omega.x*Omega.y);
scale.x=L.xy*sqrt(Omega.x*Omega.y)/(Omega.x*cL);
scale.y=L.xy*sqrt(Omega.x*Omega.y)/(Omega.y*cL);
if scale.type(1)=='p',
    scale.x=scale.x*scale.cfx;
    if data.fit==inf,
        scale.y=scale.y*scale.cfy;
    end;
end;
%%% initial point
zn=zeros(dim.z,1);
if scale.type(1)=='e',
    zn(1:dim.x)=1/dim.x;
    if data.fit==2,
        zn(dim.x+1:dim.x+dim.y)=0;
    else
        zn(dim.x+1:dim.x+dim.y)=1/dim.y;
    end;
else
    zn=zeros(dim.z,1);
end;
%%%
Best=zeros(dim.z,1);
%%%
xyfav=cell(cntr.nAverages,1);
for i=1:2*cntr.nAverages,
    xyfav{i}.x=zeros(dim.xf,1);
    xyfav{i}.Ax=zeros(dim.yf,1);
    xyfav{i}.y=zeros(dim.yf,1);
    xyfav{i}.ATy=zeros(dim.xf,1);
    xyfav{i}.xx=zeros(dim.xf,1);
    xyfav{i}.Axx=zeros(dim.yf,1);
    xyfav{i}.yy=zeros(dim.yf,1);
    xyfav{i}.ATyy=zeros(dim.xf,1);
    xyfav{i}.w=0;
end;
%%%
wf=zeros(dim.xf+dim.yf,1);
wad=zeros(dim.xf+dim.yf,1);
Awad=zeros(dim.xf+dim.yf,1);
Awf=zeros(dim.xf+dim.yf,1);
%%%
best.x=zeros(dim.xf,1);
best.y=zeros(dim.yf,1);
best.Ax=zeros(dim.yf,1);
best.ATy=zeros(dim.xf,1);
%%%
gamma_safe=1/cL;
gamma=gamma_safe;
sumgammafact=0;
sumgammabas=0;
%%%
lb=-inf;
ub=inf;
%%%
if (cntr.dgf(1)=='e')&&(cntr.restart==1),
    restartf=1-1.e-6;
else
    restartf=cntr.restart;
end;
tstart=cputime;
%%%
if cntr.rhofixed(2)=='f',
    flagrho=1;
    if data.fit==inf,
        if max(abs(data.b))<=cntr.delta,
            error('Trivial problem with zero solution');
        end;
        rho=data.L/(max(abs(data.b))-cntr.delta);
    else
        if norm(data.b)<=cntr.delta,
            error('Trivial problem with zero solution');
        end;
        rho=data.L/(norm(data.b)-cntr.delta);
    end;
else
    flagrho=0;
    rho=cntr.rho;
end;
if cntr.accuracy(1)=='a',
    tolerance=cntr.eps*rho;
else
    tolerance=cntr.eps*data.L/(1+cntr.omega);
end;

%%%
Stage=0;
stop=0;
Steps=0;
Calls=0;
rhofactor=1;
rhoub=rho;
nBsteps=0;
nProx=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (1==1)
    Stage=Stage+1;
    CallsLoc=0;
    tme=cputime-tstart;
    if cntr.printlevel,
        srho=MySprintf(rho);
        fprintf(1,'***** %6.2f: Stage=%d rho=%s Steps=%d Calls=%d\n',...
            tme,Stage,srho,Steps,Calls); %%% ai
%         disp(sprintf('***** %6.2f: Stage=%d rho=%s Steps=%d Calls=%d',...
%             tme,Stage,srho,Steps,Calls));
   end;
    %% updating bundle
    Best(1:dim.x)=Best(1:dim.x)/rhofactor;
    %%
    for i=1:cntr.nAverages,
        xyfav{i}.x=0*xyfav{i}.x;
        xyfav{i}.xx=0*xyfav{i}.xx;
        xyfav{i}.y=0*xyfav{i}.y;
        xyfav{i}.yy=0*xyfav{i}.yy;
        xyfav{i}.Ax=0*xyfav{i}.Ax;
        xyfav{i}.Axx=0*xyfav{i}.Axx;
        xyfav{i}.ATy=0*xyfav{i}.ATy;
        xyfav{i}.ATyy=0*xyfav{i}.ATyy;
        xyfav{i}.w=0;
    end;
    %% updating best and ub
    best.x=best.x/rhofactor;
    best.Ax=best.Ax/rhofactor;
    if data.fit==inf,
        ub=max(abs(best.Ax-rho*data.b))-rho*cntr.delta;
    else
        ub=norm(best.Ax-rho*data.b)-rho*cntr.delta;
    end;
    %% starting point for the stage
    if (Stage==1),
        if scale.type(1)=='e',
            zn(1:dim.x)=1/dim.x;
        else
            zn(1:dim.x)=0;
        end;
        if data.fit==2,
            zn(dim.x+1:dim.x+dim.y)=0;
        elseif scale.type(1)=='e',
            zn(dim.x+1:end)=1/dim.y;
        else
            zn(dim.x+1:end)=0;
        end;
    elseif cntr.restarttype(1)=='l'
        if scale.type(1)=='e'
            zn(1:dim.x)=restartf*zn(1:dim.x)+(1-restartf)/dim.x;
        else
            zn(1:dim.x)=restartf*zn(1:dim.x);
        end;
        if data.fit==2,
            zn(dim.x+1:end)=restartf*zn(dim.x+1:end);
        elseif scale.type(1)=='e',
            zn(dim.x+1:end)=restartf*zn(dim.x+1:end)...
                +(1-restartf)/dim.y;
        else
            zn(dim.x+1:end)=restartf*zn(dim.x+1:end);
        end;
    else
        if scale.type(1)=='e'
            zn(1:dim.x)=restartf*Best(1:dim.x)/rhofactor+(1-restartf)/dim.x;
        else
            zn(1:dim.x)=restartf*Best(1:dim.x)/rhofactor;
        end;
        if data.fit==2,
            zn(dim.x+1:end)=restartf*Best(dim.x+1:end);
        elseif scale.type(1)=='e',
            zn(dim.x+1:end)=restartf*Best(dim.x+1:end)/rhofactor...
                +(1-restartf)/dim.y;
        else
            zn(dim.x+1:end)=restartf*Best(dim.x+1:end)/rhofactor;
        end;
    end;
    for StepLoc=1:cntr.StepMax,
        Steps=Steps+1;
        iloc=0;
        w=zn;
        z=zn;
        if StepLoc==1,
            [omega,domega]=GetOmegaMostNew(w,dim,scale,data.fit);
        end;
        [Phi,xyf]=GetPhiGBSPRelease(w,dim,data,rho,cntr.dgf);
        CallsLoc=CallsLoc+1;
        Calls=Calls+1;
        %%%
        % ubo=ub;
        % lbo=lb;
        [ub,lb,best,rhoub,Best]=UpdateBoundsGBSPRelease(dim,ub,lb,xyf,best,data,rho,rhoub,cntr.dgf,cntr.delta);
        %%%
        while(1==1)
            %% cutting stuff:
            %% constraint on u: bL.cf'*u+bL.c \leq 0
            if cntr.mode(1)=='a',
                if iloc==0,
                    bL.cf=zeros(dim.z,1);
                    bL.c=-inf;
                    wad=0*wad;
                    Awad=0*Awad;
                elseif iloc==1,
                    if cntr.dgf(1)=='e',
                        if data.fit==2,
                            bL.cf=[xyf.ATyy;-xyf.ATyy;rho*data.b-xyf.Axx];
                        else
                            bL.cf=[xyf.ATyy;-xyf.ATyy;rho*data.b-xyf.Axx;-rho*data.b+xyf.Axx];
                        end;
                    else
                        bL.cf=[xyf.ATyy;rho*data.b-xyf.Axx];
                    end;
                    bL.c=-bL.cf'*w;
                    wad=wf;
                    Awad=Awf;
                end;
            elseif cntr.mode(1)=='b',
                if cntr.dgf(1)=='e'
                    if data.fit==2,
                        bL.cf=[xyf.ATyy;-xyf.ATyy;rho*data.b-xyf.Axx];
                    else
                        bL.cf=[xyf.ATyy;-xyf.ATyy;rho*data.b-xyf.Axx;-rho*data.b+xyf.Axx];
                    end;
                else
                    bL.cf=[xyf.ATyy;rho*data.b-xyf.Axx];
                end;
                bL.c=-bL.cf'*w(1:dim.x+dim.y);
                wad=wf;
                Awad=Awf;
            else
                bL.cf=zeros(dim.z,1);
                bL.c=-inf;
                wad=0*wad;
                Awad=0*Awad;
            end;
            %%%%%%%%%%%%
            iloc=iloc+1;
            [zn,omegan,domegan,lambda,dB]=ProxGBSPMostNew(gamma*Phi-domega,bL,dim,scale,data.fit);
            %%%%%%%%%%%%
            nBsteps=nBsteps+dB;
            nProx=nProx+1;
            lhs=gamma*Phi'*(w-zn);
            pos=omegan-omega-domega'*(zn-z);
            if pos<-1.e-10,
                fprintf(1,'This should be positive: %.3e\n',pos);
            end;
            if (lhs<=pos+1.e-18+cntr.postolerance*gamma*rho*cntr.eps)
                break;
            end;
            w=zn;
            [Phi,xyf]=GetPhiGBSPRelease(w,dim,data,rho,cntr.dgf);
            CallsLoc=CallsLoc+1;
            Calls=Calls+1;
            wf(1:dim.xf)=xyf.x/xyf.w;
            wf(dim.xf+1:dim.xf+dim.yf)=xyf.y/xyf.w;
            Awf(1:dim.xf)=xyf.ATy/xyf.w;
            Awf(dim.xf+1:end)=xyf.Ax/xyf.w;
            %%%
            % ubo=ub;
            % lbo=lb;
            [ub,lb,best,rhoub,Best]=UpdateBoundsGBSPRelease(dim,ub,lb,xyf,best,data,rho,rhoub,cntr.dgf,cntr.delta);
            %%%
            if iloc>2,
                gamma=max(cntr.gamma_down*gamma,gamma_safe);
            end;
            %%
            if iloc>100,
                error('too long inner loop');
            end;
        end;
        %%
        omega=omegan;
        domega=domegan;
        %%%
        sumgammafact=sumgammafact+gamma*xyf.w+lambda;
        sumgammabas=sumgammabas+gamma;
        uboc=inf;
        lboc=-inf;
        for i=1:2*cntr.nAverages,
            if (StepLoc<Starts(i))&&(i<=cntr.nAverages),
                continue;
            end;
            if i<=cntr.nAverages,
                fctr=1;
            else
                fctr=StepLoc^(i-cntr.nAverages+1);
            end;
            dnm=xyfav{i}.w+fctr*lambda;
            xyfav{i}.xx=xyfav{i}.x/dnm+(fctr/dnm)*lambda*wad(1:dim.xf);
            xyfav{i}.x=xyfav{i}.x+fctr*gamma*xyf.x+fctr*lambda*wad(1:dim.xf);
            xyfav{i}.Axx=xyfav{i}.Ax/dnm+(fctr/dnm)*lambda*Awad(dim.xf+1:end);
            xyfav{i}.Ax=xyfav{i}.Ax+fctr*gamma*xyf.Ax+fctr*lambda*Awad(dim.xf+1:end);
            xyfav{i}.yy=xyfav{i}.y/dnm+(fctr/dnm)*lambda*wad(dim.xf+1:end);
            xyfav{i}.y=xyfav{i}.y+fctr*gamma*xyf.y+fctr*lambda*wad(dim.xf+1:end);
            xyfav{i}.ATyy=xyfav{i}.ATy/dnm+(fctr/dnm)*lambda*Awad(1:dim.xf);
            xyfav{i}.ATy=xyfav{i}.ATy+fctr*gamma*xyf.ATy+fctr*lambda*Awad(1:dim.xf);
            xyfav{i}.w=xyfav{i}.w+fctr*gamma*xyf.w+fctr*lambda;
            % ubo=ub;
            % lbo=lb;
            [ub,lb,best,rhoub,Best]=UpdateBoundsGBSPRelease(dim,ub,lb,xyfav{i},best,data,rho,rhoub,cntr.dgf,cntr.delta);
            if (ub<uboc)||(lb>lboc)
                ibst=i;
            end;
            uboc=ub;
            lboc=lb;
        end;
        Scores(ibst)=Scores(ibst)+1;
        %%%
        if flagrho,
            errf=ub/rho;
        else
            errf=(ub-lb)/rho;
        end;
        while(errf<=errortarget),
            errortarget=errortarget*decrement;
        end;
        if (flagrho==0)&&(ub<=lb+tolerance),
            stop=1;
        end;
        if (flagrho==1)&&(ub<=tolerance),
            stop=1;
        end;
        if cntr.printlevel,
            if (mod(Steps,cntr.printlevel)==1)||stop,
                ireport=ireport+1;
                sUB=MySprintf(ub);
                sLB=MySprintf(lb);
                target=tolerance;
                tused=cputime-tstart;
                fprintf(1,'    %8.2f: UB=%s LB=%s Steps=%d/%d Calls=%d/%d\n',...
                    tused,sUB,sLB,StepLoc,Steps,CallsLoc,Calls);
                if flagrho,
                    if (lb>0.01*ub),
                        fprintf(1,'              UB/LB-kappa=%4.2f UB/Target=%4.2f\n',...
                            ub/lb-cntr.kappa,ub/target);
                    else
                        fprintf(1,'              UB/Target=%4.2f\n',ub/target);
                    end;
                end;
                pause(0.0001);
                if stop,
                    break;
                end;
            end;
        end;
        %%%
        
        if iloc<=2,
            gamma=gamma*cntr.gamma_up;
            if gamma/gamma_safe>ACCMAX,
                gamma=ACCMAX*gamma_safe;
            end;
        end;
        if lb<=0,
            ratio=inf;
        else
            ratio=ub/lb;
        end;
        if (ratio<=1+cntr.kappa)&&flagrho,
            break;
        end;
        tused=cputime-tstart;
        if tused>cntr.cpumax,
            break;
        end;
    end;
    if stop,
        break;
    end;
    tused=cputime-tstart;
    if tused>cntr.cpumax,
        break;
    end;
    rhofactor=(rho/rhoub)*(1+cntr.omega);
    rho=rho/rhofactor;
    if cntr.accuracy(1)=='a',
        tolerance=cntr.eps*rho;
    end;
    lb=0;
end;
if stop,
    sol.ok=0;
else
    sol.ok=1;
end;
tused=cputime-tstart;
% if cntr.printlevel,
%     sUB=MySprintf(ub);
%     sLB=MySprintf(lb);
%     srho=MySprintf(rho);
% end;
sol.rho=rho;
sol.x=best.x/rho;
sol.Ax=best.Ax/rho;
sol.y=best.y;
sol.ATy=best.ATy;
sol.ell1norm=sum(abs(sol.x));
if flagrho,
    sol.opt_lwb=1/rhoub;
else
    sol.opt_lwb=0;
end;
sol.res_lwb=(lb+rho*cntr.delta)/rho;
sol.p=data.fit;
if data.fit==inf,
    sol.res=max(abs(sol.Ax-data.b));
    sol.res_rel=sol.res/max(abs(data.b));
else
    sol.res=norm(sol.Ax-data.b);
    sol.res_rel=sol.res/norm(data.b);
end;
if flagrho,
    sol.rhofixed='off';
else
    sol.rhofixed='on';
end;
if flagrho,
    sol.err.R=sol.ell1norm*cntr.omega/(1+cntr.omega);
    sol.err.abs=ub/sol.rho;
    sol.err.rel=sol.err.abs/(data.L*sol.opt_lwb);
else
    sol.err.R=0;
    sol.err.abs=ub/sol.rho;
    sol.err.rel=inf;
end;
sol.Steps=Steps;
sol.Calls=Calls;
sol.Stages=Stage;
sol.cpu=tused;
sol.control.omega=cntr.omega;
sol.control.ItrMax=cntr.StepMax;
sol.control.cpumax=cntr.cpumax;
sol.control.kappa=cntr.kappa;
sol.control.accuracy=cntr.accuracy;
sol.control.eps=cntr.eps;
sol.control.dgf=cntr.dgf;
sol.control.rhofixed=cntr.rhofixed;
sol.control.rho=cntr.rho;
sol.control.printlevel=cntr.printlevel;
