d2 = 102+(0:15)*20;
TOTALITER = 10;
cntr.printlevel = 1;

LEN = length(d2);

n1 = 33;
n2 = 34;
d1 = 141;

l2diff = csvread('l2diff_size.csv');
rng_gauss = csvread('rng_gauss.csv');
rng_idx = csvread('rng_idx.csv');
sparsity_vect = 100;
counterVect = 100*(0:(TOTALITER-1));
counterVect2 = 150*(0:(TOTALITER-1));

res_sparsity = -ones(LEN*TOTALITER,1);
res_timer = -ones(LEN*TOTALITER,1);
res_maxerror = -ones(LEN*TOTALITER,1);
res_l1error = -ones(LEN*TOTALITER,1);

res_counter = 1;
for i = 1:LEN
    for j = 1:TOTALITER
        fprintf('NEW PROBLEM: Size %d, Iteration %d\n',d2(i),j); 
        counter = counterVect2(j)+1;
        x0 = zeros(d1*d2(i),1);
        while sum(x0)< sparsity_vect
            x0(mod(rng_idx(counter),d1*d2(i))+1) = 1;
            counter = counter+1;
        end
        A = reshape(rng_gauss((counterVect(j)+1):(counterVect(j)+n1*d1)),n1,d1)';
        A = A';
        B = reshape(rng_gauss((counterVect(j)+n1*d1+1):(counterVect(j)+n1*d1+n2*d2(i))),n2,d2(i))';
        B = B';
        U = kron(A,B);
        y = U*x0;
        %check:
        %X = vec2mat(x0,d2);
        %Y = A*X*B';
        %yprime = reshape(Y',n1*n2,1);
        %sum(abs(yprime-y))

		data.A = U;
		data.ydim = n1*n2;
		data.xdim = d1*d2(i);
		data.L = max(max(abs(U)));
		data.b = y;	
		data.delta = l2diff(j+(i-1)*TOTALITER);
		data.Mult=@MultPlain;
                data.fit = 2;

        fprintf('starting now\n===============================\n');
        sol=MirrorProx(data,cntr);
        res_timer(res_counter) = sol.cpu;        
        res_sparsity(res_counter) = d2(i);
        res_maxerror(res_counter) = max(abs(sol.x-x0));
        res_l1error(res_counter) = sum(abs(sol.x-x0))/sparsity_vect;

        fprintf('Betahat sum: %f\n, Size %f\n, Max error: %f\n, L1 error: %f\n, Solution time (seconds): %f\n',sum(sol.x),d2(i),res_maxerror(res_counter),res_l1error(res_counter),res_timer(res_counter));

        res_counter = res_counter+1;
    end

res_matrix = [res_sparsity'; res_timer'; res_maxerror'; res_l1error']';
csvwrite('mirrorprox_res_size.csv',res_matrix);

end

res_matrix = [res_sparsity'; res_timer'; res_maxerror'; res_l1error']';
csvwrite('mirrorprox_res_size.csv',res_matrix);

quit();
