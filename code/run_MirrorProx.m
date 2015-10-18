sparsity_vect = [2*(1:50),125,150];
TOTALITER = 10;

LEN = length(sparsity_vect);
LEN_GAUSS = 100000;
LEN_IDX = 1500;

cntr.printlevel = 1;

n1 = 33;
n2 = 34;
d1 = 141;
d2 = 142;

l2diff = csvread('l2diff.csv');
rng_gauss = csvread('rng_gauss.csv');
rng_idx = csvread('rng_idx.csv');
counterVect = 10000*(0:(TOTALITER-1));
counterVect2 = 150*(0:(TOTALITER-1));

res_sparsity = -ones(LEN*TOTALITER,1);
res_timer = -ones(LEN*TOTALITER,1);
res_maxerror = -ones(LEN*TOTALITER,1);
res_l1error = -ones(LEN*TOTALITER,1);

res_counter = 1;
for i = 1:LEN
    for j = 1:TOTALITER
        fprintf('NEW PROBLEM: Sparisity %d, Iteration %d\n',sparsity_vect(i),j); 
        counter = counterVect2(j)+1;
        x0 = zeros(d1*d2,1);
        while sum(x0)< sparsity_vect(i)
            x0(rng_idx(counter)) = 1;
            counter = counter+1;
        end
        A = reshape(rng_gauss((counterVect(j)+1):(counterVect(j)+n1*d1)),n1,d1)';
        A = A';
        B = reshape(rng_gauss((counterVect(j)+n1*d1+1):(counterVect(j)+n1*d1+n2*d2)),n2,d2)';
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
		data.xdim = d1*d2;
		data.L = max(max(abs(U)));
		data.b = y;	
		data.delta = l2diff(j+(i-1)*TOTALITER);
		data.Mult=@MultPlain;
                data.fit = 2;

        fprintf('starting now\n===============================\n');
        sol=MirrorProx(data,cntr);
        res_timer(res_counter) = sol.cpu;        
        res_sparsity(res_counter) = sparsity_vect(i);
        res_maxerror(res_counter) = max(abs(sol.x-x0));
        res_l1error(res_counter) = sum(abs(sol.x-x0))/sparsity_vect(i);

        fprintf('Betahat sum: %f\n, Max error: %f\n, L1 error: %f\n, Solution time (seconds): %f\n',sum(sol.x),res_maxerror(res_counter),res_l1error(res_counter),res_timer(res_counter));
        res_counter = res_counter+1;
    end

res_matrix = [res_sparsity'; res_timer'; res_maxerror'; res_l1error']';
csvwrite('mirrorprox_res.csv',res_matrix);

end

res_matrix = [res_sparsity'; res_timer'; res_maxerror'; res_l1error']';
csvwrite('mirrorprox_res.csv',res_matrix);

quit();
