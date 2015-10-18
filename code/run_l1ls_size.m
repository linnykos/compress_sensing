d2 = 102+(0:15)*20;
TOTALITER = 10;

LEN = length(d2);
LEN_GAUSS = 100000;
LEN_IDX = 150;

lambda = 0.01;      % regularization parameter
rel_tol = 10^-8;     % relative target duality gap


n1 = 33;
n2 = 34;
d1 = 141;

rng_gauss = csvread('rng_gauss.csv');
rng_idx = csvread('rng_idx.csv');
sparsity_vect = 100;
counterVect = 100*(0:(TOTALITER-1));
counterVect2 = 150*(0:(TOTALITER-1));

res_sparsity = -ones(LEN*TOTALITER,1);
res_timer = -ones(LEN*TOTALITER,1);
res_maxerror = -ones(LEN*TOTALITER,1);
res_l1error = -ones(LEN*TOTALITER,1);
res_l2error = -ones(LEN*TOTALITER,1);

res_matrix = [res_sparsity'; res_timer'; res_maxerror'; res_l1error']';
csvwrite('l1ls_res_size.csv',res_matrix);

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

        fprintf('starting now\n===============================\n');
        tic;
	[x,status] = l1_ls(U,y,lambda,rel_tol,true);
        res_timer(res_counter) = toc;
        res_sparsity(res_counter) = d2(i);
        res_maxerror(res_counter) = max(abs(x-x0));
        res_l1error(res_counter) = sum(abs(x-x0))/sparsity_vect;
	res_l2error(res_counter) = sum((y-U*x).^2);


	fprintf('Betahat sum: %f\n, Size %f\n, Max error: %f\n, L1 error: %f\n, Solution time (seconds): %f\n',sum(x),d2(i),res_maxerror(res_counter),res_l1error(res_counter),res_timer(res_counter));
        res_counter = res_counter+1;
    end

res_matrix = [res_sparsity'; res_timer'; res_maxerror'; res_l1error']';
csvwrite('l1ls_res_size.csv',res_matrix);

end

res_matrix = [res_sparsity'; res_timer'; res_maxerror'; res_l1error']';
csvwrite('l1ls_res_size.csv',res_matrix);
res_l2error = res_l2error';
csvwrite('l2diff_size.csv',res_l2error);
                                     

quit();
