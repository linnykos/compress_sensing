el1sparse = csvread('el1sparse.csv'); 
l1ls = csvread('l1ls_res.csv'); 
fhtp = csvread('fhtp_res.csv'); 
fhtpa = csvread('fhtp_res_agnostic.csv'); 
mirror = csvread('mirrorprox_res.csv'); 
el1 = csvread('el1.csv');
loqo_el1sparse = csvread('res_loqo_el1sparse.csv'); 
loqo_el1 = csvread('res_loqo_el1.csv'); 

spar = [2*(1:50),125,150];
spar2 = [2,10,20,30,40,50,60,70,80,90,100,125,150];
spar3 = [2,20,50,70,100,150];
len1 = length(spar);
len2 = length(spar2);
len3 = length(spar3);

avg_time_sparse = zeros(len1,1); %param matrix 
sd_time_sparse = zeros(len1,1);
avg_time = zeros(len1,1); %param vector 
sd_time = zeros(len1,1);
avg_time_l1ls = zeros(len1,1); %l1ls 
sd_time_l1ls = zeros(len1,1);
avg_time_loqosparse = zeros(len2,1); %loqo el1sparse
sd_time_loqosparse = zeros(len2,1);
avg_time_loqo = zeros(len3,1); %loqo el1
sd_time_loqo = zeros(len3,1);
avg_time_fhtp = zeros(len1,1); %fhtp
sd_time_fhtp = zeros(len1,1);
avg_time_fhtpa = zeros(len1,1); %fhtp agnostic
sd_time_fhtpa = zeros(len1,1);
avg_time_mirror = zeros(len1,1); %mirror
sd_time_mirror = zeros(len1,1);


counter_10 = 1;
for i=1:52
    avg_time_sparse(i) = mean(el1sparse(counter_10:(counter_10+9),2));
    sd_time_sparse(i) = std(el1sparse(counter_10:(counter_10+9),2));
       
    avg_time(i) = mean(el1(counter_10:(counter_10+9),2));
sd_time(i) = std(el1(counter_10:(counter_10+9),2));

avg_time_l1ls(i) = mean(l1ls(counter_10:(counter_10+9),2));
sd_time_l1ls(i) = std(l1ls(counter_10:(counter_10+9),2));

avg_time_fhtp(i) = mean(fhtp(counter_10:(counter_10+9),2));
sd_time_fhtp(i) = std(fhtp(counter_10:(counter_10+9),2));

avg_time_fhtpa(i) = mean(fhtpa(counter_10:(counter_10+9),2));
sd_time_fhtpa(i) = std(fhtpa(counter_10:(counter_10+9),2));

avg_time_mirror(i) = mean(mirror(counter_10:(counter_10+9),2));
sd_time_mirror(i) = std(mirror(counter_10:(counter_10+9),2));


    counter_10 = counter_10+10;
end


counter_2 = 1;
for i=1:len2
    avg_time_loqosparse(i) = mean(loqo_el1sparse(counter_2:(counter_2+1),2));
    sd_time_loqosparse(i) = std(loqo_el1sparse(counter_2:(counter_2+1),2));
    
    counter_2 = counter_2+2;
end

counter_2 = 1;
for i=1:len3
    avg_time_loqo(i) = mean(loqo_el1(counter_2:(counter_2+1),2));
    sd_time_loqo(i) = std(loqo_el1(counter_2:(counter_2+1),2));
   
    counter_2 = counter_2+2;
end
%%



fig = figure(1);
set(fig, 'Position', [1 1 900 600]);
semilogy(spar,avg_time_sparse,'bv',...
    spar,avg_time,'b^',...
spar2,avg_time_loqosparse,'gv',spar3,avg_time_loqo,'g^',...
spar,avg_time_l1ls,'rv',...
spar,avg_time_fhtp,'kv',...
spar,avg_time_fhtpa,'k^',...
	spar,avg_time_mirror,'mv',...
        'LineWidth',2);
hold on;
semilogy(spar,avg_time_sparse,'b-',...
    spar,avg_time,'b-',...
spar,avg_time_l1ls,'r-',...
spar2,avg_time_loqosparse,'g-',spar3,avg_time_loqo,'g-',...
spar,avg_time_fhtp,'k-',...
spar,avg_time_fhtpa,'k-',...
	spar,avg_time_mirror,'m-',...
'LineWidth',2);
errorbar(spar,avg_time_sparse,sd_time_sparse,'b');
errorbar(spar,avg_time,sd_time,'b');
errorbar(spar,avg_time_l1ls,sd_time_l1ls,'r');
errorbar(spar2,avg_time_loqosparse,sd_time_loqosparse,'g');
errorbar(spar3,avg_time_loqo,sd_time_loqo,'g');
errorbar(spar,avg_time_fhtp,sd_time_fhtp,'k');
errorbar(spar,avg_time_fhtpa,sd_time_fhtpa,'k');
errorbar(spar,avg_time_mirror,sd_time_mirror,'m');
title('Solution Time vs. Sparsity');
xlabel('Number of nonzeros in signal');
ylabel('Solution Time (seconds)');
legend('Simplex KCS',...
     'Simplex',...
    'IPM KCS',...
    'IPM',...
     'L_1\_L_S',...
	'FHTP Oracle',...
	'FHTP Agnostic',...
	'Mirror Prox',...
		'Location','SouthEast');
hold off;
print(fig,'-dpng','plot_time.png');
close;
quit();
