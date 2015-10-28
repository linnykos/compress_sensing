el1sparse = csvread('el1sparse_size.csv'); 
l1ls = csvread('l1ls_res_size.csv'); 
fhtp = csvread('fhtp_res_size.csv'); 
mirror = csvread('mirrorprox_res_size.csv'); 
el1 = csvread('el1_size.csv');

loqo_el1sparse = csvread('res_loqo_el1sparse_size.csv'); 
loqo_el1 = csvread('res_loqo_el1_size.csv'); 


spar = [102+(0:15)*20];
spar2 = [102,202,302,402];
len1 = length(spar);
len2 = length(spar2);

avg_time_sparse = zeros(len1,1); %simplex kcs
sd_time_sparse = zeros(len1,1);
avg_time = zeros(len1,1); %simplex
sd_time = zeros(len1,1);
avg_time_l1ls = zeros(len1,1); %l1ls 
sd_time_l1ls = zeros(len1,1);
avg_time_fhtp = zeros(len1,1); %fhtp
sd_time_fhtp = zeros(len1,1);
avg_time_mirror = zeros(len1,1); %mirror
sd_time_mirror = zeros(len1,1);


avg_time_loqosparse = zeros(len2,1); %ipm kcs
sd_time_loqosparse = zeros(len2,1);
avg_time_loqo = zeros(len2,1); %ipm
sd_time_loqo = zeros(len2,1);



counter_10 = 1;
for i=1:16
    avg_time_sparse(i) = mean(el1sparse(counter_10:(counter_10+9),2));
    sd_time_sparse(i) = std(el1sparse(counter_10:(counter_10+9),2));
       
    avg_time(i) = mean(el1(counter_10:(counter_10+9),2));
sd_time(i) = std(el1(counter_10:(counter_10+9),2));

avg_time_l1ls(i) = mean(l1ls(counter_10:(counter_10+9),2));
sd_time_l1ls(i) = std(l1ls(counter_10:(counter_10+9),2));

avg_time_fhtp(i) = mean(fhtp(counter_10:(counter_10+9),2));
sd_time_fhtp(i) = std(fhtp(counter_10:(counter_10+9),2));

avg_time_mirror(i) = mean(mirror(counter_10:(counter_10+9),2));
sd_time_mirror(i) = std(mirror(counter_10:(counter_10+9),2));


    counter_10 = counter_10+9;
end

counter_2 = 1;
for i=1:len2
    avg_time_loqosparse(i) = mean(loqo_el1sparse(counter_2:(counter_2+1),2));
    sd_time_loqosparse(i) = std(loqo_el1sparse(counter_2:(counter_2+1),2));
    
    
    avg_time_loqo(i) = mean(loqo_el1(counter_2:(counter_2+1),2));
    sd_time_loqo(i) = std(loqo_el1(counter_2:(counter_2+1),2));
   
    counter_2 = counter_2+2;
end

fig = figure(1);
set(fig, 'position', [100 100 1200 900]);
semilogy(spar,avg_time_sparse,'bv',...
    spar,avg_time,'b^',...
spar2,avg_time_loqosparse,'gv',spar2,avg_time_loqo,'g^',...
spar,avg_time_l1ls,'rv',...
spar,avg_time_fhtp,'kv',...
	spar,avg_time_mirror,'mv',...
        'LineWidth',2);
hold on;
semilogy(spar,avg_time_sparse,'b-',...
    spar,avg_time,'b-',...
spar2,avg_time_loqosparse,'gv',spar2,avg_time_loqo,'g^',...
spar,avg_time_l1ls,'r-',...
spar,avg_time_fhtp,'k-',...
	spar,avg_time_mirror,'m-',...
'LineWidth',2);
errorbar(spar,avg_time_sparse,sd_time_sparse,'b');
errorbar(spar,avg_time,sd_time,'b');
errorbar(spar2,avg_time_loqosparse,sd_time_loqosparse,'g');
errorbar(spar2,avg_time_loqo,sd_time_loqo,'g');
errorbar(spar,avg_time_l1ls,sd_time_l1ls,'r');
errorbar(spar,avg_time_fhtp,sd_time_fhtp,'k');
errorbar(spar,avg_time_mirror,sd_time_mirror,'m');
title('Solution Time vs. Sparsity');
xlabel('Size (n_2)');
ylabel('Solution Time (seconds)');
legend('Simplex KCS',...
     'Simplex',...
   'IPM KCS',...
    'IPM',...
    'L_1\_L_S',...
	'FHTP Oracle',...
	'Mirror Prox',...
		'Location','SouthEast');
hold off;
print(fig,'-dpng','plot_time_size.png');
close;
quit();
