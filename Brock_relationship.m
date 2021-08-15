clc;
clear;
Fr = 3.71;
d_h_star_min = 0.001; h_star_min_end = 0.92;h_star_min_begin = 0.1;
% So_lambda_div_hn = zeros(1,1/d_h_star_min-1);
% hn_star = zeros(1,1/d_h_star_min-1);
ind = 1;len1 = length(h_star_min_begin:d_h_star_min:(h_star_min_end-d_h_star_min));
hn_star = zeros(1,len1); So_lambda_div_hn = zeros(1,len1); c_n = zeros(1,len1);
% for h_star_min = 0.2:d_h_star_min:(1-d_h_star_min)
for h_star_min = h_star_min_begin:d_h_star_min:(h_star_min_end-d_h_star_min)
    % calc. of So_lambdaStar
    ha_star = max((1/2/Fr^2)*(1+2*Fr+sqrt(1+4*Fr)),(1/2/Fr^2)*(1+2*Fr-sqrt(1+4*Fr)));
    hb_star = min((1/2/Fr^2)*(1+2*Fr+sqrt(1+4*Fr)),(1/2/Fr^2)*(1+2*Fr-sqrt(1+4*Fr)));
    h_star_max = h_star_min*(1/2*(sqrt(1+8*(1/h_star_min)^3)-1));
    A = (1+ha_star+ha_star^2)/(ha_star-hb_star);
    B = (1+hb_star+hb_star^2)/(ha_star-hb_star);
    K1 = log((h_star_max-ha_star)/(h_star_min-ha_star));
    K2 = log((h_star_max-hb_star)/(h_star_min-hb_star));
    So_lambdaStar = (h_star_max-h_star_min)+A*K1-B*K2;

    % calc. of hn_star
    h_star_av_1 = 1/2*(h_star_max^2-h_star_min^2)+(A-B)*(h_star_max-h_star_min)+A*ha_star*K1-B*hb_star*K2;
    h_star_av_2 = (h_star_max-h_star_min)+A*K1-B*K2;
    h_star_av = h_star_av_1/h_star_av_2;
    hn_star(ind) = (((1+Fr)*h_star_av-1)/Fr)^(2/3);

    So_lambda_div_hn(ind) = So_lambdaStar/hn_star(ind);
    c_n(ind) = (1+Fr)*sqrt(1/hn_star(ind));
    ind = ind+1;
end

% h_star_min = 0.9:d_h_star_min:(1-d_h_star_min);
h_star_min = h_star_min_begin:d_h_star_min:(h_star_min_end-d_h_star_min);
h_n_max = (h_star_min./hn_star(1:length(h_star_min))).*(1/2.*(sqrt(1+8*(1.0 ./h_star_min).^3)-1));
h_n_min = h_star_min./hn_star(1:length(h_star_min));

% select the physically meaningful roots
% num_real = 1;
% for k=1:length(h_n_max)
%     if imag(h_n_max(k))==0
%         num_real=num_real+1;
%     end
% end
num_real = sum(imag(h_n_max)==0);

h_n_max_real = zeros(1, num_real);
h_n_min_real = zeros(1, num_real);
So_lambda_div_hn_real = zeros(1, num_real);
c_n_real = zeros(1, num_real);
num_real = 1;
for k=1:length(h_n_max)
    if imag(h_n_max(k))==0
        h_n_max_real(num_real)=h_n_max(k);
        h_n_min_real(num_real)=h_n_min(k);
        So_lambda_div_hn_real(num_real) = So_lambda_div_hn(k);
        c_n_real(num_real) = c_n(k);
        num_real=num_real+1;
    end
end


% plot(h_star_min, So_lambda_div_hn)
scatter(So_lambda_div_hn_real(1:length(h_n_max_real)), h_n_max_real)
title('hmax/hn')
res_mat = [So_lambda_div_hn_real', h_n_max_real'];
csvwrite(strcat(num2str(Fr),'_hmax','.csv'),res_mat)
% csvwrite(strcat(num2str(t),'.csv'),coord)
figure(2)
scatter(So_lambda_div_hn_real(1:length(h_n_max_real)),c_n_real(1:length(h_n_max_real)))
title('c/gw')
res_mat1 = [So_lambda_div_hn_real', c_n_real'];
csvwrite(strcat(num2str(Fr),'_cn','.csv'),res_mat1)

figure(3)
scatter(So_lambda_div_hn_real(1:length(h_n_max_real)),h_n_min_real(1:length(h_n_max_real)))
title('hmin/hn')
res_mat2 = [So_lambda_div_hn_real', h_n_min_real'];
csvwrite(strcat(num2str(Fr),'_hmin','.csv'),res_mat2)

% figure(3)
% scatter(So_lambda_div_hn_real(1:length(h_n_max_real)),So_lambda_div_hn_real(1:length(h_n_max_real))./c_n_real(1:length(h_n_max_real)))