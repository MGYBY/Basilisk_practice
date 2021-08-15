clc;
clear;
% long wave approximation
Fr  = 5.60;
h_star_min = (1/2/Fr^2)*(1+2*Fr+sqrt(1+4*Fr));
hmax_div_hmin = 1/2*(sqrt(1+(2/h_star_min)^3)-1);
hn_star = (((1+Fr)*h_star_min-1)/Fr)^(2/3);
hmax_div_hn = hmax_div_hmin*h_star_min/hn_star;
hmin_div_hn = h_star_min/hn_star;
c_div_gh = (1+Fr)/sqrt(h_star_min);
c_div_U = c_div_gh/Fr;
umax_div_U = (c_div_gh-1/hmax_div_hn*(c_div_gh-Fr))/Fr;
qmax_divUH = umax_div_U*hmax_div_hn;