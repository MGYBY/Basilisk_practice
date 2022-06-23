clc; clear;
delete './periodic_formatted.csv'
normal_depth = 0.00798;
slope = 0.05011;
wavelength_dimless = 50.026;
wavelength = wavelength_dimless*normal_depth/slope;
l_x = 1200*normal_depth/slope;
x = 0.0; i_periodic = 1; count = 1;
xcoord = importdata('C:\Users\Admin\Desktop\xcoord_60.txt',' ');
hsol = importdata('C:\Users\Admin\Desktop\hnsol_60.txt',' ');
total_count = floor(1200/wavelength_dimless)*size(xcoord, 1);
while count <= total_count
    if mod(i_periodic, size(xcoord, 1))~=0
        xcoord_periodic(count) = xcoord(mod(i_periodic, size(xcoord, 1)));
        h_periodic(count) = hsol(mod(i_periodic, size(xcoord, 1)));
    else
        xcoord_periodic(count) = xcoord(mod(i_periodic, size(xcoord, 1))+1);
        h_periodic(count) = hsol(mod(i_periodic, size(xcoord, 1))+1);
    end
    x = x + wavelength_dimless;
    i_periodic = i_periodic + 1;
    count = count + 1;
end

xcoord_periodic_mod = linspace(0, l_x, total_count); 
res_series = [xcoord_periodic_mod', h_periodic'];
csvwrite ('./periodic_formatted.csv', res_series)