clc; clear;
system('rm ./vertSlice');
system('rm ./gaugeLoc');
system('rm ./vert_vel_profile');

% coordinate transform
max_level = 16;
% num_slices = 2^max_level;
% for testing
num_slices = 5;
% for selected range
% num_slices = 1200; % based on single cells
xMiddleNS = 2069.5;
vert_reso = 128; % num of points to do vert. interp.
vert_interp_range = 4.251; % range for vert. interp. in terms of Nusselt Scaling
so = 0.015;
lxNS = 2400.0; % domain length in terms of Nusselt Scaling

gerris_filename = 'out-666.667.gfs';
vert_slice_filename = 'vert_vel_profile';

gerris_interp_command_head = 'gerris2D -e "OutputLocation { istep = 1 } vertSlice gaugeLoc" ';
gerris_interp_command_tail = ' > /dev/null';
gerris_interp_command = gerris_interp_command_head+" "+gerris_filename+gerris_interp_command_tail;

for k = 1:1:num_slices
    % construct interpolation coordinates for a specified location 
    gauge_mat = zeros(vert_reso,3);
    % x_coord = k*(lxNS/(2^max_level))-0.5*(lxNS/(2^max_level));
    % x_coord = k*(lxNS/num_slices)-0.5*(lxNS/num_slices);
    x_coord = xMiddleNS - 7.0/8.0*num_slices/((2^max_level)/lxNS) + k/((2^max_level)/lxNS);
    for l = 1:1:vert_reso
        y_coord = vert_interp_range/vert_reso*l;

        x_coord_trans = y_coord-0.5*lxNS;
        y_coord_trans = -1.0*(x_coord-0.5*lxNS);
        gauge_mat(l,1) = x_coord_trans;
        gauge_mat(l,2) = y_coord_trans;
        gauge_mat(l,3) = 0.0;
    end
    writematrix(gauge_mat, 'gaugeLoc','Delimiter',' ');
    system('mv gaugeLoc.txt gaugeLoc');

    % do the slice interpolation in aid of gerris2D
    system(gerris_interp_command);

    % find the max value in the slice
    % useful columns: u.z-9, f-13
    slice_matrix = readmatrix('vertSlice', 'NumHeaderLines', 1);
    vert_vel_fluid = slice_matrix(:,9).*slice_matrix(:,13);
    % vert_vel_fluid = abs(slice_matrix(:,9).*slice_matrix(:,13));
    slice_max = max(vert_vel_fluid);
    % system('rm ./vertSlice');
    % system('rm ./gaugeLoc');

    % write (append) the max vert vel to file
    fid = fopen(vert_slice_filename, 'a+');
    fprintf(fid, '%g %g \n', x_coord, slice_max);
    fclose(fid);
end