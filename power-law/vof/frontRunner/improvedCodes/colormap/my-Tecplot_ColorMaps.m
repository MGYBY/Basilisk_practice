clear;
close all;
clc; clear;

% destdir = '~/'; % linux
destdir = '.\';

cmm = load('Spectral.txt');
cm = colormap(cmm);
map = cm;
colorbar
mapname = "Spectral";
%     title(char(mapname(i)))
title(char(mapname))
drawnow;
pause(1);
outname = sprintf('%s/%s.map',destdir,char(mapname));
fid = fopen(outname,'w');
%     if (fid == -1)
%         error('cant open map file')
%     end
[map_row, map_col] = size(map);
num_control_points = round(map_row/10);
total_cp = num_control_points+1;
fprintf(fid, '#!MC 900\n$!COLORMAP\n  CONTOURCOLORMAP = USERDEF\n');
fprintf(fid, '$!COLORMAPCONTROL RESETTOFACTORY\n$!COLORMAP\n  USERDEFINED\n');
fprintf(fid, '    {\n    NUMCONTROLPOINTS = %d\n', total_cp);
    
% first cp
fprintf(fid, '    CONTROLPOINT %d\n', 1);
fprintf(fid, '      {\n      COLORMAPFRACTION = %f\n',0.0);
fprintf(fid, '      LEADRGB\n');
fprintf(fid, '        {\n');
fprintf(fid, '        R = %d\n',round(map(1,1)*255));
fprintf(fid, '        G = %d\n',round(map(1,2)*255));
fprintf(fid, '        B = %d\n',round(map(1,3)*255));
fprintf(fid, '        }\n');
fprintf(fid, '      TRAILRGB\n');
fprintf(fid, '        {\n');
fprintf(fid, '        R = %d\n',round(map(1,1)*255));
fprintf(fid, '        G = %d\n',round(map(1,2)*255));
fprintf(fid, '        B = %d\n',round(map(1,3)*255));
fprintf(fid, '        }\n');
fprintf(fid, '      }\n');    
    
% middle cps
for i = 1:1:num_control_points
%         fprintf(fid, '    CONTROLPOINT %d\n',round((i-1)/10)+1);
    if round(map_row*(i/num_control_points)) < 1
        current_row = 1;
    else
    current_row = round(map_row*(i/num_control_points));
    end
    fprintf(fid, '    CONTROLPOINT %d\n',i+1);
    fprintf(fid, '      {\n      COLORMAPFRACTION = %f\n',i/num_control_points);
    fprintf(fid, '      LEADRGB\n');
    fprintf(fid, '        {\n');
    fprintf(fid, '        R = %d\n',round(map(current_row,1)*255));
    fprintf(fid, '        G = %d\n',round(map(current_row,2)*255));
    fprintf(fid, '        B = %d\n',round(map(current_row,3)*255));
    fprintf(fid, '        }\n');
    fprintf(fid, '      TRAILRGB\n');
    fprintf(fid, '        {\n');
    fprintf(fid, '        R = %d\n',round(map(current_row,1)*255));
    fprintf(fid, '        G = %d\n',round(map(current_row,2)*255));
    fprintf(fid, '        B = %d\n',round(map(current_row,3)*255));
    fprintf(fid, '        }\n');
    fprintf(fid, '      }\n');
end
    
    % the last point
%     fprintf(fid, '    CONTROLPOINT %d\n', total_cp);
%     fprintf(fid, '      {\n      COLORMAPFRACTION = %f\n',1);
%     fprintf(fid, '      LEADRGB\n');
%     fprintf(fid, '        {\n');
%     fprintf(fid, '        R = %d\n',round(map(map_row,1)*255));
%     fprintf(fid, '        G = %d\n',round(map(map_row,2)*255));
%     fprintf(fid, '        B = %d\n',round(map(map_row,3)*255));
%     fprintf(fid, '        }\n');
%     fprintf(fid, '      TRAILRGB\n');
%     fprintf(fid, '        {\n');
%     fprintf(fid, '        R = %d\n',round(map(map_row,1)*255));
%     fprintf(fid, '        G = %d\n',round(map(map_row,2)*255));
%     fprintf(fid, '        B = %d\n',round(map(map_row,3)*255));
%     fprintf(fid, '        }\n');
%     fprintf(fid, '      }\n');

fprintf(fid, '    }\n');

fclose(fid);



% end