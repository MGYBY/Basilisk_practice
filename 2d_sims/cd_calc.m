U = 1.03774;
H = 0.00798;
W = 0.40;
topCoord = 1.70; 
bottomCoord = 1.30;
Lx = 42.0;
MAXLEVEL = 10;
dyGauge = (topCoord - bottomCoord - (Lx / (2^MAXLEVEL))) / 6.;
dx = dyGauge;

m10 = importdata('./gf',' ');
m20 = importdata('./gb',' ');
m1 = m10(:,2:end);
m2 = m20(:,2:end);


time = m10(:,1);

delta_fw = sum((0.50*9.81*(10^(3))*dx*(m1.^2-m2.^2)), 2); % take sum by cols
normal_fw = sum((0.50*(10^(3))*W*(U^2)*H), 2); % take sum by cols

cd = delta_fw./normal_fw;
cd_time_series = [time, cd];

csvwrite ('./cd_time_series.csv', cd_time_series)
