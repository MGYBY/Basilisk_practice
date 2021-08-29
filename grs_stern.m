% five-valued group for CD
cd = [36.0219 38.888 40.3493 40.8908 41.10471];

% refinement ratio
rr = 2.0;

pk = zeros(1,3);
cd_extrap = zeros(1,3);
fe = zeros(1,3);
cd_extrap_last = 0.0;
fe_mod = zeros(1,5);

% three groups, starts from the first value
for k=1:3
    pk(k) = (1.0/log(rr))*(log((cd(k+1)-cd(k))/(cd(k+2)-cd(k+1))));
    cd_extrap(k) = (rr^(pk(k))*cd(k+2)-cd(k+1))/(rr^(pk(k))-1);
    fe(k) = abs(cd(k+1)-cd_extrap(k))/cd_extrap(k);
    if(k==3)
        cd_extrap_last = cd_extrap(k);
    end
end

for k=1:5
    fe_mod(k) = abs(cd(k)-cd_extrap_last)/cd_extrap_last;
end

csvwrite('orig_fr.csv',[pk', cd_extrap', fe'])
csvwrite('mod_fr.csv',[cd', fe_mod'])