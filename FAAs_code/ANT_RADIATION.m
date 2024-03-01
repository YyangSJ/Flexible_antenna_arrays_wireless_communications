function AE_linear = ANT_RADIATION(theta,phi)
G_mE=8;
phi3db=65*pi/180;
theta3db=65*pi/180;
Am=30;
SLA=30;
AH=-min([12*(phi/phi3db)^2,Am]);
AV=-min([12*((theta-pi/2)/theta3db)^2,SLA]);
AE=G_mE-min([-AH-AV,Am]);
AE_linear=10^(AE/10);
end

