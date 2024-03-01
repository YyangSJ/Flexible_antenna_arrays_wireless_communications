clear all;
close all;
rng(2024)
f=10e9;lambda=3e8/f;
d=lambda/2;
Nh=4; % horizontal antenna number
Nv=4; % vertical antenna number
y_UPA=[-(Nh-1)/2:(Nh-1)/2]*d;
z_UPA=[-(Nv-1)/2:(Nv-1)/2]*d;
psi_max=pi/2; % maximum bending angle
PSI=pi/800:pi/800:psi_max; % numerical search
realization=333;
K=16 % number of users
L=12 % number of channel paths
sigma=1; % noise power
PP=sqrt(10.^((-15:5:20)/10)); % SNR
for reali=1:realization
    alpha=1/sqrt(2)*(normrnd(0,1,K,L)+1i*normrnd(0,1,K,L)); % path gain
    theta=unifrnd(0,pi,L,K);
    phi=unifrnd(-pi/2,pi/2,L,K);
    for i=1:length(PSI)
        psi=PSI(i);
        y_FAA=(Nh-1)*d/2/psi*sin(psi*y_UPA/((Nh-1)*d/2));  % y-position mapping
        x_FAA=-(Nh-1)*d/2/psi+(Nh-1)*d/2/psi*cos(psi*y_UPA/((Nh-1)*d/2));  % x-position mapping
        z_FAA=z_UPA;
        H_UPA=zeros(Nh*Nv,K);H_FAA=zeros(Nh*Nv,K);
        for k=1:K
            for l=1:L
                for ny=1:Nh
                    psi_ny=2*psi/(Nh-1)*(ny-1)-psi;
                    psi_ny_usca=pi/(Nh-1)*(ny-1)-pi/2;
                    for nz=1:Nv
                        r_upa(ny,nz)=(y_UPA(ny)*sin(phi(l,k))*sin(theta(l,k))+z_UPA(nz)*cos(theta(l,k)));
                        r_cyl(ny,nz)=(x_FAA(ny)*cos(phi(l,k))*sin(theta(l,k))+y_FAA(ny)*sin(phi(l,k))*sin(theta(l,k))+z_FAA(nz)*cos(theta(l,k)));
                        er_upa(ny,nz)=exp(-1i*2*pi/lambda*r_upa(ny,nz))*ANT_RADIATION(theta(l,k),phi(l,k));
                        er_cyl(ny,nz)=exp(-1i*2*pi/lambda*r_cyl(ny,nz))*ANT_RADIATION(theta(l,k),phi(l,k)-psi_ny);
                    end
                end
                r_eupav=er_upa(:); r_ecylv=er_cyl(:);
                H_UPA(:,k)=H_UPA(:,k)+1/sqrt(L)*alpha(k,l)*r_eupav;
                H_FAA(:,k)=H_FAA(:,k)+1/sqrt(L)*alpha(k,l)*r_ecylv;
            end
        end
        H_UPA=H_UPA';
        H_FAA=H_FAA';
        F_UPA=H_UPA'*inv(H_UPA*H_UPA'+1e-6); % ZF precoding
        F_FAA=H_FAA'*inv(H_FAA*H_FAA'+1e-6);
        
        for k=1:K
            F_UPA(:,k)= F_UPA(:,k)/norm( F_UPA(:,k));
            F_FAA(:,k)= F_FAA(:,k)/norm( F_FAA(:,k));
        end
        for p=1:length(PP)
            P=PP(p)^2;
            [~,rate_UPA(reali,p),rate_UPA_wo_WF(reali,p)]=WF(H_UPA',F_UPA,K,K,sigma,P); % WF power allocation
            [~,rate_FAA_temp(i,p),rate_FAA_temp_wo_WF(i,p)]=WF(H_FAA',F_FAA,K,K,sigma,P);
        end
    end
    [rate_FAA(reali,:),~]=max(rate_FAA_temp);
    [rate_FAA_wo_WF(reali,:),~]=max(rate_FAA_temp_wo_WF);
end
rate_UPA_mean=mean(rate_UPA,1)
rate_FAA_mean=mean(rate_FAA,1)
rate_UPA_wo_WF_mean=mean(rate_UPA_wo_WF,1)
rate_FAA_wo_WF_mean=mean(rate_FAA_wo_WF,1) 
figure
plot(-15:5:20,rate_UPA_mean,'-r')
hold on
plot(-15:5:20,rate_FAA_mean,'-black')
hold on
plot(-15:5:20,rate_UPA_wo_WF_mean,'--r')
hold on
plot(-15:5:20,rate_FAA_wo_WF_mean,'--black')
grid on
