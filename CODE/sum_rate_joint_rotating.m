function [sum_r_temp] = sum_rate_joint_rotating(psi1,psi2,psi3,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,Ant_type,ka)
d=lambda/2;

Rot1_local=[cos(psi1) -sin(psi1);sin(psi1) cos(psi1)];
Rot2_local=[cos(psi2) -sin(psi2);sin(psi2) cos(psi2)];
Rot3_local=[cos(psi3) -sin(psi3);sin(psi3) cos(psi3)];
p_FAA1=Rot1_local*[x_UPA;y_UPA];
p_FAA2=Rot2_local*[x_UPA;y_UPA];
p_FAA3=Rot3_local*[x_UPA;y_UPA];

x_FAA1=p_FAA1(1,:)+d*(Nh-1)*(1+sqrt(3)/3)/pi;
y_FAA1=p_FAA1(2,:);
z_FAA1=z_UPA;

x_FAA2=p_FAA2(1,:)+d*(Nh-1)*(1+sqrt(3)/3)/pi;
y_FAA2=p_FAA2(2,:);
z_FAA2=z_UPA;

x_FAA3=p_FAA3(1,:)+d*(Nh-1)*(1+sqrt(3)/3)/pi;
y_FAA3=p_FAA3(2,:);
z_FAA3=z_UPA;

Rot2_global=[cos(2*pi/3) -sin(2*pi/3);sin(2*pi/3) cos(2*pi/3)];
Rot3_global=[cos(4*pi/3) -sin(4*pi/3);sin(4*pi/3) cos(4*pi/3)];

p_FAA2_rot_g=Rot2_global*[x_FAA2;y_FAA2];
p_FAA3_rot_g=Rot3_global*[x_FAA3;y_FAA3];
x_FAA2=p_FAA2_rot_g(1,:);
y_FAA2=p_FAA2_rot_g(2,:);
x_FAA3=p_FAA3_rot_g(1,:);
y_FAA3=p_FAA3_rot_g(2,:);

H_FAA1=zeros(Nh*Nv,3*K);
H_FAA2=zeros(Nh*Nv,3*K);
H_FAA3=zeros(Nh*Nv,3*K);

psi_1=psi1;
psi_2=psi2+2*pi/3;
psi_3=psi3+4*pi/3;

rr=100;
for k=1:K*3
    for l=1:L
        for ny=1:Nh
            for nz=1:Nv
                if Ant_type==0
                    %                     r_cyl(ny,nz)=(x_FAA1(ny)*cos(phi1(l,k))*sin(theta1(l,k))+y_FAA1(ny)*sin(phi1(l,k))*sin(theta1(l,k))+z_FAA1(nz)*cos(theta1(l,k)));
                    %                     er_cyl(ny,nz)=exp(-1i*2*pi/lambda*r_cyl(ny,nz));
                    %                     r_cy2(ny,nz)=(x_FAA2(ny)*cos(phi2(l,k))*sin(theta2(l,k))+y_FAA2(ny)*sin(phi2(l,k))*sin(theta2(l,k))+z_FAA2(nz)*cos(theta2(l,k)));
                    %                     er_cy2(ny,nz)=exp(-1i*2*pi/lambda*r_cy2(ny,nz));
                    %                     r_cy3(ny,nz)=(x_FAA3(ny)*cos(phi3(l,k))*sin(theta3(l,k))+y_FAA3(ny)*sin(phi3(l,k))*sin(theta3(l,k))+z_FAA3(nz)*cos(theta3(l,k)));
                    %                     er_cy3(ny,nz)=exp(-1i*2*pi/lambda*r_cy3(ny,nz));
                    %
                    r_cyl(ny,nz)=rr-sqrt((x_FAA1(ny)-rr*cos(phi1(l,k))*sin(theta1(l,k)))^2+(y_FAA1(ny)-rr*sin(phi1(l,k))*sin(theta1(l,k)))^2+(z_FAA1(nz)-rr*cos(theta1(l,k)))^2);
                    er_cyl(ny,nz)=exp(-1i*2*pi/lambda*r_cyl(ny,nz));
                    r_cy2(ny,nz)=rr-sqrt((x_FAA2(ny)-rr*cos(phi2(l,k))*sin(theta2(l,k)))^2+(y_FAA2(ny)-rr*sin(phi2(l,k))*sin(theta2(l,k)))^2+(z_FAA2(nz)-rr*cos(theta2(l,k)))^2);
                    er_cy2(ny,nz)=exp(-1i*2*pi/lambda*r_cy2(ny,nz));
                    r_cy3(ny,nz)=rr-sqrt((x_FAA3(ny)-rr*cos(phi3(l,k))*sin(theta3(l,k)))^2+(y_FAA3(ny)-rr*sin(phi3(l,k))*sin(theta3(l,k)))^2+(z_FAA3(nz)-rr*cos(theta3(l,k)))^2);
                    er_cy3(ny,nz)=exp(-1i*2*pi/lambda*r_cy3(ny,nz));
                end
                if Ant_type==1
                    %                     r_cyl(ny,nz)=(x_FAA1(ny)*cos(phi1(l,k))*sin(theta1(l,k))+y_FAA1(ny)*sin(phi1(l,k))*sin(theta1(l,k))+z_FAA1(nz)*cos(theta1(l,k)));
                    %                     er_cyl(ny,nz)=exp(-1i*2*pi/lambda*r_cyl(ny,nz))*CosineAnt(theta1(l,k),phi1(l,k)-psi_1,ka);
                    %                     r_cy2(ny,nz)=(x_FAA2(ny)*cos(phi2(l,k))*sin(theta2(l,k))+y_FAA2(ny)*sin(phi2(l,k))*sin(theta2(l,k))+z_FAA2(nz)*cos(theta2(l,k)));
                    %                     er_cy2(ny,nz)=exp(-1i*2*pi/lambda*r_cy2(ny,nz))*CosineAnt(theta2(l,k),phi2(l,k)-psi_2,ka);
                    %                     r_cy3(ny,nz)=(x_FAA3(ny)*cos(phi3(l,k))*sin(theta3(l,k))+y_FAA3(ny)*sin(phi3(l,k))*sin(theta3(l,k))+z_FAA3(nz)*cos(theta3(l,k)));
                    %                     er_cy3(ny,nz)=exp(-1i*2*pi/lambda*r_cy3(ny,nz))*CosineAnt(theta3(l,k),phi3(l,k)-psi_3,ka);
                    r_cyl(ny,nz)=rr-sqrt((x_FAA1(ny)-rr*cos(phi1(l,k))*sin(theta1(l,k)))^2+(y_FAA1(ny)-rr*sin(phi1(l,k))*sin(theta1(l,k)))^2+(z_FAA1(nz)-rr*cos(theta1(l,k)))^2);
                    er_cyl(ny,nz)=exp(-1i*2*pi/lambda*r_cyl(ny,nz))*CosineAnt(theta1(l,k),phi1(l,k)-psi_1,ka);
                    r_cy2(ny,nz)=rr-sqrt((x_FAA2(ny)-rr*cos(phi2(l,k))*sin(theta2(l,k)))^2+(y_FAA2(ny)-rr*sin(phi2(l,k))*sin(theta2(l,k)))^2+(z_FAA2(nz)-rr*cos(theta2(l,k)))^2);
                    er_cy2(ny,nz)=exp(-1i*2*pi/lambda*r_cy2(ny,nz))*CosineAnt(theta2(l,k),phi2(l,k)-psi_2,ka);
                    r_cy3(ny,nz)=rr-sqrt((x_FAA3(ny)-rr*cos(phi3(l,k))*sin(theta3(l,k)))^2+(y_FAA3(ny)-rr*sin(phi3(l,k))*sin(theta3(l,k)))^2+(z_FAA3(nz)-rr*cos(theta3(l,k)))^2);
                    er_cy3(ny,nz)=exp(-1i*2*pi/lambda*r_cy3(ny,nz))*CosineAnt(theta3(l,k),phi3(l,k)-psi_3,ka);
                end
            end
        end
        r_ecylv=er_cyl(:);
        r_ecy2v=er_cy2(:);
        r_ecy3v=er_cy3(:);
        H_FAA1(:,k)=H_FAA1(:,k)+1/sqrt(L)*alpha1(k,l)*r_ecylv;
        H_FAA2(:,k)=H_FAA2(:,k)+1/sqrt(L)*alpha2(k,l)*r_ecy2v;
        H_FAA3(:,k)=H_FAA3(:,k)+1/sqrt(L)*alpha3(k,l)*r_ecy3v;
    end
end



H_FAA=[H_FAA1;H_FAA2;H_FAA3];

HI=inv(H_FAA'*H_FAA);
 
for k=1:K*3
    r_temp(k)=log2(1+P/K/3/(HI(k,k)));
end
sum_r_temp=-sum(r_temp); 
end

