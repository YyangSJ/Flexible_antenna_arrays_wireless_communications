function [sum_r_temp] = sum_rate_separate_rotating(psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi,theta,alpha,lambda,P,Ant_type,ka,sector_index)
d=lambda/2;

Rot1_local=[cos(psi) -sin(psi);sin(psi) cos(psi)]; 
p_FAA=Rot1_local*[x_UPA;y_UPA]; 
x_FAA=p_FAA(1,:)+d*(Nh-1)*(1+sqrt(3)/3)/pi;
y_FAA=p_FAA(2,:);
z_FAA=z_UPA;
 
switch sector_index
    case 1
        'ok';
    case 2
        Rot=[cos(2*pi/3) -sin(2*pi/3);sin(2*pi/3) cos(2*pi/3)];
        p_FAA =Rot*[x_FAA;y_FAA];
        x_FAA =p_FAA(1,:);
        y_FAA =p_FAA(2,:);
    otherwise
        Rot=[cos(4*pi/3) -sin(4*pi/3);sin(4*pi/3) cos(4*pi/3)];
        p_FAA =Rot*[x_FAA;y_FAA];
        x_FAA =p_FAA(1,:);
        y_FAA =p_FAA(2,:);
end


H_FAA=zeros(Nh*Nv,K);

            rr=100;
for k=1:K
    for l=1:L
        for ny=1:Nh
            switch sector_index
                case 1
                    psi_ny=psi;
                case 2
                    psi_ny=psi+2*pi/3;
                otherwise
                    psi_ny=psi+4*pi/3;
            end
            for nz=1:Nv
                if Ant_type==0
                    r_cy(ny,nz)=rr-sqrt((x_FAA(ny)-rr*cos(phi(l,k))*sin(theta(l,k)))^2+(y_FAA(ny)-rr*sin(phi(l,k))*sin(theta(l,k)))^2+(z_FAA(nz)-rr*cos(theta(l,k)))^2);
                    er_cy(ny,nz)=exp(-1i*2*pi/lambda*r_cy(ny,nz));
                end
                if Ant_type==1    
                    r_cy(ny,nz)=rr-sqrt((x_FAA(ny)-rr*cos(phi(l,k))*sin(theta(l,k)))^2+(y_FAA(ny)-rr*sin(phi(l,k))*sin(theta(l,k)))^2+(z_FAA(nz)-rr*cos(theta(l,k)))^2);
                    er_cy(ny,nz)=exp(-1i*2*pi/lambda*r_cy(ny,nz))*CosineAnt(theta(l,k),phi(l,k)-psi_ny,ka);
                end
            end
        end
        r_ecyv=er_cy(:); 
        H_FAA(:,k)=H_FAA(:,k)+1/sqrt(L)*alpha(k,l)*r_ecyv;
    end
end

HI=inv(H_FAA'*H_FAA);
sum_r_temp=0;
for k=1:K
    sum_r_temp=sum_r_temp+log2(1+P/K/3/(HI(k,k)));
end
sum_r_temp=-sum_r_temp;
end

