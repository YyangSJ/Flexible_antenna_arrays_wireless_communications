clear all;
close all;
rng(2024)

f=10e9;lambda=3e8/f;
d=lambda/2;
Nh=8; % horizontal antenna number
Nv=2; % vertical antenna number
x_UPA=zeros(1,Nh);
y_UPA=[-(Nh-1)/2:(Nh-1)/2]*d;
z_UPA=[-(Nv-1)/2:(Nv-1)/2]*d;
psi_max=pi/2; % maximum bending angle
realization=70;
K=16; % number of users
L=5; % number of channel paths
sigma=1; % noise power
PP=10.^((-15:5:15)/10); % SNR

%Ant_type=1 ; % 0: Omni-directional pattern, 1: Cosine pattern, 2: Ideal sector pattern
for reali=1:realization
    reali
    alpha1=1/sqrt(2)*(normrnd(0,1,3*K,L)+1i*normrnd(0,1,3*K,L)); % path gain
    alpha2=1/sqrt(2)*(normrnd(0,1,3*K,L)+1i*normrnd(0,1,3*K,L)); % path gain
    alpha3=1/sqrt(2)*(normrnd(0,1,3*K,L)+1i*normrnd(0,1,3*K,L)); % path gain
    theta1=unifrnd(pi/3,pi/3*2,L,3*K);
    theta2=unifrnd(pi/3,pi/3*2,L,3*K);
    theta3=unifrnd(pi/3,pi/3*2,L,3*K);
    
    for i = 1:L
        for j = 1:K
            if rand < 0.5
                randomMatrix(i, j) = rand * (pi/3);
            else
                randomMatrix(i, j) = 5*pi/3 + rand * (pi/3);
            end
        end
    end
    phi1_u1 = randomMatrix;
    
    phi1=[phi1_u1,unifrnd(pi/3,pi,L,K),unifrnd(pi,5/3*pi,L,K)];
    phi2=[phi1_u1,unifrnd(pi/3,pi,L,K),unifrnd(pi,5/3*pi,L,K)];
    phi3=[phi1_u1,unifrnd(pi/3,pi,L,K),unifrnd(pi,5/3*pi,L,K)];
    
    for p=1:length(PP)
        P=PP(p)
        
        Sum_rate_FIX_omni(reali,p)=-sum_rate_semi_joint_bending(1e-6,1e-6,1e-6,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,0)
        Sum_rate_FIX_ka1(reali,p)=-sum_rate_semi_joint_bending(1e-6,1e-6,1e-6,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,5)
        Sum_rate_FIX_ka2(reali,p)=-sum_rate_semi_joint_bending(1e-6,1e-6,1e-6,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,6)
%         vars = [
%             optimizableVariable('psi',[-pi/2,pi/2],'Type','real')
%             ];
%         psi_ini = 1e-6;
%         T_psi_ini = array2table(psi_ini);
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,0,1,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,0,1,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,0,1,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_bending_omni(reali,p)=sum_rate_semi_joint_bending(psi1_opt,psi2_opt,psi3_opt,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,1)
%         
%         
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,1,1,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,1,1,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,1,1,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_bending_ka1(reali,p)=sum_rate_semi_joint_bending(psi1_opt,psi2_opt,psi3_opt,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,1)
%         
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,1,2,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,1,2,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_bending(x.psi,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,1,2,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_bending_ka2(reali,p)=sum_rate_semi_joint_bending(psi1_opt,psi2_opt,psi3_opt,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,2)
%         
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,0,1,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,0,1,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,0,1,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_rotating_omni(reali,p)=sum_rate_semi_joint_rotating(psi1_opt,psi2_opt,psi3_opt,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,1)
%         
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,1,1,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,1,1,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,1,1,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_rotating_ka1(reali,p)=sum_rate_semi_joint_rotating(psi1_opt,psi2_opt,psi3_opt,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,1)
%         
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,1,2,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,1,2,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_rotating(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,1,2,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_rotating_ka2(reali,p)=sum_rate_semi_joint_rotating(psi1_opt,psi2_opt,psi3_opt,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,2)
%         
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,0,1,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,0,1,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,0,1,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_folding_omni(reali,p)=sum_rate_semi_joint_folding(psi1_opt,psi2_opt,psi3_opt,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,1)
%         
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,1,1,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,1,1,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,1,1,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_folding_ka1(reali,p)=sum_rate_semi_joint_folding(psi1_opt,psi2_opt,psi3_opt,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,1)
%         
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1(:,1:K),theta1(:,1:K),alpha1(1:K,:),lambda,P,1,2,1)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi1_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi2(:,1+K:2*K),theta2(:,1+K:2*K),alpha2(1+K:2*K,:),lambda,P,1,2,2)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi2_opt=table2array(bestPoint(results));
%         fun = @(x)sum_rate_separate_folding(x.psi,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi3(:,1+2*K:3*K),theta3(:,1+2*K:3*K),alpha3(1+2*K:3*K,:),lambda,P,1,2,3)
%         results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',33,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
%         psi3_opt=table2array(bestPoint(results));
%         Sum_rate_FAA_folding_ka2(reali,p)=sum_rate_semi_joint_folding(psi1_opt,psi2_opt,psi3_opt,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,2)
%         
    end
end

Sum_rate_FIX_omni_mean= mean(Sum_rate_FIX_omni,1)
Sum_rate_FIX_ka1_mean= mean(Sum_rate_FIX_ka1,1)
Sum_rate_FIX_ka2_mean= mean(Sum_rate_FIX_ka2,1)

Sum_rate_FAA_bending_omni_mean= mean(Sum_rate_FAA_bending_omni,1)
Sum_rate_FAA_bending_ka1_mean= mean(Sum_rate_FAA_bending_ka1,1)
Sum_rate_FAA_bending_ka2_mean= mean(Sum_rate_FAA_bending_ka2,1)


Sum_rate_FAA_rotating_omni_mean= mean(Sum_rate_FAA_rotating_omni,1)
Sum_rate_FAA_rotating_ka1_mean= mean(Sum_rate_FAA_rotating_ka1,1)
Sum_rate_FAA_rotating_ka2_mean= mean(Sum_rate_FAA_rotating_ka2,1)


Sum_rate_FAA_folding_omni_mean= mean(Sum_rate_FAA_folding_omni,1)
Sum_rate_FAA_folding_ka1_mean= mean(Sum_rate_FAA_folding_ka1,1)
Sum_rate_FAA_folding_ka2_mean= mean(Sum_rate_FAA_folding_ka2,1)


