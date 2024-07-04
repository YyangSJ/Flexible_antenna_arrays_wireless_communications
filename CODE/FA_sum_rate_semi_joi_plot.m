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
    
    % for i = 1:L
    %     for j = 1:K
    %         if rand < 0.5
    %             randomMatrix(i, j) = rand * (pi/2);
    %         else
    %             randomMatrix(i, j) = 3*pi/2 + rand * (pi/2);
    %         end
    %     end
    % end
    %     phi1_u1 = randomMatrix;
    %
    %     phi1=[phi1_u1,unifrnd(pi/3-pi/6,pi+pi/6,L,K),unifrnd(pi-pi/6,5/3*pi+pi/6,L,K)];
    %     phi2=[phi1_u1,unifrnd(pi/3-pi/6,pi+pi/6,L,K),unifrnd(pi-pi/6,5/3*pi+pi/6,L,K)];
    %     phi3=[phi1_u1,unifrnd(pi/3-pi/6,pi+pi/6,L,K),unifrnd(pi-pi/6,5/3*pi+pi/6,L,K)];
    
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
        vars = [
            optimizableVariable('psi1',[-pi/2,pi/2],'Type','real'),
            optimizableVariable('psi2',[-pi/2,pi/2],'Type','real'),
            optimizableVariable('psi3',[-pi/2,pi/2],'Type','real')
            ];
        fun = @(x)sum_rate_semi_joint_bending(x.psi1,x.psi2,x.psi3,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,0); %Create Function Handles
        psi_ini = [1e-6, 1e-6, 1e-6];
        T_psi_ini = array2table(psi_ini);
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_bending_omni(reali,p)=-results.MinObjective;
        
        fun = @(x)sum_rate_semi_joint_bending(x.psi1,x.psi2,x.psi3,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,1); %Create Function Handles
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_bending_ka1(reali,p)=-results.MinObjective;
        
        
        fun = @(x)sum_rate_semi_joint_bending(x.psi1,x.psi2,x.psi3,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,2); %Create Function Handles
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_bending_ka2(reali,p)=-results.MinObjective;
        
        Sum_rate_FIX_omni(reali,p)=-sum_rate_semi_joint_bending(1e-6,1e-6,1e-6,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,0)
        Sum_rate_FIX_ka1(reali,p)=-sum_rate_semi_joint_bending(1e-6,1e-6,1e-6,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,1)
        Sum_rate_FIX_ka2(reali,p)=-sum_rate_semi_joint_bending(1e-6,1e-6,1e-6,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,2)
        
        vars = [
            optimizableVariable('psi1',[-pi/6,pi/6],'Type','real'),
            optimizableVariable('psi2',[-pi/6,pi/6],'Type','real'),
            optimizableVariable('psi3',[-pi/6,pi/6],'Type','real') ];
        fun =  @(x)sum_rate_semi_joint_rotating(x.psi1,x.psi2,x.psi3,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,0);
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_rotating_omni(reali,p)=-results.MinObjective;
        
        fun =  @(x)sum_rate_semi_joint_rotating(x.psi1,x.psi2,x.psi3,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,1);
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_rotating_ka1(reali,p)=-results.MinObjective;
        fun =  @(x)sum_rate_semi_joint_rotating(x.psi1,x.psi2,x.psi3,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,2);
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_rotating_ka2(reali,p)=-results.MinObjective;
        
        
        vars = [
            optimizableVariable('psi1',[-pi/6,pi/6],'Type','real'),
            optimizableVariable('psi2',[-pi/6,pi/6],'Type','real'),
            optimizableVariable('psi3',[-pi/6,pi/6],'Type','real') ];
        fun =  @(x)sum_rate_semi_joint_folding(x.psi1,x.psi2,x.psi3,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,0,0);
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_folding_omni(reali,p)=-results.MinObjective;
        
        fun =  @(x)sum_rate_semi_joint_folding(x.psi1,x.psi2,x.psi3,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,1);
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_folding_ka1(reali,p)=-results.MinObjective;
        fun =  @(x)sum_rate_semi_joint_folding(x.psi1,x.psi2,x.psi3,x_UPA,y_UPA,z_UPA,Nh,Nv,K,L,phi1,phi2,phi3,theta1,theta2,theta3,alpha1,alpha2,alpha3,lambda,P,1,2);
        results = bayesopt(fun,vars,'InitialX',T_psi_ini,'MaxObjectiveEvaluations',100,'ParallelMethod','clipped-model-prediction','PlotFcn',[]);
        Sum_rate_FAA_folding_ka2(reali,p)=-results.MinObjective;
        
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
 

