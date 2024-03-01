function [p_alloc,rate,rate_woWF]= WF(H,F,K,N,sigma_2,P_max)
lambda0 = zeros(K,1);
for k=1:K
    lambda0(k,1) = abs(H(:,k)'*F(:,k))^2/sigma_2;
end
[lambda1,~] = sort(lambda0,'descend'); % sort in decreasing order
i=1;
p_alloc = zeros(K,1);
while(1)
    mu = 1/(N-i+1)*(P_max+sum(1./(lambda1(1:N-i+1,1)))); % Calculating water injection position 
    for k=1:N-i+1
      p_alloc(k,1) = mu - 1./lambda1(k);
    end
    if all(p_alloc(1:N-i+1,1)>=0)
        break
    end
    if p_alloc(N-i+1,1)<0
        p_alloc(N-i+1,1) = 0;
    end
    i = i+1;
    if i > N
        break
    end
end
rate=sum(log2(1+lambda1.*p_alloc));
rate_woWF=sum(log2(1+lambda1*P_max/K));
end

