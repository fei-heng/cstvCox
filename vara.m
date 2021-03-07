function  [beta, var, inv] = vara(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta,wipw,rhohat,delta)
% estimate the variance of time-vary coefficient in the Cox model (ipw-c)

% beta = beta_acc
% rhohat = rhohat_cc(:,ics)

U=zeros(ncovt, 1);
F=zeros(ncovt, ncovt);
FU=zeros(ncovt, ncovt);

% bz = covar * beta;
bzt = zeros(nsamp, nsamp);
for i=1:nsamp
    bzt(:,i) = covart(:,:,i) * beta;
end

wt = wkernk.*delta;
s0 = zeros(nsamp,1);
s1 = zeros(ncovt, nsamp);
s2 = zeros(ncovt, ncovt, nsamp);

for i=1:nsamp
    if wt(i)>0
    for l=1:nsamp
        if time(l) >= time(i)
            s0(i) = s0(i) + exp(bzt(l,i));
            for k=1:ncovt
                s1(k,i) = s1(k,i) + exp(bzt(l,i))*covart(l,k,i);
                for j=1:ncovt
                    s2(j,k,i) = s2(j,k,i) + exp(bzt(l,i))*covart(l,j,i)*covart(l,k,i);
                end
            end
        end
    end
    
    for k=1:ncovt
        U(k) = U(k)+wt(i)*(covar(i,k)-s1(k,i)/s0(i))*(wipw(i)*deltacs(i)+(1-wipw(i))*rhohat(i));
        for j=1:ncovt
            F(j,k) = F(j,k)+wt(i)*(s2(j,k,i)/s0(i)-s1(j,i)*s1(k,i)/s0(i)^2)...
                *(wipw(i)*deltacs(i)+(1-wipw(i))*rhohat(i));
            FU(j,k) = FU(j,k)+(wkernk(i))^2*(covar(i,j)-s1(j,i)/s0(i))*(covar(i,k)-s1(k,i)/s0(i))...
                *(wipw(i)*deltacs(i)+(1-wipw(i))*rhohat(i))^2*delta(i);
        end
    end
    end
end



try
    ncov = ncovt/2;
    inv = pinv(F(1:ncov, 1:ncov));
    var = inv*FU(1:ncov, 1:ncov)*inv; 
catch
    var = nan(ncov,ncov);
    inv = nan(ncov,ncov);
    beta = nan(ncov*2,1);
end

% try
%     ncov = ncovt/2;
%     var = zeros(ncov,ncov);
%     inv = zeros(ncov,ncov);
%     for k=1:ncov
%         for j=1:ncov
%             inv(k,j) = pinv(F(k,j));
%             var(k,j) = inv(k,j)*FU(k,j)*inv(k,j);
%         end
%     end
% catch
%     var = nan(ncov,ncov);
%     inv = nan(ncov,ncov);
%     beta = nan(ncov*2,1);
% end









