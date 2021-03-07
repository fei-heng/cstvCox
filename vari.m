function  [beta, var] = vari(nsamp,ncovt,time,covar,covart,deltacs,wkernk,beta,wipw)
% estimate the variance of time-vary coefficient in the Cox model (ipw-c)


U=zeros(ncovt, 1);
F=zeros(ncovt, ncovt);
FU=zeros(ncovt, ncovt);

% bz = covar * beta;
bzt = zeros(nsamp, nsamp);
for i=1:nsamp
    bzt(:,i) = covart(:,:,i) * beta;
end

wt = wkernk.*wipw.*deltacs;
s0 = zeros(nsamp,1);
s1 = zeros(ncovt, nsamp);
s2 = zeros(ncovt, ncovt, nsamp);

for i=1:nsamp
    if wt(i) > 0
        for l=1:nsamp
            if time(l) >= time(i) 
                s0(i) = s0(i) + exp(bzt(l,i))*wipw(l);
                for k=1:ncovt
                    s1(k,i) = s1(k,i) + exp(bzt(l,i))*covart(l,k,i)*wipw(l);
                    for j=1:ncovt
                        s2(j,k,i) = s2(j,k,i) + exp(bzt(l,i))*covart(l,j,i)*covart(l,k,i)*wipw(l);
                    end
                end
            end
        end
        
        % ignore the derivative part for beta
        for k=1:ncovt
            U(k) = U(k)+wt(i)*(covar(i,k)-s1(k,i)/s0(i));
            for j=1:ncovt
                F(j,k) = F(j,k)+wt(i)*(s2(j,k,i)/s0(i)-s1(j,i)*s1(k,i)/(s0(i))^2);
                FU(j,k) = FU(j,k)+(wt(i))^2*(covar(i,j)-s1(j,i)/s0(i))*(covar(i,k)-s1(k,i)/s0(i));
            end
        end    
    end
end




try
    ncov = ncovt/2;
    var = pinv(F(1:ncov, 1:ncov))*FU(1:ncov, 1:ncov)*pinv(F(1:ncov, 1:ncov)); 
catch
    var = nan(ncov,ncov);
    beta = nan(ncov*2,1);
end

% try
%     ncov = ncovt/2;
%     var = zeros(ncov,ncov);
%     for k=1:ncov
%         for j=1:ncov
%             var(k,j) = pinv(F(k,j))*FU(k,j)*pinv(F(k,j));
%         end
%     end
% catch
%     var = nan(ncov,ncov);
%     beta = nan(ncov*2,1);
% end
        
        
        








