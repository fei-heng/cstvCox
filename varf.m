function  [beta, var] = varf(nsamp,ncov,time,covar2,deltacs,wkernk,beta)
% estimate the variance of time-vary coefficient in the Cox model (full)

U=zeros(ncov, 1);
F=zeros(ncov, ncov);

bz = zeros(nsamp,1);
for i=1:nsamp
    for j=1:ncov
        bz(i) = covar2(i,j) * beta(j);
    end
end

wt = wkernk.*deltacs;
s0 = zeros(nsamp,1);
s1 = zeros(ncov, nsamp);
s2 = zeros(ncov, ncov, nsamp);

for i=1:nsamp
    if wt(i) > 0
        for l=1:nsamp
            if time(l) >= time(i)
                s0(i) = s0(i) + exp(bz(l));
                for k=1:ncov
                    s1(k,i) = s1(k,i) + exp(bz(l))*covar2(l,k);
                    for j=1:ncov
                        s2(j,k,i) = s2(j,k,i) + exp(bz(l))*covar2(l,j)*covar2(l,k);
                    end
                end
            end
        end
        
        for k=1:ncov
            U(k) = U(k)+wt(i)*(covar2(i,k)-s1(k,i)/s0(i));
            for j=1:ncov
                F(j,k) = F(j,k)+wt(i)*(s2(j,k,i)/s0(i)-s1(j,i)*s1(k,i)/(s0(i))^2);
            end
        end
    end
end

try
    var = pinv(F);
catch
    var = nan(ncov,ncov);
    beta = nan(ncov*2,1);
end
        
        








