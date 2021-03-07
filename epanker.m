function  [kh] = epanker(tk, tvalue, hband)
% Calculate (1/h)*K((t_k-t)/h), K is Epanechnikov kernel
% K(u)=3/4(1-u^2), -1<u<1
u = (tk-tvalue)/hband;
if abs(u)<1
    kh = 3/4/hband*(1-u^2);
else
    kh = 0;
end
