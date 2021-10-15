function [ X ] = LogSquare(D,rho)

[U,S,V] = svd(D,'econ');
sy = size(S);
S0 = diag(S);
r = length(S0);

P = [rho*ones(r,1), -1*rho*S0, (rho+1)*ones(r,1), -1*rho*S0];

rt = zeros(r,1);

for t = 1:r
    p = P(t,:);
    rts = roots(p);
    
    rts = rts(rts==real(rts));
    rts=rts(rts>=0);
    rts=[rts;0];
    L = length(rts);
    if L == 1
        rt(t) = rts;
    else
%         funval = log(1+rts.^2)+rho.*(rts-S0(t)).^2;
        funval = atan(rts.^2)+(rho/2).*(rts-S0(t)).^2;
        rttem = rts(funval==min(funval));
        rt(t) = rttem(1);
    end
end

sig = diag(rt);
for i=1:3
    S(i,:)=sig(i,:);
end
X = U*S*V';

end