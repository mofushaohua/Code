function [x,objV] = logwshrinkObj(x,tol,sX,mode)

if ~exist('mode','var')
    % mode = 1是采用lateral slice的方法
    % mode = 2是采用front slice的方法
    % mode = 3是采用top slice的方法
    mode = 1;
end

X=reshape(x,sX);
if mode == 1
    Y=X2Yi(X,3);
elseif mode == 3
    Y=shiftdim(X, 1);
else
    Y = X;
end
% 

Yhat = fft(Y,[],3);
% weight = C./sTrueV+eps;
% weight = 1;
% tau = rho*weight;
objV = 0;
if mode == 1
    n3 = sX(2);
elseif mode == 3
    n3 = sX(1);
else
    n3 = sX(3);
end

for i = 1:n3
    [uhat,shat,vhat] = svd(full(Yhat(:,:,i)),'econ');
    shat = diag(shat);
    epson = 1e-7;
    tau = tol - shat*epson;
    diagS1 = shat - epson;
    temp_s = diagS1.^2-4*(tau); 
    temp_shat = 0.5*(diagS1+sqrt(temp_s));
    temp_shat(temp_s < 0) = 0;
    temp_shat = diag(temp_shat);
    Yhat(:,:,i) = uhat*temp_shat*vhat';
end

Y = ifft(Yhat,[],3);
if mode == 1
    X = Yi2X(Y,3);
elseif mode == 3
    X = shiftdim(Y, 2);
else
    X = Y;
end

x = X(:);

end
 