function Anew = convertA(A)

[M N] = size(A);

lags = (M-1)/N;

ch = N;

Anew = permute(reshape(A(2:end,:),[ch lags ch]),[1 3 2]);

end