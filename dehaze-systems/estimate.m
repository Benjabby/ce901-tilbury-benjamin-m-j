function est_t = estimate(A, I)

[m, n, ~] = size(I);
A = reshape(A, [1 3])';
I = reshape(I,[],3);
IA = I * A / norm(A) ;
IR = sqrt(sum(I.^2,2) - IA.^2);
%IR = real(IR);
H = (norm(A) - IA) ./ IR ;
H(isinf(H)) = 0;%max(H(~isinf(H)));

C = cov(IA, H) ./ cov(IR, H) ;

est_eta = C(1,2) ;

est_t = 1 - (IA - est_eta * IR) / norm(A) ;
% est_t = real(est_t);


est_t = reshape(est_t,m,n);

end