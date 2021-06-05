function [est_t est_l est_eta] = estimate(A, I)

IA = I * A / norm(A) ;
IR = sqrt(sum(I.^2,2) - IA.^2);
%IR = real(IR);
H = (norm(A) - IA) ./ IR ;
H(isinf(H)) = max(H(~isinf(H)));

C = cov(IA, H) ./ cov(IR, H) ;

est_eta = C(1,2) ;

est_t = 1 - (IA - est_eta * IR) / norm(A) ;
est_t = real(est_t);
est_l = IR ./ est_t ;
est_l = real(est_l);