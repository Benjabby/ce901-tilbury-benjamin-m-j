function [predImage, predT, predA] = dehazeFattal(img)

NS = 256 ;
N = NS * NS ;

A = [0.8 0.8 0.9]' ;

R = [0.8 0.4 0.8]' ;

Rt = R - A * (R' * A) / (norm(A)^2);

true_eta = A' * R / (norm(A) * norm(Rt)) ;

%
% Generating synthetic transmission & shading functions
%

i=0;
for y=1:NS
    for x=1:NS
    i=i+1 ;
    
    xx = x / NS - 0.5 ;
    yy = y / NS - 0.5 ;

    rr = norm([xx yy]) ;
    
    t(i) = 0.3 + 0.2 * sin(rr * 60) + 0.1 * rand ;
    l(i) = 1 + 0.3 * sin((xx+yy) * 80) + 0.1 * rand ;
    
    end
end

t=t';
l=l';
    

%
% Generating synthetic input image
%

i=0;
for y=1:NS
    for x=1:NS
    i=i+1 ;

    I(i,:) = t(i) .* l(i) .* R + (1-t(i)) .* A ;
    
    im(x,y,1) = I(i,1) ; im(x,y,2) = I(i,2) ; im(x,y,3) = I(i,3) ;
    tim(x,y,1) = l(i) .* R(1) ; tim(x,y,2) = l(i) .* R(2) ; tim(x,y,3) = l(i) .* R(3) ; 
    end
end

% estimating eta, trans. (and shading)

disp(size(I))

[est_t est_l est_eta] = estimate(A,I) ;

assignin("base", "fuckmatlab", I);

est_t = reshape(est_t,NS,NS) ;
est_l = reshape(est_l,NS,NS) ;


t = reshape(t,NS,NS) ;
l = reshape(l,NS,NS) ;


%
% Generating output haze-free image
%

i=0;
for y=1:NS
    for x=1:NS
    i=i+1 ;

    J(i,:) = (I(i,:) - (1-est_t(i)) * A') / est_t(i) ;
    
    oim(x,y,1) = J(i,1) ; oim(x,y,2) = J(i,2) ; oim(x,y,3) = J(i,3) ;
    end
end

figure(1)
clf

ai = ones(3,3,3) ;
ai(2,2,:) = A ;
ai = uint8(ai*255) ; 
subplot(2,4,5), imagesc(ai), title(sprintf('True eta = %f \nEst. eta = %f\nTrue A=%.2f,%.2f,%.2f\n\nAirlight color',true_eta,est_eta,A(1),A(2),A(3))) ;

colormap(gray(256))
subplot(2,4,3), imagesc(est_t), title('Est. trans.') ;
subplot(2,4,4), imagesc(est_l), title('Est. shading') ;
subplot(2,4,7), imagesc(t), title('True trans.') ;
subplot(2,4,8), imagesc(l), title('True shading') ;

im = uint8(im*255) ;
subplot(2,4,1), imshow(im), title('Input image') ;
oim = uint8(oim*255) ;
subplot(2,4,2), imshow(oim), title('Output image') ;
tim = uint8(tim*255) ;
subplot(2,4,6), imshow(tim), title('True image') ;
%
% Here is the decorrelation procedure that computes eta, and then trans.
% (t) and the shading (l)
%

function [est_t est_l est_eta] = estimate(A, I)

IA = I * A / norm(A) ;
IR = sqrt(sum(I.^2,2) - IA.^2);
IR = real(IR);
H = (norm(A) - IA) ./ IR ;
H(isinf(H)) = max(H(~isinf(H)));

C = cov(IA, H) ./ cov(IR, H) ;

est_eta = C(1,2) ;

est_t = 1 - (IA - est_eta * IR) / norm(A) ;
est_l = IR ./ est_t ;