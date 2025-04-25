%filtro chebyshev
omegap = 200; 
omegas = 300; 
% ap = 2; 
% as = 8;
% 
% deltap = 10^(-ap/20); deltas = 10^(-as/20);

deltap=0.95;
deltas=0.15;

epsilon = sqrt((1/deltap^2) - 1)
n = ceil(acosh(sqrt((deltas^-2) - 1)/sqrt((deltap^-2) - 1))/acosh(omegas/omegap))

c0 = 0.5*((sqrt(1 + epsilon^-2) + epsilon^-1)^(1/n)(sqrt(1 + epsilon^-2) + epsilon^-1)^(-1/n))
num = 1;
den = 1;

if mod(n,2) == 0
    k = n/2;                    
    for i = 1:k
        bk(i) = 2*c0*sin((((2*i) - 1)*pi)/(2*n))
        ck(i) = c0^2 + (cos((((2*i) - 1)*pi)/(2*n)))^2
        num = conv(num,[0 0 (omegap^n)*ck(i)*1/sqrt(1 + epsilon^2)]);
        den = conv(den,[1 bk(i)*omegap ck(i)*omegap^2]);
    end
    
else
   
    k = (n-1)/2;
    for i = 1:k
       bk(i) = 2*c0*sin((((2*i) - 1)*pi)/(2*n))
       ck(i) = c0^2 + (cos((((2*i) - 1)*pi)/(2*n)))^2
       num = conv(num,[0 0 c0*(omegap^n)*ck(i)]);
       den = conv(den,[1 bk(i)*omegap ck(i)*omegap^2]);
    end
    den = conv(den,[1,omegap*c0]);
end

sys = tf(num,den)

bode(sys,[omegap,omegas]), grid on;

%Bode en valor abs
figure(2)
options = bodeoptions;
options.MagUnits = 'abs';
bode(sys,options), grid on;

