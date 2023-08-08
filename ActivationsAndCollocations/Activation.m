%-------------------------------------------------------------------------%
function [act, act_d1, act_d2, act_d3, act_d4] = Activation(x, w, b, ...
                                                     activationType)
%-------------------------------------------------------------------------%

%{
Decription:
This function determines the function outputs along with its first
derivative, second derivative, and third derivative for an activation
function. So far we only are using the sigmoid function. In the future we
will add more activation functions.

Author(s):
Kristofer Drozd
Enrico Schiassi

Date:
May, 09, 2019

Inputs:
x - independent variable
w - input weight
b - bias
    
Outputs:
act - function output
act_d1 - first derivative output
act_d2 - second derivative output
act_d3 - third derivative output
act_d4 - fourth derivative output
%}
    
N = length(x);
L = length(w);

x = repmat(x,1,L);
w = repmat(w',N,1);
b = repmat(b',N,1);

if  strcmp(activationType, 'logistic')

    act = 1./(exp(- b - w.*x) + 1);
    act_d1 = (w.*exp(- b - w.*x))./(exp(- b - w.*x) + 1).^2;
    act_d2 = (2.*w.^2.*exp(- 2.*b - 2.*w.*x))./(exp(- b - w.*x) + 1).^3 - (w.^2.*exp(- b - w.*x))./(exp(- b - w.*x) + 1).^2;
    act_d3 = (w.^3.*exp(- b - w.*x))./(exp(- b - w.*x) + 1).^2 - (6.*w.^3.*exp(- 2.*b - 2.*w.*x))./(exp(- b - w.*x) + 1).^3 + (6.*w.^3.*exp(- b - w.*x).*exp(- 2.*b - 2.*w.*x))./(exp(- b - w.*x) + 1).^4;
    act_d4 = (14.*w.^4.*exp(- 2.*b - 2.*w.*x))./(exp(- b - w.*x) + 1).^3 - (w.^4.*exp(- b - w.*x))./(exp(- b - w.*x) + 1).^2 + (24.*w.^4.*exp(- 4.*b - 4.*w.*x))./(exp(- b - w.*x) + 1).^5 - (36.*w.^4.*exp(- b - w.*x).*exp(- 2.*b - 2.*w.*x))./(exp(- b - w.*x) + 1).^4;

elseif  strcmp(activationType, 'tanh')

    act = (exp(b + w.*x) - exp(- b - w.*x))./(exp(b + w.*x) + exp(- b - w.*x));
    act_d1 =(w.*exp(b + w.*x) + w.*exp(- b - w.*x))./(exp(b + w.*x) + exp(- b - w.*x)) - ((exp(b + w.*x) - exp(- b - w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2;
    act_d2 = (2.*(exp(b + w.*x) - exp(- b - w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^2)./(exp(b + w.*x) + exp(- b - w.*x)).^3 - (2.*(w.*exp(b + w.*x) + w.*exp(- b - w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 - (w.^2.*exp(- b - w.*x) - w.^2.*exp(b + w.*x))./(exp(b + w.*x) + exp(- b - w.*x)) - ((exp(b + w.*x) - exp(- b - w.*x)).*(w.^2.*exp(- b - w.*x) + w.^2.*exp(b + w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2;
    act_d3 = (w.^3.*exp(- b - w.*x) + w.^3.*exp(b + w.*x))./(exp(b + w.*x) + exp(- b - w.*x)) - (6.*(exp(b + w.*x) - exp(- b - w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^3)./(exp(b + w.*x) + exp(- b - w.*x)).^4 + (6.*(w.*exp(b + w.*x) + w.*exp(- b - w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^2)./(exp(b + w.*x) + exp(- b - w.*x)).^3 + ((exp(b + w.*x) - exp(- b - w.*x)).*(w.^3.*exp(- b - w.*x) - w.^3.*exp(b + w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 - (3.*(w.^2.*exp(- b - w.*x) + w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) + w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 + (3.*(w.^2.*exp(- b - w.*x) - w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 + (6.*(exp(b + w.*x) - exp(- b - w.*x)).*(w.^2.*exp(- b - w.*x) + w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^3 ;
    act_d4 = ((w.^2.*exp(- b - w.*x) + w.^2.*exp(b + w.*x)).^2.*(6.*exp(b + w.*x) - 6.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^3 - (2.*(3.*w.^2.*exp(- b - w.*x) - 3.*w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^2)./(exp(b + w.*x) + exp(- b - w.*x)).^3 - ((6.*w.^2.*exp(- b - w.*x) - 6.*w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^2)./(exp(b + w.*x) + exp(- b - w.*x)).^3 - (w.^4.*exp(- b - w.*x) - w.^4.*exp(b + w.*x))./(exp(b + w.*x) + exp(- b - w.*x)) - (4.*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^3.*(6.*w.*exp(b + w.*x) + 6.*w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^4 + (4.*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^4.*(6.*exp(b + w.*x) - 6.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^5 - ((exp(b + w.*x) - exp(- b - w.*x)).*(w.^4.*exp(- b - w.*x) + w.^4.*exp(b + w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 + ((w.^2.*exp(- b - w.*x) + w.^2.*exp(b + w.*x)).*(3.*w.^2.*exp(- b - w.*x) - 3.*w.^2.*exp(b + w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 + ((w.^2.*exp(- b - w.*x) - w.^2.*exp(b + w.*x)).*(3.*w.^2.*exp(- b - w.*x) + 3.*w.^2.*exp(b + w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 - ((w.^3.*exp(- b - w.*x) + w.^3.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 + ((w.^3.*exp(- b - w.*x) - w.^3.*exp(b + w.*x)).*(w.*exp(b + w.*x) + w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 + ((3.*w.^3.*exp(- b - w.*x) - 3.*w.^3.*exp(b + w.*x)).*(w.*exp(b + w.*x) + w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 - ((3.*w.^3.*exp(- b - w.*x) + 3.*w.^3.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^2 + (2.*(3.*w.^2.*exp(- b - w.*x) + 3.*w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) + w.*exp(- b - w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^3 + (3.*(w.^2.*exp(- b - w.*x) + w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).*(6.*w.*exp(b + w.*x) + 6.*w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^3 - ((w.^3.*exp(- b - w.*x) - w.^3.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).*(6.*exp(b + w.*x) - 6.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^3 - (6.*(w.^2.*exp(- b - w.*x) + w.^2.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)).^2.*(6.*exp(b + w.*x) - 6.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^4 - (2.*(exp(b + w.*x) - exp(- b - w.*x)).*(w.^3.*exp(- b - w.*x) - w.^3.*exp(b + w.*x)).*(w.*exp(b + w.*x) - w.*exp(- b - w.*x)))./(exp(b + w.*x) + exp(- b - w.*x)).^3;

elseif  strcmp(activationType, 'sin')    

    act = sin(b + w.*x);
    act_d1 = w.*cos(b + w.*x);
    act_d2 = -w.^2.*sin(b + w.*x);
    act_d3 = -w.^3.*cos(b + w.*x);
    act_d4 = w.^4.*sin(b + w.*x);

elseif  strcmp(activationType, 'cos')

    act = cos(b + w.*x);
    act_d1 = -w.*sin(b + w.*x);
    act_d2 = -w.^2.*cos(b + w.*x);
    act_d3 = w.^3.*sin(b + w.*x) ;
    act_d4 = w.^4.*cos(b + w.*x) ;

elseif  strcmp(activationType, 'gaussian')

    act = exp(-(b + w.*x).^2);
    act_d1 = -2.*w.*exp(-(b + w.*x).^2).*(b + w.*x);
    act_d2 = 4.*w.^2.*exp(-(b + w.*x).^2).*(b + w.*x).^2 - 2.*w.^2.*exp(-(b + w.*x).^2);
    act_d3 = 12.*w.^3.*exp(-(b + w.*x).^2).*(b + w.*x) - 8.*w.^3.*exp(-(b + w.*x).^2).*(b + w.*x).^3;
    act_d4 = 4.*w.^4.*exp(-(b + w.*x).^2).*(4.*(b + w.*x).^4 - 12.*(b + w.*x).^2 + 3);

elseif  strcmp(activationType, 'arctan')

    act = atan(b + w.*x);
    act_d1 = w./((b + w.*x).^2 + 1);
    act_d2 = -(2.*w.^2.*(b + w.*x))./((b + w.*x).^2 + 1).^2;
    act_d3 = (8.*w.^3.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^3 - (2.*w.^3)./((b + w.*x).^2 + 1).^2 ;
    act_d4 = 24.*w.^4.*(b+w.*x).*(1-(b+w.*x).^2)./((b+w.*x).^2+1).^4;

elseif  strcmp(activationType, 'sinh')

    act = sinh(b + w.*x);
    act_d1 = w.*cosh(b + w.*x);
    act_d2 = w.^2.*sinh(b + w.*x);
    act_d3 = w.^3.*cosh(b + w.*x);
    act_d4 = w.^4.*sinh(b + w.*x);

elseif  strcmp(activationType, 'soft_plus')

    act = log(exp(b + w.*x) + 1);
    act_d1 = (w.*exp(b + w.*x))./(exp(b + w.*x) + 1);
    act_d2 = (w.^2.*exp(b + w.*x))./(exp(b + w.*x) + 1) - (w.^2.*exp(2.*b + 2.*w.*x))./(exp(b + w.*x) + 1).^2;
    act_d3 =  (w.^3.*exp(b + w.*x))./(exp(b + w.*x) + 1) - (3.*w.^3.*exp(2.*b + 2.*w.*x))./(exp(b + w.*x) + 1).^2 + (2.*w.^3.*exp(b + w.*x).*exp(2.*b + 2.*w.*x))./(exp(b + w.*x) + 1).^3 ;
    act_d4 =  w.^4.*exp(b+w.*x).*(-4.*exp(b+w.*x) + exp(2.*(b+w.*x)) + 1)./(exp(b+w.*x)+1).^4;

elseif  strcmp(activationType, 'bent')    

    act = b + ((b + w.*x).^2 + 1).^(1./2)./2 + w.*x - 1./2;
    act_d1 = w + (w.*(b + w.*x))./(2.*((b + w.*x).^2 + 1).^(1./2));
    act_d2 = w.^2./(2.*((b + w.*x).^2 + 1).^(1./2)) - (w.^2.*(b + w.*x).^2)./(2.*((b + w.*x).^2 + 1).^(3./2));
    act_d3 = (3.*w.^3.*(b + w.*x).^3)./(2.*((b + w.*x).^2 + 1).^(5./2)) - (3.*w.^3.*(b + w.*x))./(2.*((b + w.*x).^2 + 1).^(3./2));
    act_d4 = 3.*w.^4.*(4.*b.^2 + 8.*b.*w.*x + 4.*w.^2.*x.^2 - 1)./(2.*(b+w.*x).^2+1).^(7./2);

elseif  strcmp(activationType, 'inv_sinh')

    act = log(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x);
    act_d1 = (w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x);
    act_d2 = (w.^2./((b + w.*x).^2 + 1).^(1./2) - (w.^2.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^(3./2))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x) - (w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)).^2./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^2;
    act_d3 =   (2.*(w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)).^3)./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^3 - ((3.*w.^3.*(b + w.*x))./((b + w.*x).^2 + 1).^(3./2) - (3.*w.^3.*(b + w.*x).^3)./((b + w.*x).^2 + 1).^(5./2))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x) - (3.*(w.^2./((b + w.*x).^2 + 1).^(1./2) - (w.^2.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^(3./2)).*(w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^2;
    act_d4 =   (((3.*w.^3.*(b + w.*x))./((b + w.*x).^2 + 1).^(3./2) - (3.*w.^3.*(b + w.*x).^3)./((b + w.*x).^2 + 1).^(5./2)).*(w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^2 - ((3.*w.^4)./((b + w.*x).^2 + 1).^(3./2) - (18.*w.^4.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^(5./2) + (15.*w.^4.*(b + w.*x).^4)./((b + w.*x).^2 + 1).^(7./2))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x) - ((w.^2./((b + w.*x).^2 + 1).^(1./2) - (w.^2.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^(3./2)).*((3.*w.^2)./((b + w.*x).^2 + 1).^(1./2) - (3.*w.^2.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^(3./2)))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^2 - (6.*(w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)).^4)./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^4 + (((9.*w.^3.*(b + w.*x))./((b + w.*x).^2 + 1).^(3./2) - (9.*w.^3.*(b + w.*x).^3)./((b + w.*x).^2 + 1).^(5./2)).*(w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)))./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^2 + (6.*(w.^2./((b + w.*x).^2 + 1).^(1./2) - (w.^2.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^(3./2)).*(w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)).^2)./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^3 + (2.*((3.*w.^2)./((b + w.*x).^2 + 1).^(1./2) - (3.*w.^2.*(b + w.*x).^2)./((b + w.*x).^2 + 1).^(3./2)).*(w + (w.*(b + w.*x))./((b + w.*x).^2 + 1).^(1./2)).^2)./(b + ((b + w.*x).^2 + 1).^(1./2) + w.*x).^3;

elseif  strcmp(activationType, 'soft_sign')

    act = (b + w.*x)./(abs(b + w.*x) + 1);
    act_d1 = w./(abs(b + w.*x) + 1) - (w.*sign(b + w.*x).*(b + w.*x))./(abs(b + w.*x) + 1).^2;
    act_d2 =(2.*w.^2.*sign(b + w.*x).^2.*(b + w.*x))./(abs(b + w.*x) + 1).^3 - (2.*w.^2.*sign(b + w.*x))./(abs(b + w.*x) + 1).^2 - (2.*w.^2.*dirac(b + w.*x).*(b + w.*x))./(abs(b + w.*x) + 1).^2;
    act_d3 = (6.*w.^3.*sign(b + w.*x).^2)./(abs(b + w.*x) + 1).^3 - (6.*w.^3.*dirac(b + w.*x))./(abs(b + w.*x) + 1).^2 - (6.*w.^3.*sign(b + w.*x).^3.*(b + w.*x))./(abs(b + w.*x) + 1).^4 - (2.*w.^3.*(b + w.*x).*dirac(1, b + w.*x))./(abs(b + w.*x) + 1).^2 + (12.*w.^3.*dirac(b + w.*x).*sign(b + w.*x).*(b + w.*x))./(abs(b + w.*x) + 1).^3;
    act_d4 = (24.*w.^4.*dirac(b + w.*x).^2.*(b + w.*x))./(abs(b + w.*x) + 1).^3 - (8.*w.^4.*dirac(1, b + w.*x))./(abs(b + w.*x) + 1).^2 - (24.*w.^4.*sign(b + w.*x).^3)./(abs(b + w.*x) + 1).^4 + (24.*w.^4.*sign(b + w.*x).^4.*(b + w.*x))./(abs(b + w.*x) + 1).^5 + (48.*w.^4.*dirac(b + w.*x).*sign(b + w.*x))./(abs(b + w.*x) + 1).^3 - (2.*w.^4.*(b + w.*x).*dirac(2, b + w.*x))./(abs(b + w.*x) + 1).^2 + (16.*w.^4.*sign(b + w.*x).*(b + w.*x).*dirac(1, b + w.*x))./(abs(b + w.*x) + 1).^3 - (72.*w.^4.*dirac(b + w.*x).*sign(b + w.*x).^2.*(b + w.*x))./(abs(b + w.*x) + 1).^4;

elseif  strcmp(activationType, 'elu')
    
    if x <= 0
        act = 0.1*(exp(b + w.*x) - 1);
        act_d1 = 0.1*w.*exp(b + w.*x);
        act_d2 = 0;
        act_d3 = 0;
        act_d4 = 0;
    else
        act = x;
        act_d1 = 1;
        act_d2 = 0;
        act_d3 = 0;
        act_d4 = 0;
    end
    
    
end

end