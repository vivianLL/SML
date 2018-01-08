function b = gauss(x,mu,sigma)
d = length(mu);
c = (2*pi)^(d/2);
c = c*(det(sigma)^(1/2));
b= exp((-1/2)*((x-mu)'*inv(sigma)*(x-mu)));
b = b/c;
return
end