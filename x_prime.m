function [x_output]=x_prime(a)
x=1;

x_prime=zeros(1,10);
for z=1:a
    x_prime(z)=x*z;
end

