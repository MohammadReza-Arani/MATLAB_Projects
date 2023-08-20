function y_hat = lm_func(t,p,c)
% global example_number
example_number=c;
%  example 1:  easy for LM ... a poor initial guess is ok
if example_number == 1
   y_hat = p(1)*exp(-t/p(2)) + p(3)*t.*exp(-t/p(4));
end

%  example 2: medium for LM ... local minima
if example_number == 2
   mt = max(t);
   y_hat = p(1)*(t/mt) + p(2)*(t/mt).^2 + p(3)*(t/mt).^3 + p(4)*(t/mt).^4;
end

% example 3: difficult for LM ... needs a very good initial guess
if example_number == 3
   y_hat = p(1)*exp(-t/p(2)) + p(3)*sin(t/p(4));
end

