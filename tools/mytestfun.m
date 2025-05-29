function f = mytestfun(x,y,z)
global num_evals
%f = tanh(5*(x+y+z));
%f = cosh(5*(x+y+z)).^(-2);

f = 1./(1 + 25*(x.^2 + y.^2 + z.^2));

%f = tanh(5*(x+z)) + exp(y);
%f = cos(10*(x+y+z));
num_evals = num_evals + numel(f);
end 