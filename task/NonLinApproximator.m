function [P, sgP] = NonLinApproximator (y,r,fun, P_0)
% P_0 - first linear approximation
N = size(y, 2);
delta = 1e^-10;

p_length = length(P_0);

for iterations = 1 : 1000

C = zeros(1, N);
for ii = 1 : N
    C(ii) = fun(r(:, ii), P_0); 
end

% find jacobian of fun
J = jacobian(fun, r, p_length, P_0);

% G*(P - P_0) = J*(y - fun(r, P_0))
Y = J'*(y' - C');
G = J'*J;
P = P_0' + G\Y;

if isnan(P)
    break
end

if iterations > 1
    if abs(P_0 - P_last) < delta
        break
    end
end

P_last = P;
P_0 = P';
end

end


function [J] = jacobian(fun, r, p_length, P_0)
step = 1e-6; 
N = length(r(1, :));

J = zeros(N, p_length);

for ii = 1 : N
    for jj = 1 : p_length
        delta = zeros(1, p_length);
            delta(jj) = step*abs(P_0(jj));
            J(ii, jj) = (fun(r(:, ii),P_0 + delta) - fun(r(:, ii),P_0 - delta))/(2*step);
    end   
end
end
