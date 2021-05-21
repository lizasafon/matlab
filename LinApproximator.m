
function [P,sgP] = LinApproximator(y,r,funcs)
 
N = length(y);
M  = size(funcs, 1);

%g - matrix with f_i(r_j) scalars
g = zeros(N, M);
for ii = 1 : N
    for jj = 1 : M
        f = cell2mat(funcs(jj));
        vec = num2cell(r(:, ii));
        g(ii, jj) =  f(vec{:});
    end
end


Y = (y*g)';
G = g'*g;
P = G\Y; % G*a = Y (from seminar, P=a)

%difference function
E_rr = norm(y - g*a)^2/N;

%from seninar
Gkk = trace(G);
sgP = sqrt(E_rr/(2*Gkk)); %sgP = Delta_a
end
