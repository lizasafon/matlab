function E = KronigPenney(k) 


%initializing constants
global h k_n m a b U0 Emax e0;
h = 1.054571817 * 10^(-34);
m = 9.109383701528 * 10^(-31);
a = 0.5e-9;
b = 2e-9;
U0 = -4 * e0;
Emax = 26 * e0;
e0 = 1.602176634e-19;

k = k / (a+b);
E = zeros(size(k));


%cutting points, where F(x) > 3
energy_range = [];
for ii = U0/e0: (Emax - U0)/(1000*e0): Emax/e0
%     F(ii*e0) - cos(k_n * (a+b))
    if abs(F(ii*e0)) <= 3
        energy_range = [energy_range, 0];
        energy_range(length(energy_range)) = ii;
    end
end

%producing some plots if needed

% f = zeros(length(energy_range));
% for ii=1:length(energy_range)
%     f(ii) = delta(energy_range(ii));
% end
% 
% figure; hold on; grid on; plot(energy_range, f)


figure; grid on; hold on; xlabel('k*(a+b)'); ylabel('Energy, eV')
for n = 1 : length(k) 
    k_n = k(n);
    column = []; %column for solutions of F(E) = cos(k(a+b)) for the given k
    %cycle that finds all solutions 
    for ii = 1 : length(energy_range)-1
        if delta(energy_range(ii))*delta(energy_range(ii+1)) < 0
         solution = fzero(@delta, [energy_range(ii) energy_range(ii+1)]);
         if isempty(column(abs(column - solution) < 1e-10)) 
                column = [column; 0];
                column(length(column)) = energy_range(ii); %adds solution to column
         end   
         end
    end
    %plots the results upon the previous plots
    plot(k_n * (a+b), column, 'linestyle', 'none', 'marker', '.')
    E(1:length(column), n) = column; %inserts the column into E matrix
end

end


%left function
function x = F(E)
global h m a b U0;
    mu = sqrt(2 * m * E) / h;
    lambda = sqrt(2 * m * (E - U0)) / h;
    x = cos(mu*a) * cos(lambda*b) - (lambda^2 + mu^2)/(2*mu*lambda)* sin(mu*a) * sin(lambda*b);
end

%right function
function x = delta(E)
global k_n a b e0;
E = E *e0;
    x = F(E) - cos(k_n * (a+b));
end
