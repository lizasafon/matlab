function Q = ElectroStaticBalls(XYZ,R,F)
%R - ball's linear size
% XYZ - massiv of radius-vectors
% F - vector of ball's potentials

XYZ = XYZ';% for example with quadrupole 

N = length(R);
C = zeros(N, N);% potential coefficients
%checking if balls do not overlap
for ii = 1:(N - 1)
    for jj = (ii + 1):N
        if norm(XYZ(:,ii)-XYZ(:,jj)) <= R(ii) + R(jj)
            error('Too close')
        end
    end
end

for ii = 1:N
    for jj = ii:N
        if jj == ii
            C(ii, ii)= 1/R(ii);
        else
            C(ii, jj) = 1/norm(XYZ(:, ii) - XYZ(:,jj));
            C(jj, ii) = C(ii, jj);
        end
    end
end

Q = C \ F; %solving CQ=F system
end
