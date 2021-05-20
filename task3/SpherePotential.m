function [F,X,Y,P] = SpherePotential(XYZ,Q,R,r0,a,b,Dx,Dy,Nxy)
XYZ = XYZ';

%a - local x, b - local y

%forming user's basis
global e1 e2 e3;
e1 = a/norm(a);
e2 = b - dot(e1, b)*e1;
e3 = cross(e1, e2)/(norm(cross(e1, e2)));
S = [e1, e2, e3]; %transformation martix



%forming user's plane
xy = S\(XYZ - r0); % S*xy + r0 = XYZ

X = zeros(Nxy(1), Nxy(2));
Y = zeros(Nxy(1), Nxy(2));

for ii=1 : Nxy(1) %y-coordinates of points on user's plane
    Y(ii, :) = Dy(1): (Dy(2) - Dy(1))/(Nxy(1) - 1) : Dy(2); 
end
for jj = 1 : Nxy(2) % x-coordinates of points on user's plane
X(:, jj) = Dx(1) : (Dx(2) - Dx(1))/(Nxy(2) - 1) : Dx(2); 
end




F = zeros(Nxy(1), Nxy(2)); %potentials of points in plane


for ii = 1 : Nxy(1)
    for jj = 1 : Nxy(2) %running through all points
        for qq = 1 : size(Q, 1) %running through all charges
            if norm([X(ii, jj); Y(ii, jj); 0] - xy(:, qq)) < R(qq)
                F(ii, jj) = F(ii, jj) + Q(qq)/R(qq); %if point is in ball
            else
                F(ii, jj) = F(ii, jj) + Q(qq)/norm([X(ii, jj); Y(ii, jj); 0] - xy(:, qq));               
            end% phi = Q/r
        end
    end
end

figure; hold on; grid on; mesh(X,Y,F);
end

