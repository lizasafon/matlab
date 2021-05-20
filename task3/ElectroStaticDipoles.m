function [Q,D] = ElectroStaticDipoles(XYZ,R,F)
XYZ = XYZ';
N = length(R);
F = [F; zeros(3 * N, 1)];
C = zeros(4*N, 4*N);% potential coefficients

%checking if balls do not overlap
for y_global_index = 1:(N - 1)
    for x_global_index = (y_global_index + 1):N
        if norm(XYZ(:,y_global_index)-XYZ(:,x_global_index)) <= R(y_global_index) + R(x_global_index)
            error('Too close')
        end
    end
    
%producing potential matrix (high blocks)
for y_global_index = 1 : N
    for x_global_index = 1 : 4 * N
        block_local_index = mod(x_global_index, N);
        if(block_local_index == 0)
            block_local_index = N; 
        end %for not to do fourh cycle
        blocknumber_x = (x_global_index - block_local_index) / N;
        delta_r = XYZ(1:3,y_global_index)-XYZ(1:3,block_local_index);
        if(y_global_index == x_global_index)
            C(y_global_index, x_global_index) = 1 / R(y_global_index);
        else
            if(blocknumber_x == 0)
                C(y_global_index, x_global_index) = 1 / norm(delta_r);
            else
                if(block_local_index == y_global_index)
                    C(y_global_index, x_global_index) = 0;
                else
                    C(y_global_index, x_global_index) = delta_r(blocknumber_x) / norm(delta_r) ^ 3;
                end
            end
        end  
    end
end

%producing other blocks of matrix
for y_global_index = (N + 1) : (4 * N)
    for x_global_index = 1 : (4 * N)
        block_local_index_x = mod(x_global_index, N);
        block_local_index_y = mod(y_global_index, N);
        if(block_local_index_x == 0)
            block_local_index_x = N;
        end
        if(block_local_index_y == 0)
            block_local_index_y = N;
        end
        block_number_x = (x_global_index - block_local_index_x) / N;
        block_number_y = (y_global_index - block_local_index_y) / N;
        delta_r = XYZ(1:3,block_local_index_y)-XYZ(1:3,block_local_index_x);
        if(block_local_index_y == x_global_index)
            M(y_global_index, x_global_index) = 0;
        else
            if(y_global_index == x_global_index)
                M(y_global_index, x_global_index) = 1 / R(block_local_index_y) ^ 3;
            else
                if(block_number_x == 0)
                    M(y_global_index, x_global_index) = delta_r(block_number_y) / norm(delta_r) ^ 3;
                else
                    if(block_local_index_x == block_local_index_y)
                        M(y_global_index, x_global_index) = 0;
                    else
                        if(block_number_y == block_number_x)
                            M(y_global_index, x_global_index) = (3 * delta_r(block_number_y) ^ 2 - norm(delta_r)^2 ) / norm(delta_r) ^ 5;
                        else
                            M(y_global_index, x_global_index) = 3 * delta_r(block_number_y) * delta_r(block_number_x)/ norm(delta_r) ^ 5;
                        end
                    end
                end
            end
        end 
    end
end
Final_vector = M \ F;
Q = Final_vector(1:N, 1);
Dx = Final_vector(N+1:2*N,1)';
Dy = Final_vector(2*N+1:3*N,1)';
Dz = Final_vector(3*N+1:4*N,1)';
D = [Dx; Dy; Dz]';
D = -D;



end