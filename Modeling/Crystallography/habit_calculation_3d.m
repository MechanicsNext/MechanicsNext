function [cb0, cm0, varargout] = habit_calculation_3d(U_matrix_cell_array)
%% M-M and A-M interface calculation
% Input type: E = Transformation strain is specified for variants
%             U = Stretch tensor is specified for variants
U = U_matrix_cell_array;
%
number_of_cv = size(U, 1);
%
twin_num = 2*nchoosek(size(U, 1), 2);
twin_elements = [];
am_elements = [];
corresponding_twin = [];
% Stretch tensors (Assuming strains are Greene-Lagrange)
for nt=1:number_of_cv
    for mt=1:nt-1
        if(nt ~= mt)
            U1 = U{nt};
            U2 = U{mt};
            %% Martensite-martensite interface
            % Using Bhattacharya 5.16
            C = transpose(inv(U2))*transpose(U1)*U1*inv(U2);
            [eigvec, eigval] = eigs(C);
            lambda1 = eigval(1,1);
            lambda2 = eigval(2,2);
            lambda3 = eigval(3,3);
            %
            e1 = eigvec(:, 1);
            e2 = eigvec(:, 2);
            e3 = eigvec(:, 3);
            %
            if(lambda1 <= 1.0 && lambda3 >= 1.0)
                a1 = sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1 + sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3;
                a2 = sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1 - sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3;
                n1 = ((sqrt(lambda3) - sqrt(lambda1))/sqrt(lambda3 - lambda1)*(-sqrt(1-lambda1)*U2'*e1 + sqrt(lambda3-1)*U2'*e3));
                n2 = ((sqrt(lambda3) - sqrt(lambda1))/sqrt(lambda3 - lambda1)*(-sqrt(1-lambda1)*U2'*e1 - sqrt(lambda3-1)*U2'*e3));
                %
                a1 = a1*norm(n1); % First solution
                n1 = n1/norm(n1);
                a2 = a2*norm(n2); % Second solution
                n2 = n2/norm(n2);
                % Calculate rotation
                Q1 = (U2 + kron(a1, n1')) * inv(U1);
                Q2 = (U2 + kron(a2, n2')) * inv(U1);
                %
                twin_elements(end+1, :) = [nt mt a1' n1' reshape(Q1, 1, 9)];
                twin_elements(end+1, :) = [nt mt a2' n2' reshape(Q2, 1, 9)];
            else
                disp('No solution to twinning equation')
            end
        end
    end
end
%% Martensite-austenite interface
debug_info = [];
for mt=1:size(twin_elements, 1)
    U1 = U{twin_elements(mt, 2)};
    %
    a = twin_elements(mt, 3:5);
    n = twin_elements(mt, 6:8);
    %
    delta = a*((U1*inv(U1^2 - eye(3)))*n');
    eta = trace(U1^2) - det(U1^2) - 2 + (norm(a)^2)/(2*delta);
    debug_info(end+1, :) = [twin_elements(mt, 1) twin_elements(mt, 2) delta eta 0 0 0 0];
    %
    if(delta < -2 && eta >= 0)
        % First pair of solutions
        lambda = 0.5*(1 - sqrt(1 + 2.0/delta));
        % Make sure 0 <= lambda <= 1
        if(lambda >= 0 && lambda <= 1)
            debug_info(end, 5) = lambda;
            C = (U1 + lambda*kron(n, a'))*(U1 + lambda*kron(a, n'));
            [eigvec, eigval] = eigs(C);
            % Make sure that the middle eigenvalue of C equals 1
            if(abs(eigval(2,2) - 1.0) < 0.01)
                lambda1 = eigval(1,1);
                lambda2 = eigval(2,2);
                lambda3 = eigval(3,3);
                %
                debug_info(end, 6:8) = [eigval(1,1) eigval(2,2) eigval(3,3)];
                %
                e1 = eigvec(:, 1);
                e2 = eigvec(:, 2);
                e3 = eigvec(:, 3);
                %
                b1 = sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1 + sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3;
                b2 = sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1 - sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3;
                m1 = (sqrt(lambda3) - sqrt(lambda1))/sqrt(lambda3-lambda1)*(-sqrt(1-lambda1)*e1 + sqrt(lambda3-1)*e3);
                m2 = (sqrt(lambda3) - sqrt(lambda1))/sqrt(lambda3-lambda1)*(-sqrt(1-lambda1)*e1 - sqrt(lambda3-1)*e3);
                %
                b1 = b1*norm(m1);
                m1 = m1/norm(m1);
                b2 = b2*norm(m2);
                m2 = m2/norm(m2);
                %
                Q1 = (eye(3, 3) + kron(b1, m1')) * inv( (eye(3, 3) + lambda * kron(a, (inv(U1)*n' ))) * U1 );
                Q2 = (eye(3, 3) + kron(b2, m2')) * inv( (eye(3, 3) + lambda * kron(a, (inv(U1)*n' ))) * U1 );
                % Save b, m
                am_elements(end+1, :) = [twin_elements(mt, 1) twin_elements(mt, 2) b1' m1' lambda reshape(Q1, 1, 9)];
                am_elements(end+1, :) = [twin_elements(mt, 1) twin_elements(mt, 2) b2' m2' lambda reshape(Q2, 1, 9)];
                corresponding_twin(end+1) = mt;
                corresponding_twin(end+1) = mt;
            end
        end
        % Second pair of solutions
        lambda = 1.0 - 0.5*(1 - sqrt(1 + 2.0/delta));
        % Make sure 0 <= lambda <= 1
        if(lambda >= 0 && lambda <= 1)
            C = (U1 + lambda*kron(n, a'))*(U1 + lambda*kron(a, n'));
            [eigvec, eigval] = eigs(C);
            % Make sure that the middle eigenvalue of C equals 1
            if(abs(eigval(2,2) - 1.0) < 0.01)
                lambda1 = eigval(1,1);
                % lambda2 = eigval(2,2);
                lambda3 = eigval(3,3);
                %
                e1 = eigvec(:, 3);
                % e2 = eigvec(:, 2);
                e3 = eigvec(:, 1);
                %
                b3 = sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1 + sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3;
                b4 = sqrt(lambda3*(1-lambda1)/(lambda3-lambda1))*e1 - sqrt(lambda1*(lambda3-1)/(lambda3-lambda1))*e3;
                m3 = (sqrt(lambda3) - sqrt(lambda1))/sqrt(lambda3-lambda1)*(-sqrt(1-lambda1)*e1 + sqrt(lambda3-1)*e3);
                m4 = (sqrt(lambda3) - sqrt(lambda1))/sqrt(lambda3-lambda1)*(-sqrt(1-lambda1)*e1 - sqrt(lambda3-1)*e3);
                %
                b3 = b3*norm(m3);
                m3 = m3/norm(m3);
                b4 = b4*norm(m4);
                m4 = m4/norm(m4);
                %
                Q3 = (eye(3, 3) + kron(b3, m3')) * inv( (eye(3, 3) + lambda * kron(a, (inv(U1)*n' ))) * U1 );
                Q4 = (eye(3, 3) + kron(b4, m4')) * inv( (eye(3, 3) + lambda * kron(a, (inv(U1)*n' ))) * U1 );
                % Save b, m
                am_elements(end+1, :) = [twin_elements(mt, 1) twin_elements(mt, 2) b3' m3' lambda reshape(Q3, 1, 9)];
                am_elements(end+1, :) = [twin_elements(mt, 1) twin_elements(mt, 2) b4' m4' lambda reshape(Q4, 1, 9)];
                corresponding_twin(end+1) = mt;
                corresponding_twin(end+1) = mt;
            end
        end
    end
end

cb0 = am_elements(:, 3:5)';
cm0 = am_elements(:, 6:8)';

if(nargout == 3)
    varargout{1} = twin_elements;
elseif(nargout == 4)
    varargout{1} = twin_elements;
    varargout{2} = am_elements;
elseif(nargout == 5)
    varargout{1} = twin_elements;
    varargout{2} = am_elements;    
    varargout{3} = corresponding_twin;    
end
end