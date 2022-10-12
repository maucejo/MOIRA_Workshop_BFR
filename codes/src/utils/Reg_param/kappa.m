function [k, n, e] = kappa(X, f)
%     [k, n, e] = KAPPA(X, f) Returns the discrete curvature (1D), normals 
%     and evolute of a curve 'X' in the Euclidean space. The parameter 'f'
%     determines wether rotating frame corrections are performed.
    % Parse input arguments
    [N, p] = size(X);
    if N < p
        X = X';
        [N, p] = size(X);
    end
    if nargin < 2
        f = 0;
    end
    % Obtain pairs of combinations between adjacent triplets
    X = [X, zeros(N, 3 - p)];
    A = X(1 : N - 2, :) - X(3 : N, :);
    B = X(2 : N - 1, :) - X(3 : N, :);
    X = X(:, 1 : p);
    % Pre-calculate cross-products
    I = cross(A, B, 2);
    J = cross(B, I, 2);
    K = cross(A, I, 2);
    % Obtain relative circumcentre
    A = A .* A;
    I = I .* I;
    O = (sum(A, 2) .* J - sum(B .* B, 2) .* K) ./ sum(I, 2) / 2;
    O = O(:, 1 : p);
    % Initialise extrapolating coefficients
    I = (X(2, :) - X(1, :)) ./ (X(3, :) - X(2, :));
    J = (X(N, :) - X(N - 1, :)) ./ (X(N - 1, :) - X(N - 2, :));
    
    % Calculate end-point extrapolated values
    i = I .* (O(2, :) - O(1, :));
    j = J .* (O(N - 2, :) - O(N - 3, :));
    C = [O(1, :) - i; O; O(N - 2, :) + j];
    if norm(X(N, :) - X(1, :)) < 1e-6
        C([1, N], :) = [1; 1] .* (C(1, :) + C(N, :)) / 2;
    end
    % Derive 1-D curvature from radii
    R = sqrt(sum(C .* C, 2));
    i = isnan(R);
    R(i) = Inf;
    k = 1 ./ R;
    if nargout < 2
        return
    end
    % Calculate normal vectors
    n = C .* k;
    if any(i)
        t = [X(2, :) - X(1, :); B(:, 1 : p)];
        t = [t; zeros(1, p)] + [zeros(1, p); t];
        if p == 1
            n = NaN * k;
        end
        if p == 2
            n = [-t(:, 2), t(:, 1)];
        end
        if p == 3
            n(i, :) = cross(t(i, :), rand(sum(i), 3), 2);
        end
        n = n ./ sqrt(sum(n .* n, 2));
    end
    
    % Frenet corrections via rotating frame
    if f
        [~, idx] = max(k);
        val = n(idx, 1);
        
        s = sum(n(1 : N - 1, :) .* n(2 : N, :) < 0, 2) == p;
        s = s .* (1 : N - 1)';
        s = s(s > 0);
        while any(s)
            if numel(s) == 1
                if s > N / 2
                    if s == N
                        s = s - 1;
                    end
                    n(s + 1 : N, :) = -n(s + 1 : N, :);
                else
                    n(1 : s, :) = -n(1 : s, :);
                end
                break
            end
            n(s(1) + 1 : s(2), :) = -n(s(1) + 1 : s(2), :);
            s = sum(n(1 : N - 1, :) .* n(2 : N, :) < 0, 2) == p;
            s = s .* (1 : N - 1)';
            s = s(s > 0);
        end
        if n(idx, 1) ~= val(1)
            n = -n;
        end
    end
    if nargout < 3
        return
    end
    % Calculate evolute with end-point extrapolation
    O = O + X(3 : N, :);
    I = I .* (O(2, :) - O(1, :));
    J = J .* (O(N - 2, :) - O(N - 3, :));
    e = [O(1, :) - I; O; O(N - 2, :) + J];
    if p > 1
        e(isnan(e)) = Inf;
    end
end