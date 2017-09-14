function dy_int = derivative_interior(a,b,N, data, k,method)

% grid : the points on real line
% data : the polluted function value corresponding to grid points
% k : the number of points on the substencils 
G = Coeff_Matrix(k);

[C,GG] = Coeff_Emma(k,method);

%initialization
dy_int = zeros(N+1,1);

y_stencil = zeros(2*k-1,1);
	
for i=k:N-k+2
    y_stencil = data(i-k+1:i+k-1,1);
    dy_int(i) = C'*GG*y_stencil;
end

dy_int = dy_int/((b-a)/N);

