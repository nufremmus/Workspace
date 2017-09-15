function [Coeff,GG] = Coeff_yang(k,method)

%call Coeff_Matrix
G = Coeff_Matrix(k);

%Initialization
SizeInt =2*k-1;
Coeff = zeros(k,1);
b = ones(k,1);
GG  = zeros (k,SizeInt); %%GG is the augmented big matrix with rows from G

for i =1 : k
GG(i,i:i+k-1) = G(i,:);
end

%----------get coefficients-------------------
switch method
	case 'Penrose'
		%solve for GG^T*GG = \lambda*[1,1,1...1]^T
		GenInv_GG = pinv(GG*GG');
		Coeff = GenInv_GG*b;
        Coeff = Coeff/sum(Coeff);
	case 'Uniform'
		Coeff = b/sum(b);
end
