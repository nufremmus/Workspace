function distance = d_unnorm_k(x,y,k)
 
% C=ones(9571,1);
% if ~(isvector(x)&isvector(y))
    % error('Input must be a vector')
% end
% dd = (diag(x'*y))';

% distance = (x-dd).*(x-dd)+(y-dd).*(y-dd);
% %distance = nthroot(((x-dd)*C)^k+((y-dd)*C)^k,k);
% distance = distance/((x+y-dd).*(x+y-dd));
% toc
%distance = bitxor(x,y) / bitor(x,y); % does not apply to sparse matrices
distance = nthroot(sum(max(x-y,0).^k+max(y-x,0).^k),k);
distance = full(distance);

end
