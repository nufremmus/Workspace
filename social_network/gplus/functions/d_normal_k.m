function distance = d_normal_k(x,y,k)
 
% C=ones(9571,1);
% if ~(isvector(x)&isvector(y))
    % error('Input must be a vector')
% end
% dd = (diag(x'*y))';

% distance = (x-dd).*(x-dd)+(y-dd).*(y-dd);
% %distance = nthroot(((x-dd)*C)^k+((y-dd)*C)^k,k);
% distance = distance/((x+y-dd).*(x+y-dd));
% toc
distance = 0;
%distance = bitxor(x,y) / bitor(x,y); % does not apply to sparse matrices
if (sum(abs(x-y))+x*y')~=0
distance = nthroot(sum(max(x-y,0))^k+sum(max(y-x,0))^k,k)/ (sum(abs(x-y))+x*y');
%fprintf(1,'%d\n',sum(sum(max(full(x-y),0))))
%fprintf(1,'%d\n',sum(sum(max(full(y-x),0))))

end

distance = full(distance);

end
