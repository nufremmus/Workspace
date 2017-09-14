function distance = d_normal(x,y,k,ias)
x=full(x);
y=full(y);

x = reshape(x,[1,length(x)]);
y = reshape(y,[1,length(y)]);
ias = reshape(ias,[length(ias),1]);

if ~(isvector(x)&isvector(y))
    error('Input must be a vector')
end
%for jovanas data
% dd = (diag(x'*y))';
% distance = nthroot((((x-dd)*ias)^k+((y-dd)*ias)^k),k);
% distance = distance/((x+y-dd)*ias);

%for yuxiang's data;x and y are logic vectors
dd = x & y;
distance = nthroot((((x-dd)*ias)^k+((y-dd)*ias)^k),k);

%avoid NaNs
if distance==0
	return;
end

distance = distance/((x+y-dd)*ias);
distance = full(distance);

end
