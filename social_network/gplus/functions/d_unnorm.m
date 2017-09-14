function distance = d_unnorm(x,y,k,ias)

x = reshape(x,[1,length(x)]);
y = reshape(y,[1,length(y)]);
ias = reshape(ias,[length(ias),1]);

if ~(isvector(x)&isvector(y))
    error('Input must be a vector')
end
dd = x & y;
distance = nthroot((((x-dd)*ias)^k+((y-dd)*ias)^k),k);
distance = full(distance);
end
