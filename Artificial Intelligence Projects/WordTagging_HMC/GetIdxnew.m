function id = GetIdxnew(x,y)
%find the indices of the strings that are in the cell array x in y
%e.g. y=state or y = vocabulary
[m n] = size(x);
M = max(m,n);%

id = zeros(M,1);
k=1;

%----------for words having digits---------
for i =1:M
   tf = isstrprop(x{i}, 'digit'); %tf(i) = 1 iff  x{i}(i) is a digit

	if sum(tf) > 0  % there is at least one digit in the string x
		tf1 = isstrprop(x{i},'alpha');
		if sum(tf1) > 0 % there is at least one letter in the string x
			index = 2;
			return;
		else % w has digits but no letters
			index = 1;
			return;
		end
	end

end

   [b id] = ismember(lower(x),y);
