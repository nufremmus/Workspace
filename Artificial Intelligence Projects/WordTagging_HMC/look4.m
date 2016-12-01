function index = look4(x,y,k) %x is a string; y is a 1-d cell array storing strings
% k =1 for any cell array y; k =2 for cell array Words

% s = 'a b c d e f g h i j k l m n o p q r s t u v w x y z';
% s = strsplit(s,' ');
%stage stores the starting index of each letter in the cell array Words
global alpha stage

[m n ] = size(y);
M = max(m,n);

index = 0;

start = 1;
ending = M;

stage = [stage;M] ; % add stage(27) = M

%-------------special treatments to strings containing digits----------
tf = isstrprop(x, 'digit'); %tf(i) = 1 iff  x(i) is a digit

if sum(tf) > 0  % there is at least one digit in the string x
	tf1 = isstrprop(x,'alpha');
	if sum(tf1) > 0 % there is at least one letter in the string x
		index = 2;
		return;
    else % w has digits but no letters
		index = 1;
		return;
	end
end
%-----------------------------------------------------------------------

%------------------special treatment to y = Words -----------------------
if k == 2 % y is the cell array 'words'
	if isletter(x(1))
	   ii = strmatch(lower(x(1)),alpha); % find the index of letter x(1) in size
	   start = stage(ii);
	   ending = stage(ii+1);
	 end
	 
	for i = start:ending
	  if strcmpi(x,y(i)) == 1 % if string x and string y(i) match; letter case doesn't matter
		index = i;
		break;
	  end
	 end
elseif k ==1  %for general cell array y
	for i = 1:M
	  if strcmpi(x,y(i)) == 1 % if string x and string y(i) match; letter case doesn't matter
		index = i;
		break;
	  end
	end
end