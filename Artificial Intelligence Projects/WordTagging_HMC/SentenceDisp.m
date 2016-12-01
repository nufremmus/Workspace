function [] = SentenceDisp(y) %display a cell array with tags

[m n ] = size(y);
M = max(m,n);

% for i = 1:M
	% if isempty(y{i})
		% fprintf(1,'The sentence has a blank\n');
		% return;
	% end
% end

for i = 1:M
	fprintf(1,'%5s ',y{i});
end

	fprintf(1,'\n',y{i});