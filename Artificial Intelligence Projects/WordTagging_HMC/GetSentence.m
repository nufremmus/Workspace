function [y,i_start,i_end]=GetSentence(start,x)
%x is a cell array with strings
%i is the starting index for the sentence in x
%we are trying to get a sentence starting at x(start)
%y returns the cell array
%i_start and i_end store the starting and ending positions in x

[m,n] = size(x);
M = max(m,n);

%----initialization
y=cell(1);
i_end= start;

while isempty(x{start}) & start < M
	start= start+1;
end

i_start = start;
i_end= i_start;

if i_start==M
	fprintf(1,'We have reached the end of the array.\n');
	return;
end

i_end= start;

while ~isempty(x{i_end}) & i_end < M
	i_end = i_end+1;
end

i_end = i_end-1;

if i_end>i_start
y=x(start:i_end);
end

if i_end==i_start
y=x(start);
end