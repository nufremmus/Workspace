% to find initial probability from training data
%--add paths to the directories
p = genpath(pwd);
addpath(p);

load FinalData;
load state;

% count(i) stores the number of sentences that start from state i.
count = zeros(12,1);

%'state' stores the tags in order
% state = cell(12,1);
% state=['ADJ '; 'ADV ' ;'ADP ' ;'CONJ'; 'DET ' ;'NOUN'; 'NUM '; 'PRON'; 'PRT ' ;'VERB' ;'X   '; '.   ' ];
% state = cellstr(state);


[m n] = size(tag);
M = max(m,n);% M is the length of 'tag'

previous = tag(1);
i = look4(previous,state);
count(i) = count(i)+1; %the first tag among all tags is an initial state

for l=1:M
	current = tag(l);

	j = look4(current,state); % Find the number of column where the tag is in

	%the previous is space and the current is non-space means the current is an initial state; i==0 means that the corresponding element in the tag is space
	if i==0 & j~=0
	   count(j) = count(j) + 1;
	end
	i=j; %for the next round, the index for current state becomes that for the previous state

	%save the variables count and l every 10000 steps
	if mod(l,10000)==0
	     l
		save Count_initial count l
	end
end

save Count_initial count l

totalcount = sum(count); %i-th entry stores the sum of the numbers of transitions from a specific state i to different state j's;

p_initial = count/totalcount;

save InitialProbability p_initial
