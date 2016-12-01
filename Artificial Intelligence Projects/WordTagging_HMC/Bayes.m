function [new_tags Imax] = Bayes(new_words)
%The idea here is to use Viterbi algorithm, which is essentially dynamic programming. It'll be described in details in the report.

%new_words is the test sentence;it should be a cell array of strings.
%new_tags is the output tags corresponding to the new_words; it is a cell array of strings.

%vocabulary stores sorted training vocabulary; words is a cell array
%state stores the 12 tags;state is a cell array
%p_s stores the distribution of the tags in the training data; p_s(i) stores the percentage of state(i)
%p_ws stores the conditional probability p(w|s); p_ws(i,j) = p(wj|si).
%p1 stores the distribution of tags at the start of a sentence
%trans(i,j) = p(S_{t} = state(j)|S_{t-1} = state(i))
global state p_ws p_s p1 trans vocabulary alpha stage
 
[m n] = size(new_words);
M = max(m,n); % M is the length of the test sentence

%initialization
f  = zeros(M,12); % f(i,j) stores the maximum probability for the case where the final state is Si = state(j) with the sentence(1:i). This will be explained more in details in the report.
new_tags = cell(M,1);
Imax = zeros(M,1);

p_temp = ones(1,12)/12; %substitute p_ws(j,id(t)) by pp_temp(j) when j-th word is new

id = zeros(M,1);% id(i) stores the index of word new_words{i} in vocabulary

if isempty(new_words{1})
   fprintf('\n Err: the test sentence is empty\n');
   return;
end

%find the index of the words in the sentence new_words
% the index for a new word is 0
id = GetIdxnew(new_words,vocabulary);

% for i=1:M
    % if id(i) == 0 % there is a new word
		% %fprintf(1,'Viterbi: Cannot handle new words \n');
		% return;
	% end
% end

%initialization of f
if id(1)~=0 % the 1st word is in vocabulary
	for i = 1:12
		f(1,i) = p1(i)*p_ws(i,id(1));
		% fprintf(1,'i = %d : %f\n',i, p_ws(i,id(1)));
	end
else
	for i = 1:12
		f(1,i) = p1(i)*p_temp(i);
		% fprintf(1,'i = %d : %f\n',i, p_ws(i,id(1)));
	end
end


    
for t = 2:M
	
	if id(t) ~=0 %t-th word is in the vocabulary
	
		for j = 1:12
			temp = zeros(12,1);	
			for i = 1:12
			   temp(i) = f(t-1,i)*trans(i,j)*p_ws(j,id(t));
			end
			  f(t,j) = max(temp);
		end
	else % t-th word is not in the vocabulary
		for j = 1:12
			temp = zeros(12,1);	
			for i = 1:12
			   temp(i) = f(t-1,i)*trans(i,j)*p_temp(j);
			end
			  f(t,j) = max(temp);
		end
	end
	
	
end

% recursively recover the tags that give rise to the max probability of p(w|s)


	% fprintf(1,'\n%5s\n',new_word); % display the test word
	% fprintf(1,'%5s',new_tag); % display the optimal tag
	
[MM Imax(M)] =max(f(M,:)) ;

for t = M-1:-1:1
	temp = zeros(12,1);
	for i = 1:12
		temp(i) = f(t,i)*trans(i,Imax(t+1));
    end
	[MM Imax(t)] = max(temp);
end

for i = 1:M
	new_tags{i,1} = state{Imax(i)};
end
	
	
	
	