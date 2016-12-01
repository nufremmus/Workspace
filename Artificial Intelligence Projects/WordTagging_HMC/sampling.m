function [new_tags Itag]= sampling(new_words)
%The idea here is to generate the tags in order with the knowlege of previous selected tags following the distribution p(state(i)|w)

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

p_temp = ones(1,12)/12; %substitute p_ws(j,id(t)) by pp_temp(j) when j-th word is new

%-------initialization--------------
new_tags = cell(M,5);
%Itag(i,:) stores the output tags' indices in state for i-th sample
Itag = zeros(M,5);

id = zeros(M,1);% id(i) stores the index of word new_words{i} in vocabulary

if isempty(new_words{1})
   fprintf('\n Err: the test sentence is empty\n');
   return;
end

%find the index of the words in the sentence new_words
id = GetIdxnew(new_words,vocabulary);

for sample = 1:5 %generate the sample 5 times

	%-------------find the output tag for the 1st word---------------
	%dist stores the distribution p(state(i)|w1)
	dist = zeros(1,12);
	for i = 1:12
	
	    if id(1)~= 0 % if the first word is in the vocabulary
			dist(i) = p_ws(i,id(1))*p1(i)*10000;
		else %if the first word is new
			dist = p1*10000;%take the distribution of the tags at the start of a sentence
		end		
	end  
	dist = dist/sum(dist); %normalize dist


	r = mnrnd(1,dist); %mnrnd(12,dist) is a multinomial sampling with n variables and distribution 'dist'
	Itag(1,sample) = find(r>0);
	
	%------recursively find the optimal tags for the rest words--------- 
	%the optimal tag for the (i-1)-th word plays a role in choosing that for i-th word
	 for j = 2:M
		dist = zeros(1,12);
		%fprintf(1,'Itag(%d) is %d id is %d\n',j,Itag(j-1,sample),id(j));
		if id(j) ~= 0 %if the j-th word is in the vocabulary
			for i = 1:12				
					dist(i) = p_ws(i,id(j))*trans(Itag(j-1,sample),i)*10000;
			end 
		else % if the j-th word is new, then use the transitional probability distribution from previous chosen tag to other tags
			for i = 1:12				
					dist(i) = p_temp(i)*trans(Itag(j-1,sample),i)*10000;
			end 
		end
		
		%normalize dist 
		if sum(dist) == 0 % if the entries are 0 in dist after rounding
		   dist =p_temp;
		 else
			dist = dist/sum(dist);
		end
		
		r = mnrnd(1,dist); %mnrnd(12,dist) is a multinomial sampling with n variables and distribution 'dist'
	    Itag(j,sample) = find(r>0);
	 end
		
	 for i = 1:M
		new_tags{i,sample} = state{Itag(i,sample)};
	 end
	
end
	
	