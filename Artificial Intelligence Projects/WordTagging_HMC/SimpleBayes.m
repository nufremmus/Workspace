function [new_tag answer] = SimpleBayes(new_word)

%new_word is the test word;it should be a string
%new_tag is the output tag corresponding to the new_word; it is a string.

%---------------load test data ----------------------

 % pp_ws stores the probability p(w,s) for the same w but different s; pp_ws(i) = p(w,state(i))
 
%declare global variables
%vocabulary : stores sorted training vocabulary; words is a cell array
%state : stores the 12 tags;state is a cell array
%p_s : stores the distribution of the tags in the training data; p_s(i) stores the percentage of state(i)
%p_ws : stores the conditional probability p(w|s); p_ws(i,j) = p(wj|si).
global state p_ws p_s vocabulary alpha stage
 
%initialization
pp_ws  = zeros(12,1);
 %tag(answer) is the optimal tag for new_word
answer =0;
new_tag='';

if isempty(new_word)
   fprintf('\n Err: the test word is empty\n');
   return;
end

i_w = look4(new_word,vocabulary,2); %find the index of new_word in Words
if i_w == 0 % new_word is not in the training vocabulary
   %fprintf(1,'\n "%s" is not in the vocabulary.\n',new_word);
   
	pp_ws = p_s; %sample the tags from the distribution p_s
	[M,answer] = max(pp_ws); % find argmax_i p(w,s(i)) ; 
	new_tag = state{answer};
	return;
end;

	
for j = 1:12  % calculate p(w,s) for different tags s
	s = state(j);
	pp_ws(j) = p_ws(j,i_w)*p_s(j);
end
			
	[M,answer] = max(pp_ws); % find argmax_i p(w,s(i)) ; 
	new_tag = state{answer}; %new_tag is the optimal tag
	
	% fprintf(1,'\n%5s\n',new_word); % display the test word and the optimal tag
	% fprintf(1,'%5s',new_tag); 
	


	
	
	
	
	
	