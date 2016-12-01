%--------------declare global variables---------------------
global state p_ws p_s p1 trans vocabulary 

p = genpath(pwd);
addpath(p);

%load FinalData % stores : tag and word
load State  % stores: state
load P_ws % stores: p_ws
load P_s  % stores: p_s
load P_Initial % stores: p1 = p(S1 = state(1))
load P_Trans % stores: trans, which is transitional probability matrix
load Words %stores:vocabulary
% load Letters % stores: alpha = {'a';'b'; ...;'z'}
% load Stage  %stores: stage

method = cell(1,7);
method{1}= 'Naive';
method{2}= 'Bayes';
method{3}= 'Sample 1';
method{4}= 'Sample 2';
method{5}= 'Sample 3';
method{6}= 'Sample 4';
method{7}= 'Sample 5';
% method_name = cell(1,4);
% method_name{1}= 'Naive';
% method_name{2}= 'Bayes';
% method_name{3}= 'Bayes';
% method_name{4}= 'Sampling';

method_name = cell(1,3);
method_name{1}= 'Naive';
method_name{2}= 'Bayes';
method_name{3}= 'Sampling';

%-------------load testing data---------------
 % load the variables tag_test and word_test
%load('bc.tiny.test.mat');
%load('bc.small.test.mat');
load('bc.test.mat');

%-------------Initializations-----------------
[m n] = size(word_test);
M = max(m,n);

display = 1; %display = 1:show the output_tags

total_word=0;
total_sentence=0;

right_word=zeros(7,1); 

right_sentence = zeros(7,1);
correct_sentence = zeros(3,1);

accuracy_word=zeros(3,1);
accuracy_sentence = zeros(3,1);

start = 1;
[test_sentence,i_start,i_end] = GetSentence(start,word_test); % search for the first sentence

while start <=M
     %-----------reset right_sentence for each loop-----------
	 right_sentence = zeros(7,1);%stores the number of correct sentences for 8 output sentences in each run
	
    
	[test_sentence,i_start,i_end] = GetSentence(start,word_test);% search for a sentence
	 

	 if isempty(test_sentence{1})
	   % fprintf(1,'We have tested all sentences\n');
	   %------------------performance summary----------------------
	   	fprintf(1,'%15s','----------------');
		fprintf(1,'\nPERFORMANCE SUMMARY\n');
		accuracy_word(1:2) = right_word(1:2)/total_word;
		accuracy_word(3) = sum(right_word(3:7))/(5*total_word);
		
		accuracy_sentence =  correct_sentence/total_sentence;
		
		%-----------------display the performance summary ------------------
		fprintf(1,'Total words:%d total sentences :%d \n',total_word,total_sentence);
		for l=1:3	
			fprintf(1,'%20s: \n',method_name{l});
			fprintf(1,'%20s: %f \n','Words accuracy',accuracy_word(l));
			fprintf(1,'%20s: %f \n','Sentence accuracy',accuracy_sentence(l));			
		end
		
		break; % we have reached the end of the testing data
	 end
	 
	 if ~isempty(test_sentence{1})
	 	M_s = max(size(test_sentence));%get the length of the test sentence
	   total_sentence = total_sentence+1;
		total_word = total_word+M_s;
	 end
	 
	 %output_tag stores the output indices of the tags for each method
	 %output_tag(:,1) : indices for Naive Bayes
	 %output_tag(:,2) : indices for Bayes
	 %output_tag(:,3) : indices for Bayes2
	 %output_tag(:,4:8) : indices for sampling
	
	output_tag = zeros(M_s,7);
		      s = cell(M_s,7);

	 %-----get the optimal indices for Naive Bayes------
	 for j=1:M_s
		[ss I]=SimpleBayes(test_sentence{j});%s is a string here
		s{j,1}=ss;
		output_tag(j,1) =  I;	       
	 end 
	 
	 
	 
	 %-----get the optimal indices for Bayes------
	 [ss I] = Bayes(test_sentence);
	 s(:,2) =ss;
	 output_tag(:,2) = I; 

	 %-----get the optimal indices for sampling--------
	 [ss I]= sampling(test_sentence);	 
	 output_tag(:,3:7) = I;
	 
	  s(:,3:7)  = ss; 
	 
	 %-----------------compute the accuracy-related terms----------------------
	 true_tag = tag_test(i_start:i_end);
	 output_true = GetIdxnew(true_tag,state); %tag_test store the true tags in the testing data; output_true are the indices for the true tags
	 
	 for k=1:7	 % compare each row of output_tag to output_true
	     I_cmp = find(output_tag(:,k)==output_true);
		 
		 n_right = max(size(I_cmp));
		 right_word(k) = right_word(k)+n_right;
		 if n_right == M_s% get the whole sentence right 
			right_sentence(k) = 1;
		 end
	  end
	%-----update the correct number of sentences for each method
	 correct_sentence(1:2) = correct_sentence(1:2)+right_sentence(1:2);
	 
	 tf = sum(right_sentence(3:7))>0 ;%tf = 1 if there is at least one correct sentence from sampling
	 correct_sentence(3) =  correct_sentence(3)+tf;
		
	 %----------------------display the output tags---------------------------------------
	 if display
			fprintf(1,'%15s','----------------');
			fprintf(1,'\n%15s:','Test sentence:');%show the method name
			SentenceDisp(test_sentence);
			
			fprintf(1,'%15s:','Ground truth');%show the method name
			SentenceDisp(true_tag);
			
		for k=1:7
			fprintf(1,'%15s:',method{k});%show the method name
			SentenceDisp(s(:,k));
		end
	 end
	 
	 %----------------------updata i_start and get a new sentence--------------------------------------
	 start = i_end+1;	 
	 %-------------------save data ------------------------	 
	 if mod(total_sentence,1000)==0
		save results total_sentence total_word accuracy_sentence accuracy_word correct_sentence right_word
	 end
	 
end 

save results total_sentence total_word accuracy_sentence accuracy_word correct_sentence right_word
	 
	 