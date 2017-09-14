p=genpath(pwd);
addpath(p);

load qualified;

for Id = 1:length(qualified)
%------------------------------load data-------------------------------------
	fileID= qualified(Id);
	fprintf(1,'\n-----------------------');
	fprintf(1,'\n working on file %s\n',char(fileID));
	load(sprintf('matlab_data/feat/%s_feat',char(fileID)),'feat');
	load(sprintf('matlab_data/ias/%s_ias',char(fileID)),'ias');
	load(sprintf('matlab_data/circles/%s_circles',char(fileID)),'circles');
	load(sprintf('matlab_data/disjoint/%s_disjoint',char(fileID)),'disjoint');

	total_len = 0;
	for i=1:length(disjoint)
		i1 = disjoint(i);%index for a circle
		total_len = total_len+length(i1);
	end
	%now compute the distances within the circles
	distance= zeros(total_len);
	current_len = 0;
	prev_len =0;
	
	for i=1:length(disjoint)
		i1 = disjoint(i);%index for a circle
		circle=circles{i1}; %indices for the members in a circleircl
		current_len = current_len+length(circle);
		distance((prev_len+1):current_len,(prev_len+1):current_len) = 1;
		prev_len = current_len;
	end

	clear disjoint feat circles

end%for Id
