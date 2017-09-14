addpath(genpath(pwd));
load qualified;

for i =1:length(qualified)
	fileID = qualified{i};
	load(sprintf('matlab_data/new_feat/%s_feat',fileID),'feat')
	fprintf(1,'\n working on file %s',fileID);
	features = feat;
	features = full(features);
	[m n] = size(features);
	ias = zeros(n-1,1); %the first column is indices
	
	%total = sum(sum(features()));

	for i=1:n-1
		s = sum(features(:,i+1)); %features's first columns store user ids
		ias(i) = log((m+1)/(s+1));
	end

	save(sprintf('matlab_data/new_ias/%s_ias',fileID),'ias');
	clear features m n ias fileID
end