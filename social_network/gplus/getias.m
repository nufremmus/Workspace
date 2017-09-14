addpath(genpath(pwd));
load data;% stores data; data is the cell that stores the file names
ss=data;


for i =1:length(ss)
	fileID = ss(i);
	load(sprintf('matlab_data/feat/%s_feat',char(fileID)),'feat')
	fprintf(1,'\n working on file %s',char(fileID));
	features = feat;
	features = full(features);
	[m n] = size(features);
	ias = zeros(n-1,1); %the first column is indices
	
	total = sum(sum(features));

	for i=1:n-1
		s = sum(features(:,i+1)); %features's first columns store user ids
		ias(i) = log((total+1)/(s+1));
	end

	save(sprintf('matlab_data/ias/%s_ias',char(fileID)),'ias');
	clear features m n ias fileID
end