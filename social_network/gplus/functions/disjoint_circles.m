addpath(genpath(pwd));

load data;% stores data; data is the cell that stores the file names
s=data;
qualified = [];


for i=1:length(s)
	fileID = s(i);
	fprintf(1,'\n working on file %s',char(fileID));

	load(sprintf('%s_circles',char(fileID)),'circles')
	[m n] = size(circles);
	
	disjoint=[];
	union=[];
	
	for i =1:max(m,n)
		if (~isempty(circles)) && (length(circles{i})>9)
			if isempty(intersect(union,circles{i}))
				disjoint = vertcat(disjoint,i);
				union = horzcat(union,circles{i});
			end
		end
	end
	
	fprintf(1,'\n number of disjoint circles are  %d',length(disjoint));
	
	if length(disjoint)>1
		save(sprintf('%s_disjoint',char(fileID)),'disjoint');
		qualified=vertcat(qualified,fileID);
	end
    clear disjoint union circles
end

save qualified qualified