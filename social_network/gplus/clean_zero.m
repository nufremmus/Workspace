p=genpath(pwd);
addpath(p);

load qualified

difference = zeros(2,2,length(qualified)); %(i_choice, p, Id)



 for Id = 1:length(qualified)
 % for Id = 3
	%------------------------------load data-------------------------------------);
		list=[];
		fileID= char(qualified{Id});
		fprintf(1,'\n-------------------------------');
		fprintf(1,'\nworking on file           %s\n',fileID);
		load(sprintf('matlab_data/feat/%s_feat.mat',fileID),'feat');
		
		[mf,nf] = size(feat);
		%check if there is zero vectors in features
		i=1;
		while(i<=mf)
			if(sum(feat(i,2:end))==0)
				list = [list;feat(i,1)];
				feat(i,:)=[];
				
				i=i-1;
				mf = mf-1;
			end%ifcl
			i=i+1;
		end%while

		load(sprintf('matlab_data/circles/%s_circles.mat',fileID),'circles');
		%clear out the members that have zero feature vectors
		for i=1:length(circles)
			ToDel = intersect(circles{i},list);
			if ~isempty(ToDel)
				fprintf(1,'Something is there\n');
			end
			circles{i}(ismember(ToDel,circles{i}))=[];
		end
		
		save(sprintf('matlab_data/new_feat/%s_feat',fileID),'feat');
		save(sprintf('matlab_data/new_circles/%s_circles',fileID),'circles');
		

end %Id

