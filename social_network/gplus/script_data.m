addpath(genpath(pwd));

load data% stores data; data is the cell that stores the file names
s=data;

for l=1:length(s)

	FileID = s(l);
	fprintf(1,'\n working on file %s',char(FileID));
	% fid = fopen(sprintf('gplus_data/%s.feat',char(FileID)),'r');
	% tline = fgetl(fid );
	% Nleng = length(tline);
	% feat = [];
	% while ischar(tline)	
		% feat = vertcat(feat,str2num(tline));
		% %disp(tline)
		% %fprintf(1,'%d\n',length(str2num(tline)));
		% tline = fgetl(fid);
	% end

	% save(sprintf('%s_feat',char(FileID)),'feat');
	% fclose(fid);
	
	% [m,n] = size(feat); %m is the number of people;n is the number of features

	%--------get the adjacency matrix ----------
	% fid = fopen(sprintf('gplus_data/%d.edges',FileID),'r');

	% formatSpec = '%d';
	% edges_raw = fscanf(fid,formatSpec);
	% fclose(fid)

	% edges(:,1) =(edges_raw(1:2:length(edges_raw)))'; 

	% edges(:,2) =(edges_raw(2:2:length(edges_raw)))';


	% index = features(:,1);

	% A = zeros(m,m);
	% N = length(edges);
	% for i = 1:N
		% j = edges(i,1);
		% k = edges(i,2);
		% if (j>0) && (k>0)
			% jj = find(index==j);
			% kk = find(index==k);
			% if (~isempty(jj)) && (~isempty(kk))
				% A(jj,kk) = 1;
				% A(kk,jj) = 1;
			% end
		% end
	% end

	% save(sprintf('%d_edges',FileID),'A');
	%---------------read the circles data------------
	fid = fopen(sprintf('gplus_data/%s.circles',char(FileID)),'r');
	tline = fgetl(fid);

	%###############count the lines in circles#######
	count = 0;
	while ischar(tline)
		count = count+1;
		tline = fgetl(fid);
	end

	%set-up
	circles = cell(count,1); i=1;

	frewind(fid); %go to the start of file

	%expression = 'circle\d+'; %regular expression that starts with 'circle' followed by digits
	tline = fgetl(fid);
	tline = tline(12:end);
	while ischar(tline)
		%tline = regexprep(tline,expression,''); 
		circles{i,1} = str2num(tline);
		tline = fgetl(fid);
		tline = tline(12:end);
		i=i+1;
	end
    save(sprintf('%s_circles',char(FileID)),'circles');
	
end %end for l

fclose(fid);
