addpath(genpath(pwd));

load qualified;
%change directory
directory = 'normal_ias';


for i=1:length(qualified);
	FileString=qualified{i};
	InputFile = sprintf('%s/%s_%s_1.mat',directory,'distances',FileString);
	load(InputFile,'d_between','d_within');
	
	FileName_between = sprintf('%s_txt/%s_between.txt',directory,FileString);
	FileID_between = fopen(FileName_between,'w');
	fprintf(FileID_between,'%6.6E\n',d_between);
	fclose(FileID_between);
	

	FileName_within = sprintf('%s_txt/%s_within.txt',directory,FileString);
	FileID_within =fopen(FileName_within,'w');
	fprintf(FileID_within,'%6.6E\n',d_within);
	fclose(FileID_within);
end