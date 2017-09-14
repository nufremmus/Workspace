qualified <- scan("F:/R_data_gplus/qualified.txt",what="character", sep="\t")

dir[1] = "unnorm_ias_txt"
dir[2] = "unnorm_simple_txt"
p_unnorm_wc = matrix(data = NA,nrow = length(qualified),ncol = 2)
p_unnorm_t = matrix(data = NA,nrow = length(qualified),ncol = 2)

for(i in 1:length(qualified)){
  fileID = qualified[i]
  for(j in 1:2){
    #read in d_within and d_between
    sprintf("i=%d, j= %d",i,j)
    d_within <- scan(sprintf("F:/R_data_gplus/%s/%s_within.txt",dir[j],fileID))
    d_between <- scan(sprintf("F:/R_data_gplus/%s/%s_between.txt",dir[j],fileID))
    
    #carry out one-sided test on if the mean(d_within) < mean(d_between)
    p_unnorm_wc[i,j] = wilcox.test(d_between,d_within,alternative = "greater")$p.value
    p_unnorm_t[i,j] = t.test(d_between,d_within,alternative = "greater")$p.value
  }
}