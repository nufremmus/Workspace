Please refer to “Project_Description.pdf” and “report.pdf” for detailed problem statement and solutions.

The following is a brief description of the files in this project.


P_s: stores the distribution of the tags including those not appearing at the start of a sentence. P_s(i) stores the percentage of state(i)

P_initial: stores the distribution of the tags only appearing at the start of a sentence

P_ws

P_trans

count_s

State

Words



new_words1: it's based on new_words; the strings involving numbers are
combined;'111' represent the strings with numbers but no letters; '111AAA'
represent those that have both

pp_ws stores the probability p(w,s) for the same w but different s; pp_ws(i) =
p(w,state(i))



state : stores the 12 tags; state is a cell array

p_ws : stores the conditional probability p(w|s); p_ws(i,j) = p(wj|si).

vocabulary (Training Data, processed from “bc.train”): stores sorted training vocabulary; words is a cell array

