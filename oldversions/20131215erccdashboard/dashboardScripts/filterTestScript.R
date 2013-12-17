filterIdx = NULL
for (i in 1:length(Lab5ERCConly$Feature)){
 idxAdd  = which(Lab5ERCConly[i,-c(1)] == 0)
 if (length(idxAdd) >= 2){
   filterIdx = append(filterIdx, i) 
 }
 
}
print(filterIdx)