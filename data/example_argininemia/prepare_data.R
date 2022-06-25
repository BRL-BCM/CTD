
print(sprintf("Current directory is %s: ", getwd()))

# Load the Miller2015_Heparin dataset
data(Miller2015)
fil.rate=as.numeric(Miller2015$`Times identifed in all 200 samples`[-1])/200
names(fil.rate) = rownames(Miller2015)[-1]
data_mx = Miller2015[,grep("IEM_", colnames(Miller2015))]
data_mx = data_mx[which(fil.rate>0.90), ]
rmMets = names(which(apply(data_mx, 1, function(i) any(i>20))))
if (length(rmMets)>0) {
  data_mx = data_mx[-which(rownames(data_mx) %in% rmMets),]
}
#data_mx = as.matrix(data_mx[-c(which(rownames(data_mx) == "diagnosis")),grep("IEM_", #colnames(data_mx))])
arg_data = data_mx[,which(diags=="Argininemia")]    
# To experimental (disease) CSV 
write.table(arg_data, file = file.path('data', 'example_argininemia', 'experimental.csv'), row.names=TRUE, sep=',')
# To CSV control
ind = which(diags=="No biochemical genetic diagnosis")
write.csv(data_mx[,ind], file = file.path('data', 'example_argininemia', 'control.csv'), row.names=TRUE)
