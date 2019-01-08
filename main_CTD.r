require(CTD)

# 1. DATA LOAD ---------------------------------------------------------------
# Load your patient data (p features as rows x n observations as columns)
# data_mx=read.table("/your/own/data.txt", sep="\t", header=TRUE)
data(testData)
data_mx=t(testData)
rownames(data_mx)=tolower(rownames(data_mx))


# IMPORT OR LEARN BACKGROUND KNOWLEDGE GRAPH ---------------------------------------------------------
# Read in any network via its adjacency matrix
tmp=matrix(1, nrow=100, ncol=100)
for (i in 1:100) {
  for (j in 1:100) {
    tmp[i, j]=rnorm(1, mean=0, sd=1)
  }
}
colnames(tmp)=sprintf("MolPheno%d", 1:100)
ig=graph.adjacency(tmp, mode="undirected", weighted=TRUE, add.colnames="name")
V(ig)$name=tolower(V(ig)$name)


# SET GLOBAL VARIABLES, TUNING PARAMETERS ---------------------------------
# This is the percentage of probability that gets divided uniformly between the variables in the background knowledge graph ("ig").
p0=0.1
# This is the percentage of probability that gets divided preferentially (depending on the connectivity) between the variables in
# the background knowledge graph.
p1=0.9
# This is the threshold of probability for which the probability diffusion algorithm stops recursively diffusing.
thresholdDiff=0.01
# This is the size of the largest subset you will encode against your background knowledge graph.
kmx=5
# Global placeholder for the probabilities that are assigned to each variable in the background graph.
G=vector(mode="list", length=length(V(ig)$name))
names(G)=V(ig)$name
# Global copy of the background knowledge graph's adjacency matrix.
adjacency_matrix=list(as.matrix(get.adjacency(ig, attr="weight")))

# NOTE: It is strongly advised that you save your progress up until this point into a smartly named RData file,
#       so you can begin at stage 3 at any given time.
path= paste(getwd(), "/projectFolder", sep="")
save.image(sprintf("%s/experimentPrep.RData", path))






# GET NODE PERMUTATIONS FOR THE GRAPH -------------------------------------
# Note: If you are quantifying probabilities for thousands of subsets at a time,
#       learning all permutations up front is advised. If you only want to quantify
#       only a few subsets, you only need to process the permutations starting with
#       the variables in the subsets of interest. The below code processes all possible
#       starting nodes and propagates diffusion until all nodes have been drawn.

# Get node permutations, starting with all possible nodes in background knowledge graph
perms=list()
for (n in 1:length(G)) {
 print(sprintf("Generating node permutation starting with node %s", names(G)[n]))
 perms[[n]]=mle.getPermN(n, G)
}
names(perms)=names(G)



# GET BITSTRINGS ASSOCIATED WITH PATIENTS TOP PERTURBED VARIABLES ---------
# Get bitstrings associated with each patient's top kmx variable subsets
ptBSbyK=list()
for (pt in 1:ncol(data_mx)) {
 ptID=colnames(data_mx)[pt]
 ptBSbyK[[ptID]]=mle.getPtBSbyK(data_mx, ptID, perms, kmx)
}



# GET ENCODING LENGTHS OF PATIENT BITSTRINGS ------------------------------
# Identify the most significant subset per patient, given the background graph
data_mx.pvals=t(apply(data_mx, c(1,2), function(i) 2*pnorm(abs(i), lower.tail=FALSE)))
for (pt in 1:ncol(data_mx)) {
   ptID=colnames(data_mx)[pt]
   res=mle.getEncodingLength(ptBSbyK[[ptID]], data_mx.pvals, ptID, G)
   res=res[which.max(res[,"d.score"]),]
   print(res)
}



# GET PATIENT SIMILARITY SCORES USING MUTUAL INFORMATION ------------------
# Get patient distances
data_mx.pvals=apply(data_mx, c(1,2), function(i) 2*pnorm(abs(i), lower.tail=FALSE))
res=list()
t=list(ncd=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)),
        dir=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)),
        jac=matrix(NA, nrow=ncol(data_mx), ncol=ncol(data_mx)))
rownames(t$ncd)=colnames(data_mx)
colnames(t$ncd)=colnames(data_mx)
rownames(t$dir)=colnames(data_mx)
colnames(t$dir)=colnames(data_mx)
rownames(t$jac)=colnames(data_mx)
colnames(t$jac)=colnames(data_mx)
for (i in 1:(kmx-1)) {
 res[[i]]=t
}
for (pt in 1:ncol(data_mx)) {
 print(pt)
 ptID=colnames(data_mx)[pt]
 for (pt2 in pt:ncol(data_mx)) {
   ptID2=colnames(data_mx)[pt2]
   for (k in 1:(kmx-1)) {
     tmp=mle.getPatientSimilarity(ptBSbyK[[ptID]][k], ptID, ptBSbyK[[ptID2]][k], ptID2, data_mx, perms)
     res[[k]]$ncd[ptID, ptID2]=tmp$NCD
     res[[k]]$dir[ptID, ptID2]=tmp$dirSim
     res[[k]]$ncd[ptID2, ptID]=tmp$NCD
     res[[k]]$dir[ptID2, ptID]=tmp$dirSim

     p1.sig.nodes=rownames(data_mx)[order(abs(data_mx[,ptID]), decreasing=TRUE)][1:k]
     p2.sig.nodes=rownames(data_mx)[order(abs(data_mx[,ptID2]), decreasing=TRUE)][1:k]
     p1.dirs=data_mx[p1.sig.nodes, ptID]
     p1.dirs[which(!(p1.dirs>0))]=0
     p1.dirs[which(p1.dirs>0)]=1
     p2.dirs=data_mx[p2.sig.nodes, ptID2]
     p2.dirs[which(!(p2.dirs>0))]=0
     p2.dirs[which(p2.dirs>0)]=1
     p1.sig.nodes=sprintf("%s%d", p1.sig.nodes, p1.dirs)
     p2.sig.nodes=sprintf("%s%d", p2.sig.nodes, p2.dirs)
     res[[k]]$jac[ptID, ptID2]=1-(length(intersect(p1.sig.nodes, p2.sig.nodes))/length(union(p1.sig.nodes, p2.sig.nodes)))
     res[[k]]$jac[ptID2, ptID]=1-(length(intersect(p1.sig.nodes, p2.sig.nodes))/length(union(p1.sig.nodes, p2.sig.nodes)))
   }
 }
}



# VISUALIZATIONS ----------------------------------------------------------

# Multi-dimensional scaling
# if you have diagnostic labels associated with the colnames(data_mx), send them using diagnoses parameter
diagnoses=colnames(data_mx)
diagnoses[1:25]="diseased"
diagnoses[26:50]="neg_control"
patientSim=0.8*res[[k]]$ncd+0.2*res[[k]]$dir
p=plot.mdsSim(patientSim, diagnoses, k=2, diag="diseased")
p
p=plot.mdsSim(patientSim, diagnoses, k=3, diag="diseased")
p


# Hierarchical clustering
plot.hmSim(patientSim, path=getwd(), diagnoses)

