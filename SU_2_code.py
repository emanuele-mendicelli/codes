#Importing for use dwave system
from collections import defaultdict
from dwave.system import DWaveSampler, EmbeddingComposite

#Import for python
import numpy as np
import sys


#Function that generates the tables associated with the use on k numbers of qubits per coefficient
def qubit_table(k):
    TQ=np.zeros((k, k))

    for i in range(1, k+1):
        for j in range(1, k+1):
            if (i != k and j != k):
                TQ[i-1][j-1] = (pow(2, i + j - 2 * k))
            elif i != k and j == k:
                TQ[i-1][j-1] = (-(pow(2, i - k)))
            elif i == k and j != k:
                TQ[i-1][j-1] = (-(pow(2, j - k)))
            else:
                TQ[i-1][j-1] = 1
    return TQ



#Reads the values of the variables from the line on the termina after program_name.py
k = int(sys.argv[1])
y = float(sys.argv[2])
chain = float(sys.argv[3])
lam = float(sys.argv[4])
n_reads = int(sys.argv[5])



#Input SU(2) matrix 1 plaquette
M=np.zeros((3, 3))

M[0][0]=-lam
M[0][1]=4*y
M[0][2]=0

M[1][0]=4*y
M[1][1]=3-lam
M[1][2]=y

M[2][0]=0
M[2][1]=y
M[2][2]=3-lam
#---------------------------------------------


P=np.kron(M,qubit_table(k))

Q = defaultdict(float)
Q=np.triu(P, k=1) + np.triu(P, k=0)




sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample_qubo(Q,num_reads=n_reads,chain_strength=chain)

print("k=",k,"y=",y,"chain=",chain,"lam=",lam,"n_reads=",n_reads)
print(sampleset)


#Print the results in a file
out_file ="RES_k_{}_y_{}_ch_{}_l_{}_nr_{}.txt".format(k,y,chain,lam,n_reads)
with open(out_file,'a') as file:
    sys.stdout = file # Change the standard output to the file we created.
    print("\n RUN """)
    print("k=",k,"y=",y,"chain=",chain,"lam=",lam,"n_reads=",n_reads)
    print(sampleset)

