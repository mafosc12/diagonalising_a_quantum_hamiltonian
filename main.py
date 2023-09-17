import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npl
import functools as ft
import scipy.sparse as sparse
import scipy.sparse.linalg as spsl

# Define our Pauli Matrices

I = np.array([[1, 0], [0, 1]])
sigmaX = np.array([[0, 1],[1, 0]])
sigmaZ = np.array([[1, 0], [0, -1]])

def CompleteHamiltonian(N, h=1.):
    """
    Constructs the Hamiltonian matrix of the transverse field Ising model for a system of N spins using periodic boundary conditions and in a transverse field h.
    Does this using Kronecker products of the Pauli matrices.
    Stores the Hamiltonian matrix in a numpy.ndarray object.
    
    Parameters:
    -----------
    N: integer
    h: float type, optional
    
    Returns:
    -----------
    H: ndarray
    
    """
    Matrices1 = np.zeros((N, 2, 2))
    Matrices1[:] = I
    Matrices2 = np.copy(Matrices1)
    Matrices1[0] = sigmaZ
    Matrices1[1] = sigmaZ

    Matrices2[0] = h * sigmaX
    
    
    H = np.zeros((2**N, 2**N))
    for ii in range(N):
        H -= ft.reduce(np.kron, Matrices1)
        Matrices1 = np.roll(Matrices1, 1, axis=0)
        
        H -= ft.reduce(np.kron, Matrices2)
        Matrices2 = np.roll(Matrices2, 1, axis=0)
        
    
    return H

def LocalHamiltonian(h):
    Hloc = -(np.kron(I, sigmaX)+np.kron(sigmaX,I)*(h/2)+np.kron(sigmaz, sigmaz))/r    # Your code goes here
    return Hloc

def CompleteHamiltonian2(N, h=1.):
  
    I = np.identity(N-k-1), dtype=float)
    H = np.zeros((2**N, 2**N))
    for k in range(N):
        H+=np.kron(identity(k), np.kron(LocalHamiltonian(h), (np.kron.I)*(h/2))
    return H  

# Check that the two methods do the same thing.
np.all(np.isclose(CompleteHamiltonian(10, h = 1.23), CompleteHamiltonian2(10, h = 1.23)))

hh = np.logspace(-2, 2, 10)
Ns = np.arange(2, 10, 1)

for N in Ns:
    energy = []
    for i in hh:
        E, Psi = npl.eigh(CompleteHamiltonian(N, h=h))
        energy.append(E[0])
    plt.plot(hh, energy)
    plt.yscale('symlog')
    plt.xscale('log')

Ns = np.arange(2, 11, 1)
h = 1

timingCompact = []
for N in Ns:
    print(f"N={N}")
    result = %timeit -o npl.eigh(CompleteHamiltonian2(N))
    timingCompact.append(result)
timeavg = [time.average for time in timingCompact]
timesd = [time.stdev for time in timingcompact]

# The averages can be put in a list by doing: [time.average for time in timingCompact]
# Same for the standard deviation: [time.stdev for time in timingCompact]

plt.errorbar(Ns, timeavg, yerr=timesd, fmt='.',capsize=5.)
plt.show()

# Define our Pauli Matrices

sparseI = sparse.csr_matrix(I)
sparseX = sparse.csr_matrix(sigmaX)
sparseZ = sparse.csr_matrix(sigmaZ)

# Use these to construct the Hamiltonian. 

def SparseHamiltonian(N, h=1):
    
    
    for k in range(N-1):
        H += sparse.kron(sparse.identity(2**k), sparse.kron(sparse.crs_matrix(LocalHamiltonian(h)), sparse.identity(2**(N-k-2)))
    
    return H

# Check that your Hamiltonian is really a sparse object:
type(SparseHamiltonian(10))

sparsesize = [SparseHamiltonian(N).nnz for N in np.arrange(2, 17, 1)}
completesize = [2**(2*N) for N in np.arrange(2, 17, 1)]
plt.plot(np.arange(2,17,1), sparsesize, '.', label = 'sparse hamiltonian')
plt.plot(np.arange*2,17,1), completesize, '+', label= 'complete hamiltonian')
plt.yscale('log')
plt.legend()
plt.show()

def ActWithH(psi, N, h=1.):
    
    
    psiout = np.zeros(psi.size)
    for k in range(N-1):
        psiout += np.tensordot(LocalHamiltonian(h).reshape(4,4),
        axes[[1], [1]]).transpose(1,0,2).reshape(2**N)
    psiout += np.tensordot(LocalHamiltonian(h).reshape(2,2,2,2),psireshape(2,2**(N-2),2),axes=[[2,3],[2,0]]).transpose(1,2,0).reshape(2**N)
    
    return psiOut 

def LinearHamiltonian(N, h=1):
    
    H=spsl.LinearOperator((2**N,2**N, matvec = ActWithH) 
    return H

Ns2 = np.arange(3,17,1)
h=1
timcompsp = []
timcomplo = []
for N in Ns2:
    print(f"N={N}")
    resultsp = %timeit -o spsl.eigsh(SparseHamiltonian(N))
    timcomsp.append(resultsp)
    resultslo = %timeit -o spsl.eigsh(LinearHamiltonian(N))
    timcomlo.append(resultlo)
timeavsp = [time.average for time in timcompsp]
timeavglo = [time.average for time in timcomplo]
plt.yscale('log')
plt.legend()
plt.show()



    
