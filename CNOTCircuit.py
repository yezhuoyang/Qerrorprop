from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit import QuantumCircuit, transpile
from qiskit_aer.noise import (NoiseModel, QuantumError, ReadoutError,
    pauli_error, depolarizing_error, thermal_relaxation_error)
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import random
from scipy.stats import bernoulli

class CNOTCircuit:


    def __init__(self, num_qubit, T, prob) -> None:
        self._num_qubit=num_qubit
        self._T=T
        self._gateList=[]
        self._qiskitcircuit=QuantumCircuit(num_qubit)
        self._EPSTG=nx.DiGraph()
        self._EPSTGPos={}
        self._RSG=nx.DiGraph()
        self._RSGPos={}
        self._p=prob
        self._noise_model=NoiseModel()
        self._noise_gate=pauli_error([('X',  self._p), ('I', 1 - self._p)])
        self._noise_model.add_all_qubit_quantum_error(self._noise_gate.tensor(self._noise_gate), ["cx"])
        self._qiskit_result=None
        self._M=None
        self._dist={}
        self._qiskit_dist={}
        self._entropy=None

    #Add a CNOT gate CNOT(control,target) at time
    def add_CNOT(self,control,target,time):
        self._gateList.append((control,target,time))


    def sort_gate_by_time(self):
        self._gateList=sorted(self._gateList, key=lambda x: x[2])



    def calculate_entropy(self):
        non_zero_probs = np.array(list(filter(lambda p: p > 0, self._dist.values())))
        self._entropy =  -np.sum(non_zero_probs * np.log2(non_zero_probs))
        return self._entropy


    def construct_all(self):
        self.construct_EPSTG()
        self.construct_RSG()
        self.construct_matrix()
        self.construct_qiskit_circuit()


    def construct_EPSTG(self):
        for t in range(0,self._T):
            for qindex in range(0,self._num_qubit):
                self._EPSTG.add_node("Q"+str(qindex)+"["+str(t)+"]")
                self._EPSTGPos["Q"+str(qindex)+"["+str(t)+"]"]=(qindex,t)
        for t in range(0,self._T-1):
            for qindex in range(0,self._num_qubit):
                filtered_gate = list(filter(lambda x: (x[0]==qindex) and (x[2]==t), self._gateList))
                self._EPSTG.add_edge("Q"+str(qindex)+"["+str(t)+"]","Q"+str(qindex)+"["+str(t+1)+"]") 
                for (control,target,time) in filtered_gate:
                    self._EPSTG.add_edge("Q"+str(control)+"["+str(time)+"]","Q"+str(target)+"["+str(time+1)+"]")                    
            

    def construct_RSG(self):
        for t in range(0,self._T-1):
            for qindex in range(0,self._num_qubit):
                self._RSG.add_node("Q"+str(qindex)+"["+str(t)+"]")    
                self._RSGPos["Q"+str(qindex)+"["+str(t)+"]"]=(0,t*self._num_qubit+qindex)    
        for qindex in range(0,self._num_qubit):
            self._RSG.add_node("Q"+str(qindex)+"["+str(self._T-1)+"]")              
            self._RSGPos["Q"+str(qindex)+"["+str(self._T-1)+"]"]=(4,qindex)  
        for t in range(0,self._T-1):
            for qindex1 in range(0,self._num_qubit):
                for qindex2 in range(0,self._num_qubit):
                    if self.EPSTG_has_path(qindex1,t,qindex2):
                        self._RSG.add_edge("Q"+str(qindex1)+"["+str(t)+"]","Q"+str(qindex2)+"["+str(self._T-1)+"]") 



    #Whether there is a path from Q1[t1] to Q2[T-1]
    def EPSTG_has_path(self,qindex1,t1,qindex2):
        return nx.has_path(self._EPSTG,"Q"+str(qindex1)+"["+str(t1)+"]","Q"+str(qindex2)+"["+str(self._T-1)+"]")


    def RSG_has_edge(self,qindex1,t1,qindex2):
        return self._RSG.has_edge("Q"+str(qindex1)+"["+str(t1)+"]","Q"+str(qindex2)+"["+str(self._T-1)+"]")


    def show_EPSTG(self):
        nx.draw(self._EPSTG,with_labels = True,pos=self._EPSTGPos)
        plt.show()
    
    def show_RSG(self):
        nx.draw(self._RSG,with_labels = True,pos=self._RSGPos)
        plt.show()

    
    def construct_matrix(self):
        rownum=self._num_qubit
        colnum=self._num_qubit*(self._T-1)
        self._M=np.zeros((rownum,colnum))
        for t in range(0,self._T-1):
            for qindex1 in range(0,self._num_qubit):
                for qindex2 in range(0,self._num_qubit):
                    rowindex=qindex2
                    colindex=t*self._num_qubit+qindex1
                    if self.RSG_has_edge(qindex1,t,qindex2):
                        self._M[rowindex][colindex]=1
        

    #Convert a binary integer to a vector 
    def vec_from_integer(self, integer):
        dim=self._num_qubit*(self._T-1)
        array=np.zeros((dim,1))
        for i in range(0,dim):
            array[i][0]=integer%2
            integer=integer>>1
        return array
    
    def vec_from_bitString(self,bitstring):
        dim=self._num_qubit*(self._T-1)
        array=np.zeros((dim,1))
        for i in range(0,dim):
            if(bitstring[i]==1):
                array[i][0]=1
        return array

    #Calculate the probability of a certain error pattern    
    def calc_prob(self,vector):
        numones=np.sum(vector)
        return self._p**(numones)*(1-self._p)**(self._num_qubit*(self._T-1)-numones)


    def show_matrix(self):
        print(self._M)


    def calculate_syndrome(self,vector):
        syndrome_vec=np.matmul(self._M,vector)%2
        return syndrome_vec


    def construct_qiskit_circuit(self):
        self.sort_gate_by_time()
        currenttime=0
        for (control,target,time) in self._gateList:
            if(time>currenttime):
                self._qiskitcircuit.barrier(label=str(time))
                currenttime=time
            self._qiskitcircuit.cx(control,target)
        self._qiskitcircuit.measure_all()


    def show_circuit(self):
        self._qiskitcircuit.draw(output="mpl")
        plt.show()


    def calculate_distribution_exact(self):
        distribution={i:0 for i in range(0,self._num_qubit+1)}
        num_source=self._num_qubit*(self._T-1)
        for i in range(0,1<<num_source):
            tmpvec=self.vec_from_integer(i)
            syndrome_vec=self.calculate_syndrome(tmpvec)
            syndromecount=int(np.sum(syndrome_vec))
            distribution[syndromecount]+=self.calc_prob(tmpvec)
        self._dist=distribution
        return distribution


    def calculate_distribution_sample(self,Nsample):
        distribution={i:0 for i in range(0,self._num_qubit+1)}
        num_source=self._num_qubit*(self._T-1)
        colnum=self._num_qubit*(self._T-1)
        for n in range(0,Nsample):
            bit_string = list(bernoulli.rvs(self._p, size=colnum))
            tmpvec=self.vec_from_bitString(bit_string)
            syndrome_vec=self.calculate_syndrome(tmpvec)
            syndromecount=int(np.sum(syndrome_vec))
            distribution[syndromecount]+=1
        distribution={i:distribution[i]/Nsample for i in range(0,self._num_qubit+1)}
        self._dist=distribution
        return distribution
    
    def run_qiskit(self,shots):
        simulator=AerSimulator(noise_model=self._noise_model,shots=shots)
        
        self._qiskit_result = simulator.run(self._qiskitcircuit).result()
        counts = self._qiskit_result.get_counts()
        self._qiskit_dist={i:0 for i in range(0,self._num_qubit+1)}
        for key in counts.keys():
            self._qiskit_dist[key.count('1')]+=counts[key]
        return self._qiskit_dist


    def show_distribution(self):
        # Extract keys and values
        keys = list(self._dist.keys())
        values = list(self._dist.values())

        # Plotting the distribution
        plt.figure(figsize=(8, 6))
        plt.bar(keys, values, color='skyblue')
        plt.xlabel('Number of errors')
        plt.ylabel('Probability')
        plt.title('Distribution of error numbers')
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        # Display the plot
        plt.show()


    def get_expectation_and_std(self):
        expectation=0
        expectation_square=0
        for i in range(0,self._num_qubit+1):
            expectation+=i*self._dist[i]
            expectation_square+=i*i*self._dist[i]
        std=np.sqrt(expectation_square-expectation*expectation)
        return (expectation,std)



#Generate a random circuit. In each time window, there are exactly gate_pair_each_T pairs of CNOT gates chosen randomly
def random_circuit(n_qubits,T,gate_pair_each_T,p):
    circuit=CNOTCircuit(n_qubits,T,p)
    numbers = list(range(0, n_qubits))
    for t in range(0,T):
        random.shuffle(numbers)
        for i in range(0,gate_pair_each_T):
            control=numbers[2*i]
            target=numbers[2*i+1]
            circuit.add_CNOT(control,target,t)
    return circuit



def transversal_circuit(n_qubits,T,gate_pair_each_T,p):
    assert n_qubits%2==0
    circuit=CNOTCircuit(n_qubits,T,p)
    r=n_qubits//2
    current_index=0
    for t in range(0,T):
        if(current_index+r<n_qubits):
            for i in range(0,gate_pair_each_T):
                if(current_index+r<n_qubits):
                    circuit.add_CNOT(current_index,current_index+r,t)
                    current_index=current_index+1
    return circuit
    
    

def transversal_parallel_circuit(n_qubits,p):
    assert n_qubits%2==0
    circuit=CNOTCircuit(n_qubits,3,p)
    pass

    
    
    
    
def fully_connected_circuit(n_qubits,T,p):
    circuit=CNOTCircuit(n_qubits,T,p)
    for t in range(0,T-1):
        for ctrl_index in range(0,n_qubits-1):
            for targ_index in range(ctrl_index+1,n_qubits):
                circuit.add_CNOT(ctrl_index ,targ_index ,t)
        for ctrl_index in range(n_qubits-1,0,-1):
            for targ_index in range(ctrl_index-1,-1,-1):
                circuit.add_CNOT(ctrl_index ,targ_index ,t)           
    return circuit
    
    


    
    


