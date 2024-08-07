import qiskit
from qiskit import QuantumCircuit
import matplotlib.pyplot as plt
import networkx as nx



class CNOTCircuit:


    def __init__(self, num_qubit, T) -> None:
        self._num_qubit=num_qubit
        self._T=T
        self._gateList=[]
        self._qiskitcircuit=QuantumCircuit(num_qubit)
        self._EPSTG=nx.Graph()
        self._RSG=nx.Graph()
        

    #Add a CNOT gate CNOT(control,target) at time
    def add_CNOT(self,control,target,time):
        self._gateList.append((control,target,time))


    def sort_gate_by_time(self):
        self._gateList=sorted(self._gateList, key=lambda x: x[2])


    def construct_EPSTG(self):
        for t in range(0,self._T):
            for qindex in range(0,self._num_qubit):
                self._EPSTG.add_node("Q"+qindex+"["+t+"]")
        pass


    def construct_RSG(self):
        pass

    
    def construct_matrix(self):
        pass


    def construct_qiskit_circuit(self):
        self.sort_gate_by_time()
        currenttime=0
        for (control,target,time) in self._gateList:
            if(time>currenttime):
                self._qiskitcircuit.barrier(label=str(time))
            self._qiskitcircuit.cx(control,target)


    def show_circuit(self):
        self._qiskitcircuit.draw(output="mpl")
        plt.show()


    def calculate_distribution():
        pass

    


    

