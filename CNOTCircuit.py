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
        self._EPSTG=nx.DiGraph()
        self._EPSTGPos={}
        self._RSG=nx.DiGraph()
        self._RSGPos={}

    #Add a CNOT gate CNOT(control,target) at time
    def add_CNOT(self,control,target,time):
        self._gateList.append((control,target,time))


    def sort_gate_by_time(self):
        self._gateList=sorted(self._gateList, key=lambda x: x[2])


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


    def show_EPSTG(self):
        nx.draw(self._EPSTG,with_labels = True,pos=self._EPSTGPos)
        plt.show()
    
    def show_RSG(self):
        nx.draw(self._RSG,with_labels = True,pos=self._RSGPos)
        plt.show()

    
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

    


    
