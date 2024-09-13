# SPICE-v1.0

"SPICE v1.0" is a circuit simulator that handles passive elements (Resistors - Capacitors - Inductors - Diodes) and active elements (Voltage controlled current sources) as well as independents sources (Voltage sources - Current sources). "SPICE v1.0" supports different type of analyses (Operating point - DC Sweep - AC) analysis. The following flowchart desribes the main main program workflow. 
<p align="center">
  <img src="https://i.postimg.cc/7LqSLRKZ/main-workflow.png" alt="Main Flow" width="400"/>
</p>
First the netlist in the program directory is parsed into a map where keys represents the symbols like "C1" for caps the value of each element in the map is a vector holding cap nodes and capacitance value. Then the system matricies are formulated according to each element stamp in the modified nodal analysis approach.


The following flowchart represents how OP Analysis is conducted.

<p align="center">
  <img src="https://i.postimg.cc/rmCGjkzW/OP-Analysis.png" alt="OP Analysis" width="400"/>
</p>

Since we can only solve linear matrices, if the circuit contains diodes, the diodes has to be linearized first i.e., replacing the diode with its companion model shown in the following figure. The Newton-Raphson is then used to find the operating point by linearizing the diode using the initial guess, then solve the resultting linear circuit to get a better solution point and finally iterate until convergence.  
<p align="center">
  <img src="https://i.postimg.cc/hGGh6RCT/image.png" alt="Diode Companion Model" width="400"/>
</p>



$ID=I_{sat}[exp(\frac{V_D}{nV_T})-1]$ 
    


$Geq=I_D'=\frac{I_{sat}}{nV_T}exp(\frac{V_D}{nV_T})$

Newton Approximation $f(x)=f(x_0)+f(x_0)'\delta x$

Then 

$I_D=I_D^{k} + G_{eq}^{k} (V-V_D^{k})$ 

$I_D=G_{eq}^{k} V + I_{eq}^{k}$

Then 

$I_{eq}^{k}=I_D^{k} - G_{eq}^{k} V_D^{k}$


From the initial guess, we have initial $V_D$ that will be used to determine initial $G_{eq}$ and $I_{eq}$ then solving the resultatnt linear circuit, getting a better guess for $V_D$ and consequently for $G_{eq}$ , $I_{eq}$.

For DC Sweep Analysis, the symbolic solution is obtained from OP Analysis and then the source value is sweeped according to user input. If the circuit contains non-linear elements the newton methods uses the previous sweeped solution point as the new guess for the next sweep point. The following flowchart illustrates how DC Sweep is conducted.

<p align="center">
  <img src="https://i.postimg.cc/5ywfWMFB/DC-SWEEP.png" alt="Diode Companion Model" width="400"/>
</p>


Finally, The AC Analysis is similar to DC Sweep except that there is no need to run Newton Method for every sweeped frequency, instead the Newton method is applied once so that the ac-model of the diode which is just $G_{eq}$ is obtained and then linearized circuit is sweeped with frequency.


