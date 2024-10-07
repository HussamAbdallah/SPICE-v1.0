# SPICE-v1.0

"SPICE v1.0" is a circuit simulator that handles passive elements (Resistors - Capacitors - Inductors - Diodes) and active elements (Voltage controlled current sources) as well as independents sources (Voltage sources - Current sources). "SPICE v1.0" supports different type of analyses (Operating point - DC Sweep - AC) analysis. The following flowchart desribes the main main program workflow. 
<p align="center">
  <img src="https://github.com/user-attachments/assets/67d876e0-1ba7-4b05-9ae2-e125daf7dff2" alt="Main Flow" width="400"/>
</p>



First the netlist in the program directory is parsed into a map where keys represent the symbols like "C1" for caps the value of each element in the map is a vector holding cap nodes and capacitance value. Then the system matricies are formulated according to each element stamp in the modified nodal analysis approach.


The following flowchart represents how OP Analysis is conducted.

<p align="center">
  <img src="https://github.com/user-attachments/assets/c68bf71a-0f38-4a5a-b68f-e867f625743f" alt="OP Analysis" width="400"/>
</p>

Since we can only solve linear matrices, if the circuit contains diodes, the diodes has to be linearized first i.e., replacing the diode with its companion model shown in the following figure. The Newton-Raphson is then used to find the operating point by linearizing the diode using the initial guess, then solve the resultting linear circuit to get a better solution point and finally iterate until convergence.  
<p align="center">
  <img src="https://github.com/user-attachments/assets/a10614f9-e497-4d7d-be86-fe577a731f46" alt="Diode Companion Model" width="400"/>
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
  <img src="https://github.com/user-attachments/assets/5f7d6180-5ea0-4b3f-a02b-cd3a3931ed25" alt="DC_Sweep" width="400"/>
</p>


Finally, The AC Analysis is similar to DC Sweep except that there is no need to run Newton Method for every sweeped frequency, instead the Newton method is applied once so that the ac-model of the diode which is just $G_{eq}$ is obtained and then linearized circuit is sweeped with frequency.

![Running Program](https://github.com/user-attachments/assets/45698d34-a5e2-4d4f-837b-7441b4bb3871)

To verify the Solution, The following circuit is simulated using both LTSPICE and our SPICE v1.0. The results are perfectly matched.

VD 1 0 AC 1

R1 1 2 1e3

C1 2 3 100e-3

L1 3 4 10e-3

C2 4 0 100e-3

D1 4 0 1mA_diode

* diode model statement

.model 1mA_diode D (Is=100pA n=1.679)

![Running Program](https://github.com/user-attachments/assets/7350a3a9-b287-4ea7-a32f-9e4a759ccc3c)



