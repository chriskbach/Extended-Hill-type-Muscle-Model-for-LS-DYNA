Provided is a simple arm model with muscle beam elements using an Extended Hill-type Muscle Material (EHTM) implemented as lsdyna user material.The aim of the model us to allow users to test the functionalities of their self-compiled version of the EHTM v2.0 (https://github.com/chriskbach/Extended-Hill-type-Muscle-Model-for-LS-DYNA/tree/v2.0). 
The muscle elements are routed using Part_Averaged keyword and are controlled via Equilibrium Point control. After initialization, the model holds Equilibrium position No.1 until at t=1.5s the model is moved to a more flexed position (EP-pos No.2). 


