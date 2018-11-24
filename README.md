The goal of this project is to find DC equilibrium points of circuits using interval verification algorithms like the Krawczyk operator. Our main example is Rambus ring oscillator challenge problem - http://www.cs.um.edu.mt/gordon.pace/Workshops/DCC2008/Presentations/10.pdf.
We also find DC equilibrium points for the Schmitt Trigger.

All of the implementations of our programs have been done in python2.7 and require numpy and cvxopt packages.

# Overall Structure

Our main main method is implemented in **prototype.py**. It gives a list of hyperrectangles each of which contains a unique DC equlibrium point.

Note that this method will only work properly with well-conditioned systems where the jacobians are not singular or badly conditioned at the solutions.

Our main algorithm is implemented in the function `solverLoopNoLp` in **prototype.py**. The two main parameters for this function is uniqueHypers (which should be an empty list when the function is called) and a model object. uniqueHypers will hold the list of hyperrectangles containing unique solution to a function f. The model class should have a method called f which gives point and interval evaluation of the function we are trying to solve. f should take a numpy array of variable values. If any entry in the array is an interval, f returns an interval value. Otherwise it returns a point. Similarly, the model class should also have a method called jacobian that returns jacobian and jacobian interval depending on the arguments. The model class should also have a field called bounds that specifies bounds for each variable. Look at `RambusTanh` class in **circuitModels.py** for more details. After having defined the model class, in prototype.py it is enough to write

`uniqueHypers = []`

`solverLoop(uniqueHypers, model)` 

to find all the hyperrectangles containing unique solutions. uniqueHypers will contain them.

# Defining Custom Circuit

It is also possible to define a custom circuit with long channel and short channel MOSFET models and use our method to find DC equilibrium points. 

The `Circuit` class in **circuit.py** takes a list of transistors. Each transistor can be a short channel transistor or long channel transistor. The transistor object will take in indices for source, drain and gate from an array of voltages. We can specify these requirements in a model class and extract the required function and jacobian information from the `Circuit` object. An example of creating and using a circuit with long channel or short channel MOSFET in a model object can be seen in `RambusMosfet` class in **circuitModels.py**.

# Running Our Experiments:

To run the experiments related to rambus oscillator and schmitt trigger, `cd` into experiments directory. 

To run the experiments related to our method, run **solverExp.py**.

To run the dReal experiments, run **dRealExp.py**. Follow the instructions in https://github.com/dreal/dreal4 to install dReal with python bindings.

To run the z3 experiments, run **z3Exp.py**. Follow the instructions in https://github.com/Z3Prover/z3 to install Z3 with python bindings.

Note that we kill the process running an example if it takes more than 10 hours to run. This is what we consider to be timeout.