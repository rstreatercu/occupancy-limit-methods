# occupancy-limit-methods
 Simple simulation of binding and unbinding crosslinking proteins to sites (receptors, etc.) where both sites and proteins may only bind once. In this test, three sites are placed along a circle (radius 1) at 10 degree increments, while two crosslinkers are placed at the xy coordinates (0.9, 0) and (0.7, 0.2):
 
```
|
|
|
|           x
|
|______________x_s
|
|               s
|
|            s
|
```
 
 The sites and crosslinkers have fixed positions throughout the simulation. A binding event between a crosslinker i and a site j increases the potential energy of the system by E<sub>ij</sub>=k r<sub>ij</sub><sup>2</sup>-E<sub>0</sub>, where k is the crosslinker spring constant, r<sub>ij</sub> is the distance between crosslinker and site, and E<sub>0</sub> is a binding energy. The binding probabilities follow a Boltzmann distribution, where β=1 and probability of binding is p<sub>bind</sub> = k<sub>on</sub> dt exp(-E<sub>ij</sub>/2) and probability of unbinding is p<sub>unbind</sub> = k<sub>on</sub> dt exp(E<sub>ij</sub>/2).
 
 In the simulation, sites and crosslinkers bind and unbind according to a uniform random roll, but cannot bind if the given site or crosslinker is already bound. There is some question of how to handle situations where multiple binding events are probable in a single timestep. This simulation is intended to measure the accuracy of several methods of collision resolution as timestep becomes large.
