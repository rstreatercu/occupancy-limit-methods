# occupancy-limit-methods
## Overview
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
 
## Build and run
 The simulation code is minimal and just exists in one small file, so there is not currently an installer or input file. To build and run the simulation, you will need a compiler that supports the C++11 standard, and compile and link test_occupancy.cpp with the flag:
 ```
  -std=c++11
 ```
 Running as is should give five output text files for the five collision resolution methods. To analyze the output data, you will need Wolfram Mathematica. Load and evaluate the notebook test_occupancy.nb, changing the variable "dtRange" to import data with different dt values.

## Description of the five collision resolution methods

 The methods tested are determined when initializing a simulation by the simType parameter. The options are:
 * simType 0: Collisions allowed. This mode allows multiple crosslinkers and sites to bind **in a single timestep**. The collisionCount, or number of crosslinkers and sites that bound multiple times in a single timestep, is printed. Single occupancy is still enforced for events in different timesteps. This simType has the following uses:
    * It is a negative control- when comparing simulation outputs, this simType should behave incorrectly as timestep increases.
    * The collision count provides a measure of "large" vs. "small" timestep dt. If the number of collisions is much lower than the number of simulation runs, the dt value in question is too small to show any effects of the different collision resolution methods.
    * This simType provides a "ground truth" output when dt is lowered enough to make the collision count 0. This can be compared to theory to determine whether the simulation length tmax_ is large enough that the simulation is essentially at equilibrium. Using this simType with a small enough dt that collision count is 0 also provides a ground truth output for non-equilibrium simulations, where there is no simple theoretical solution.
 * simType 1: First-come-first-serve. If a crosslink or site is bound in a timestep, that crosslink or site may not bind in the same timestep.
    * This is the naive collision-resolution solution, and is the simplest to implement. It may introduce a bias towards binding the first crosslinkers and sites listed.
 * simType 2: Random shuffle. In this simType, the sites are shuffled in order each timestep to avoid first-come-first-serve bias.
    * This simType requires the greatest amount of random numbers generated (N<sub>i</sub>N<sub>j</sub> for the event rolls + 1 for the shuffle each dt)
 * simType 3: Whether/which method. In this method, two random numbers are used to roll **whether** a binding event occurs, and then **which** binding event occurs. All of the probabilities of binding in a timestep are saved and summed up. Total probability of having some binding event is (1-exp(-p<sub>total</sub>)). If a binding event is rolled, then the individual binding event probabilities are normallized and a roll determines **which** to bind to. A major flaw with this current implementation (and using whether/which in two directions) is that only one binding event occurs, when there really could be up to two in a timestep. It could be used properly in combination with other methods.
    * This simType requires the least random numbers, but its accuracy decreases when total probability ~ 1
 * simType 4: Knockout. In this method, collisions within a timestep are allowed temporarily. Then, if there are collisions, the code will go through conflicting sites and crosslinkers and roll for a "winner" based on probability. The remaining conflicting crosslinkers/sites will unbind.
    * This simType can generally be implemented efficiently by tracking conflicts with crosslinker and site objects themselves, where a knockout event occurs infrequenctly. However, the knockout loop has a disadvantage when dealing with sites and crosslinkers **at once**. It can give incorrectly low % time bound values, for example, in the following situation:
    <table>
       <tr>
        <td></td>
        <th>Site1</th>
        <th>Site2</th>
        <th>Site3</th>
       </tr>
       <tr>
        <th>Xlink1</th>
        <td>Wants to bind</td>
        <td>Wants to bind</td>
        <td></td>
       </tr>
       <tr>
        <th>Xlink2</th>
        <td>Wants to bind</td>
        <td>Wants to bind</td>
        <td></td>
       </tr>
      </table>
      In the first knockout loop, the xlink1/site1 binding event and the xlink2/site1 binding event may win:
      <table>
       <tr>
        <td></td>
        <th>Site1</th>
        <th>Site2</th>
        <th>Site3</th>
       </tr>
       <tr>
        <th>Xlink1</th>
        <td>Wants to bind</td>
        <td><strike>Wants to bind</strike></td>
        <td></td>
       </tr>
       <tr>
        <th>Xlink2</th>
        <td>Wants to bind</td>
        <td></td>
        <td></td>
       </tr>
      </table>
      Then, say the xlink1/site1 binding event wins:
      <table>
       <tr>
        <td></td>
        <th>Site1</th>
        <th>Site2</th>
        <th>Site3</th>
       </tr>
       <tr>
        <th>Xlink1</th>
        <td>Wants to bind</td>
        <td></td>
        <td></td>
       </tr>
       <tr>
        <th>Xlink2</th>
        <td><strike>Wants to bind</strike></td>
        <td></td>
        <td></td>
       </tr>
      </table>
      Only one binding event happens, where there could have been two.
