
# McCabe-Thiele Method with Murphree Efficiency


## McCabe-Thiele Method
The method was first published by Warren L. McCabe and Ernest Thiele in 1925, and is a graphical procedure of determining the number of trays within a distillation collumn. The Murphree tray efficiency modifies the equilibrium curve of the two components in order to better account for non-idealities. 

### Current implimentation
The current implementation presented here allows for a Murphree efficiency to be applied to the classical McCabe-Thiele plot. Whilst it is relatively easy online to find McCabe-Thiele 'calculators', there are few that have implementations accounting for the Murphree plate efficiency. 

### Limitations
* The molar heats of vaporization of the feed components are equal
* For every mole of liquid vaporized, a mole of vapor is condensed
* Heat effects such as heats of solution are negligible

## Function Use
```
    DESCRIPTION: 
    Performs the McCabe-Thiele construction in order to calculate
    optimum number of stages, and optimum feed stage. Also taking into 
    account the Murphree Efficiency of the system. 

    INPUTS: 
    PaVap       :Vapour pressure of component a (more volatile)
    PbVap       :Vapour pressure of component b (less volatile)
    R_factor    :Amount Rmin is scaled by to obtain the actual reflux ratio
    xf          :Feed composition 
    xd          :Distillate composition 
    xb          :Bottoms composition 
    q           :Liquid fraction of feed
    nm          :Murphree Efficiency

    OUTPUTS: 
    A McCabe-Thiele plot, displaying optimum number of equilibrium stages, 
    optimum feed stage, actual reflux ratio, actual bottoms composition. 
```
## Example
Runnning the following: 
```
PaVap=179.2     # Vapour pressure of a
PbVap=74.3      # Vapour pressure of b 
xd=0.975        # Distillate composition 
xb=0.025        # Bottoms composition 
xf=0.5          # Feed composition 
q=0.5           # q (liquid fraction of feed)
R_factor=1.8    # Reflux ratio = R_min* R_factor
nm=0.75         #Murphree tray efficiency

McCabeThiele(PaVap,PbVap,R_factor,xf,xd,xb,q,nm)
```
Produces the following output: 

<img align='center' src="https://github.com/TomRSavage/McCabeThiele/blob/master/McCabeThielePlot.png" width="500">

## Authors

* **Tom Savage** - *Initial work* - [TomRSavage](https://github.com/TomRSavage)
