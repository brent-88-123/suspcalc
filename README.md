# suspcalc
Vehicle suspension kinematic estimation.

Inputs:
* Double a-arm suspension joint locations [x,y,z]
* Static tire centre of pressure (force application point) [x,y,z]
* Tire force vector
* Sprung mass movements:
    * Roll
    * Jounce

Outputs:
* Link points 
* Roll centre (2D-centreline)
* KPI
* Caster
* Scrub radius
* Mechancial trail


Next steps:

* Add reference frame definition
* Mirror points to LHS and properly estimate RC
* Introduce reaction force estimation:
    * Tire Fz load reaction
    * Tire Fy load reaction
    * Tire Fx load reaction
