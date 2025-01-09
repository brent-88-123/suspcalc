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
* Scrub radius
* Mechanical trail
* Roll centre (2D-centreline)

Next steps:

* Add reference frame definition
* Mirror points to LHS and properly estimate RC
* Calculate KPI/Caster/SAL/ANTI GEO, etc
* Introduce reaction force estimation:
    * Tire Fz load reaction
    * Tire Fy load reaction
    * Tire Fx load reaction
