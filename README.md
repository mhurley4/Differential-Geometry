# Differential-Geometry

The curvatures.py file contains several functions, each of which takes an input a parameterization x and parameters u, v.

The functions are:
tangent_vectors()
  Computes tangent vectors of surface
unit_normal
  Compute unit normal vector (Gauss map) of surface
first_ff
  Computes first fundamental form of surface
second_ff
  Compute second fundamental form of surface
dN_p
  Computes differential of the Gauss map in matrix form
Gaussian_curvature
  Computes Gaussian curvature of surface
mean_curvature
  Computes mean curvature of surface
principal_curvature
  Computes principal curvatures of surface
compute_curvatures
  All-in-one function that computes (and prints) first and second fundamental forms, Gaussian and mean curvatures,
    and principal curvatures.
   
Compute_Curvatures.ipnyb is a Jupyter notebook that can be used to implement functions in curvatures.py.
It imports curvatures.py and then uses functions from there.
Currently computes curvatures for three surfaces (two helicoids and sphere).
