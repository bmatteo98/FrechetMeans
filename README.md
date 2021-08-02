# Curvature:
We investigate the curvature of the tropical metric space.
  ## Curvature2D: 
    This code computes the curvature of a tropical triangle with bidimensional coordinates in the tropical projective torus R^3/R1.
    
  ## CurvatureND:
    This code computes the curvature of a tropical triangle with N-dimensional coordinates.
    
  ## typeFunction:
    This code returns the type of a given tropical triangle with bidimensional coordinates.
    
  ## sampleUltrametric_Rudy:
    This code uses Rudy's code (see References) to generate Ultrametric trees then I sample from the generated batch and compute the curvature.
    
  ## samplingSimplex:
    Sampling 3 bidimensional points from the simplex N times and compute the type of sampled triangles.

  ## samplingTriangles:
    Sampling triangles with random coordinates, integer coordinates and also triangles from both the tree space and ultrametric tree space computing the curvature and comparing results of curvatureND and curvature2D which should give the same result on R^3/R1.
    
  ## samplingTypes:
     Sampling type 1 and type 5 triangles with bidimensional coordinates and compute the curvature.

# Bo'sCode:
We translate the code provided by Bo Lin to compute Fréchet means of a given set of points in tropical metric space into R.

# SturmBacak:
We implement Sturm's and Bacak's algorithms to show that the do not converge in the tropical metric space.
