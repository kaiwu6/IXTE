Introduction of the structure

BSE_module: 
  Flow: 
  1. band strucrture of single layer and its k-space representation
  2. Generate and diagonalize the BSE matrix to get the IX energy levels and wavefunctions

  Functions:
  1. Class for single layer band structure
  2. Generate and diagonalize the matrix for each Q of IX
  3. Plot wavefunctions to show the localization of the images in real- or k- space
  4. Plot high-symmetrical points of band structure for illustration
  5. dump eigenstates for record and further use?

Fit_module:
  Purposes: fit the IX band to get tight-banding model parameters.

  Functions:
  1. band_fitting: get hopping parameters
    input: energies from BSE euqation
    ouput: hopping parameters
  2. band_expansion: 
    input: k_points (list of array), hopping parameter(unpacked)
    output: evaluated energies(list)

Transport_module:
    Purposes: calculate all the transport elements for the IX

    Flow: 

    Functions:

Instruction of usages:
1. model_only.py 
    Diagonalize BSE to get exciton band structure and fit for the tight-binding model.
    output: 2 bands in bands.dat, hopping parameters in t_1.dat and t_2.dat
    figure: 
        scatter plot for first and second band
        fitting plot + scatter plot
        do we need 3-d plot ?

2. tight_binding.py 
    We can run tight_binding directly to get the hopping parameters. It requires the calculated band structure. 

3. trans_only.py
    Calculate all transport coefficients with certain hopping parameters.
    The chemical potential is hard-coded in and need to tune by hand to be fitted in the range for plotting. 
    figures need to desgin




