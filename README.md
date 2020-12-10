# Voids_fR_public
This github hosts the void catalogs, and analysis codes used in the work by Chris Wilson and Rachel Bean at "ARXIV LINK HERE".


All analysis code was developed by Wilson and Bean. The void catalogs were created using the VIDE void finder on the ELEPHANT simulations, provided courtesy of Baojiu Li. If using these void catalogs, please cite the following papers, as well as our paper listen above.

Sutter et al. https://arxiv.org/abs/1406.1191

Cautun et al. https://arxiv.org/abs/1710.01730

Li et al. https://arxiv.org/abs/2011.05771


#####################################################################################################################


these voids are saved as numpy object array files. use np.load("path_to_file", allow_pickle=True) to load them. sometimes, encoding="latin1" must also be passed depending on the version of python.
Example:

np.load(GR_partVoidPath_z0p5,allow_pickle=True,encoding="latin1")

where here I have saved GR_partVoidPath_z0p5 as a string variable

file naming convention:
GR,F6,F5 - self explanitory

haloVoids - voids whose density and Velocty profile information was constructed from the halo data
partVoids - voids whose density and Velocty profile information was constructed from the particle data
z0p0 - at redshift 0
z0p5 - at redshift 0.5

void centers are always identified with halos. the number of a certain void are the same across the corresponding "haloVoids" and "partVoids" catalog. For example, F6_partVoids_z0p5[345] and F6_haloVoids_z0p5[345] refer to the same physical void, with velocity and density information coming from either the particles or halos respectively.

voids have the following attributes.
if cat=np.load(GR_partVoidPath_z0p5,allow_pickle=True,encoding="latin1") and cat[i]=void, then

void[0][0]=realization number
void[0][1]= void ID number (largely un-used)
void[1]="R" or "S". this field gets changed when the reclassifyVoidRS() method gets run. since these are object arrays, no issue with classifying with a number instead, although many other methods will then require updating.

void[2]= R_{eff} with units of Mpc/h
void[3][0] = density profile, bin centers, in units of r/R_{eff}
void[3][1] = density profile, values for delta

void[4][0] = radial velocity profile, bin centers in units of r/R_{eff}
void[4][1] = radial velocity profile, mean value in each bin, in units of km/s


if the void is a "haloVoid" then it will also have
void[0][2] = np.array([realization's mean number density]) literally an array of just one number cause i didnt catch this before i saved everything and too many things
void[0][3] = void macrocenter - array of 3 numbers , x,y,z positions.

We have also saved the fifth force calculated in every void, using the same naming convention as above but with "fifthConverged" as the name. The "fifthConverged" numpy arrays were too large to save in github withut complaining, so they are broken down into real_1, real_2, etc... for each of the 5 realizations. A code to re-construct them into a single numpy array is included in the example calculations notebook

We additionally saved the associated screening factors (alpha) for each void, evaluated at the peak of the fifth force profile, as well as the maximum error between the reconstructed density profile and the real density profile using the name "alphaErrors".


example
F6_fifthConverged_z0p0[i][0] = fifth force for the ith void in F6 at z=0.5. the associated radial values were not saved, but are always given as "xpoints=np.arange(0,5,.02)"
F6_fifthConverged_z0p0[i][0] = fifth field "                                                                                                                 "

F6_alphaErrors_z0p5[i][0] = screening factor (alpha) evaluated at the peak of the fifth force in the ith void in F6 at z=0.5
F6_alphaErrors_z0p5[i][1] = error between the re-constructed density and the true density in the ith void in F6 at z=0.5. These errors give an indication to the trustworthyness of the fifth force calculation. they can be used in conjunction with the "filterForces" method to filter through the "fifthConverged" arrays, selecting out only those forces with a low enough error.



Information about the useful methods should be provided in the "exampleCalculations" file.


For any questions please do not hesitate to email me at cww78@cornell.edu. I would be happy to help in any way possible
