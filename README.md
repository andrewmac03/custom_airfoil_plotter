# custom_airfoil_plotter
Used to plot, export, and evaluate the performance of custom, user specified airfoils.


xfoil_runner is a modified version of xfoil-runner made by JARC99 who's github is linked below:
https://github.com/JARC99/xfoil-runner/tree/main/xfoil_runner


## To use the airfoil plotter import airfoil from airfoil_plotter


call airfoil to generate an airfoil object, Specify the following parameters:

-num_points: roughly the number of points along the airfoil's curve (xfoils can only do less than 320 points)

-chord_length: the length from the leading edge to the trailing edge

-max_camber: the max distance of the camber line to the chord line as a percentage of the chord length

-max_camber_pos: where the maximum camber is located along the chord as a percentage of the chord length

-max_thickness: the max thickness as a percentage of the chord length

-max_thickness_pos: where the maximum thickness is located along the chord as a percentage of the chord length

-leading_edge_radius: the radius of the leading edge as a percentage of the chord length

-TE_sharpness: 1 for a zero thickness TE, and 0 for a 45 degree TE



## To plot the airfoil call the plot_airfoil method with no arguments


To export your airfoil to a .txt file call the export_airfoil method with the following arguments:

-name: the name of the airfoil which will be used to name the file

-x_foils=False: if true, the airfoil will be exported in a format best suited for xfoils, if fales the airfoil will be exported to best suit solidworks


## To run xfoils on your airfoil import run_xfoil from xfoil_runner


callrun_xfoil with the following arguments:

-name: the name of the file (same as name in export_airfoil)

-alpha_i: the first angle of attack you want to simulate in degrees

-alpha_f: the last angle of attack you want to simulate in degrees

-alpha_step the angle of attack step in degrees to sweep between the first and last angle of attack

-Re: the reynolds number to be simulated at

-n_iter: the number of iterations for the simulation to use


x_foils will  create a polar file with your results.

If you call the function as the definition of a variable then it will set the variable to a numpy array with the results

