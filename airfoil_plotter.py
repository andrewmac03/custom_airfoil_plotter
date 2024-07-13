import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sp


class airfoil:
    '''
    allows users to specify airfoil parameters, then generates an airfoil 
    profile that meets the paramters allowing for users to plot the airfoils,
    and export the profile in a format acceptable for soildworks or xfoils. 
    '''
    def __init__(self, num_points, chord_length, max_camber, max_camber_pos, 
                 max_thickness, max_thickness_pos, leading_edge_radius, 
                 TE_sharpness):
        self.chord_length = chord_length
        #since there are 2 surfaces divide number of points by 2
        self.num_points_export = int(num_points / 2)
        #set initial num points to be higher resolution then remove res later
        self.num_points = num_points * 10
        self.max_camber = max_camber
        self.max_camber_pos = max_camber_pos
        self.max_thickness = max_thickness
        self.max_thickness_pos = max_thickness_pos
        self.leading_edge_radius = leading_edge_radius
        #between 0 and 1 where 1 is a zero thickness TE and 1 is a 45 degree TE
        self.TE_sharpness = TE_sharpness
        self.x = np.arange(0, 1+1/self.num_points, 1/self.num_points)
        return
    
    def generate_camber_line(self):
        '''
        Generates camber line using NACA 4 digit equation for camber 
        line as given by: https://en.wikipedia.org/wiki/NACA_airfoil
        '''
        x = self.x
        p = self.max_camber_pos
        m = self.max_camber

        front_camber = (m / p**2) * (2 * p * x - x**2)
        back_camber = (m / (1 - p)**2) * ((1-2*p) + 2 * p * x - x**2)
        self.camber_y = np.append(front_camber[:int(self.num_points*p)], 
                                  back_camber[int(self.num_points*p):])
        front_camber_grad = 2 * m / p**2 *(p-x) 
        back_camber_grad = 2 * m / (1-p)**2 * (p-x)
        self.camber_grad = np.append(front_camber_grad[:int(self.num_points*p)], 
                                     back_camber_grad[int(self.num_points*p):])
        return


    def generate_thickness_profile(self):
        '''
        Generates thickness profile using the 
        equation given by equation on page 29 of:
        https://www.collectionscanada.ca/obj/s4/f2/dsk1/tape4/PQDD_0010/MQ59866.pdf
        Creates a system of equations for the coefficients in the 
        afformentioned equation that satisfy the boundary conditions 
        then solves the system of equations and plugs the coefficients 
        back in to generate the thickness profile

        Boundary conditions are:
        the max thickness location
        the thickness at the max thickness location
        the TE thickness
        the slope at the TE
        The radius of the LE
        '''
 
        t = self.max_thickness
        LE_r = self.leading_edge_radius
        x = sp.Symbol('x')
        a0, a1, a2, a3, a4 = sp.symbols('a:5')
        yt = 5*(a0*x**(1/2) - a1*x - a2*x**2 + a3*x**3 - a4*x**4)
    
        #BC of TE thickness of 0
        #specify first solution -> there should only be one solution
        eq1 = sp.Eq(a0, sp.solve(yt.subs(x,1), a0)[0])
    
        
        #BCs of LE shape
        #set thickness at a quarter circle to be that of a quarter 
        #circle w/ leading edge radius
        θ = np.pi/4
        x_quarter = LE_r * (1-np.cos(θ))
        y_quarter = x_quarter/(1-np.cos(θ))*np.sin(θ)
        #define a func that is y_t - y_quarter so it can be 0 at x_quarter
        y_solve_quarter = yt - y_quarter
        eq2 = sp.Eq(a1, sp.solve(y_solve_quarter.subs(x,x_quarter), a1)[0])

        y_solve_max_thickness = yt - t
        eq3= sp.Eq(a2, sp.solve(y_solve_max_thickness.subs(x, 
                                self.max_thickness_pos), a2)[0])
    
        
        derivative = yt.diff(x)
        
        #BC of sharp TE
        TE_derivative_solve = derivative + (1-self.TE_sharpness)
        eq4 = sp.Eq(a3, sp.solve(TE_derivative_solve.subs(x, 1), a3)[0])
        
        #BC of max thickness location
        eq5 = sp.Eq(a4, sp.solve(derivative.subs(x, self.max_thickness_pos),
                                 a4)[0])
    
        #find coefficients that satisfy all BCs
        #nsolve numerically solves so requires an input vector close to the 
        #solution, we can use the default thickness profile coefficients for this
        solutions = sp.solve([eq1, eq2, eq3, eq4, eq5], [a0, a1, a2, a3, a4])
        lam_f = sp.lambdify(x, yt.subs(solutions))
        vfunc = np.vectorize(lam_f)
        self.yt = vfunc(self.x)
        return
    

    def generate_airfoil(self):
        '''
        generates the airfoils by placing the upper and lower surfaces
        at some thickness away from the camber line.
        '''
        
        #need to get rid of overlapping points at LE or it 
        #won't import to SOLIDWORKS (whatever edge the curve doesn't start at)
        self.generate_camber_line()
        self.generate_thickness_profile()
        self.x_upper = self.chord_length * (self.x - self.yt * 
                                        np.sin(np.arctan(self.camber_grad)))[1:]
        self.x_lower = self.chord_length * (self.x + self.yt * 
                                        np.sin(np.arctan(self.camber_grad)))
        self.y_upper = self.chord_length * (self.camber_y + self.yt * 
                                        np.cos(np.arctan(self.camber_grad)))[1:]
        self.y_lower = self.chord_length * (self.camber_y - self.yt * 
                                        np.cos(np.arctan(self.camber_grad)))
        #want to remove resolution so the number of points is good for exporting, 
        #want to favour resolution at the LE
        #a chat gpt classis here, added 1.2 times the num points as it later removes 
        #non unique values so it will end up being less points

        # Create logarithmically spaced indices
        indices = np.linspace(1, np.log(self.num_points + 1), int(self.num_points_export*1.2))  
        # Convert back from log scale to original scale and correct offset
        values = np.exp(indices).astype(int) - 1
        # Ensure unique values
        keep_indices = sorted(set(values))  
        self.x_upper = [self.x_upper[i] for i in sorted(keep_indices) if i < self.num_points]
        self.x_lower = [self.x_lower[i] for i in sorted(keep_indices) if i < self.num_points]
        self.y_upper = [self.y_upper[i] for i in sorted(keep_indices) if i < self.num_points]
        self.y_lower = [self.y_lower[i] for i in sorted(keep_indices) if i < self.num_points]
        return
        
        
    def plot_airfoil(self):
        '''
        plots geerated airfoils :D
        '''
        
        self.generate_airfoil()
        plt.plot(np.append(self.x_lower[::-1], self.x_upper), 
                 np.append(self.y_lower[::-1], self.y_upper))
        plt.axis('equal')
        return
    
    def export_airfoil(self, name, xfoil=False):
        '''
        Exports airfoils to a txt file with a name specified by the user
        if xfoils is False it will export in a tab delimited format with
        points best suited for SOLIDWORKS useage. With the airfoil oriented
        such that it is normal to the y axis with the z axis as up
        
        if xfoils is True it will export in a space delimited file format
        and will be modified to best work with xfoils
        '''
        
        #curve will appear smoother in solidworks than in python :D
        #start at TE of airfoil so it has a point in solidworks
        self.generate_airfoil()
        x_col = np.append(self.x_lower[::-1], self.x_upper)
        if xfoil == False:
            y_col = np.zeros_like(x_col)
            z_col = np.append(self.y_lower[::-1], self.y_upper)
            formatted = np.array([x_col, y_col, z_col], dtype=float).transpose()
        else:
            #xfoils needs airfoils to be normalized and needs the start and end points to de different
            y_col = np.append(self.y_lower[::-1], self.y_upper)
            formatted = np.array([x_col[:-1], y_col[:-1]], dtype=float).transpose() / self.chord_length
        np.savetxt(f'{name}.txt', formatted, fmt='%f')   # X is an array
        return