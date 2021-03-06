Managed to simulate a simple pendulum of unit mass and of length l = g using gsl_odeiv2. Plotted theta and d(theta)/dt against time for 1000 oscillations for a small initial displacement and got sin(wt) or cos(wt) as expected. Investigated how well solution conserved energy, found it doesn't conserve energy as well as I would expect, with 0.5% change after 10,000 oscillations. Plotted how period varies with initial displacement theta0, and found when theta0 = pi/2, T =  7.4144 s

Turned on damping (q) and found that system behaved as expected; underdamped for q = 1, close to critical damping for q = 5 and overdamped for q = 10.

Now setting q = 0.5 and also turning on driving force (F), for a small force (F = 0.5) the period of the oscillations matches the period of the forced oscillations (3*pi s in this case). Of note is that for very large values of theta0, the simulation breaks and gives a period of -nan, caused due to the pendulum going 'over the top'.

For larger values of F (F = 1.44, 1.465), the program can't handle the fact that the pendulum is going past theta = pi, and the solutions stop making physical sense; the amplitude continues to decrease. For F = 1.2 we get wildly varying behaviour, with some sinusoidal trend, though it diverges for large t.

The solution is incredibly sensitive to initial conditions. Varying the initial displacement by 0.00001 completes changes the solution for large t. This is an example of a chaotic system.
