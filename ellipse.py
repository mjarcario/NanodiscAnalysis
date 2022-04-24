import gc
import sys
import numpy as np

def fit_ellipse(x, y):
	#print('Fitting ellipse')
	# Converts the data passed from the Tcl script
	# into a useable format in this python script
	x_float = [float(x_val) for x_val in x.split()]
	y_float = [float(y_val) for y_val in y.split()]
	x_array = np.asarray(x_float)
	y_array = np.asarray(y_float)
	D1 = np.vstack([x_array**2, x_array*y_array, y_array**2]).T
	D2 = np.vstack([x_array, y_array, np.ones(len(x_array))]).T
	S1 = D1.T @ D1
	S2 = D1.T @ D2
	S3 = D2.T @ D2
	T = -np.linalg.inv(S3) @ S2.T
	M = S1 + S2 @ T
	C = np.array(((0, 0, 2), (0, -1, 0), (2, 0, 0)), dtype=float)
	M = np.linalg.inv(C) @ M
	eigenval, eigenvec = np.linalg.eig(M)
	con = 4 * eigenvec[0] * eigenvec[2] - eigenvec[1]**2
	ak = eigenvec[:, np.nonzero(con > 0)[0]]
	fits = np.concatenate((ak,T @ ak)).ravel()
	del x
	del y
	#del x_float
	#del y_float
	#del x_array
	#del y_array
	del D1
	del D2
	del S1
	del S2
	del S3
	del T
	del M
	del C
	del ak
	gc.collect()

	a = fits[0]
	b = fits[1] /2
	c = fits[2]
	d = fits[3] / 2
	f = fits[4] / 2
	g = fits[5]

	den = b**2 - ( a * c )
	if den > 0:
		raise ValueError('coeffs do not represent an ellipse')

	numerator = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
	fac = np.sqrt((a - c)**2 + 4*b**2)
	a1 = np.sqrt(numerator / (den * (fac - a - c)))
	a2 = np.sqrt(numerator / (den * (-fac - a - c)))
	#print("Axes are : %.2f and %.2f" % (a1, a2))
	return str(a1) + " " + str(a2)

if __name__ == '__main__':
	#x, y = np.loadtxt('nanodiscCoordinates.txt', usecols=(0,1), unpack=True)
	#axes = fit_ellipse(x, y)
	axes = fit_ellipse(*sys.argv[1:])
	print(axes)
