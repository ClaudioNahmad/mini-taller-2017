#recomendacion: bajar el editor de texto Sublime

#I. Seccion cosmologica del codigo

import numpy as np
import scipy.integrate as integrate

#utilizando las ecuaciones que estan en el cuaderno:

def integrando(z,omega_lambda):
	return 1/np.sqrt((omega_lambda+(1-omega_lambda)*(1+z)**3)) 

def dl(z,omega_lambda,h):
	integral=integrate.quad(integrando,0,z,args=(omega_lambda)) #la sintaxis es: (a quien integra, de donde, a donde, que argumentos no va a integrar)
	c=300000
	H0=100*h
	return c/H0*(1+z)*integral[0] #sintaxis: integral[0] solo devuelve la primera entrada de la funcion

#vamos a leer los valores de mu y z del archivo mu_vs_z
zmu=np.loadtxt('SCPUnion2_mu_vs_z.txt',usecols=(1,2),skiprows=5)
covmat=np.loadtxt('SCPUnion2_covmat_sys.txt')

def loglike(entrada):
	omega_lambda=entrada[0]
	h=entrada[1]
	distancia=dl(zmu[:,0],omega_lambda,h)
	mu_teo=5*np.log10(distancia)+25
	mu_obs=zmu[:,1]	
	delta=mu_teo-mu_obs
	chisquare=np.dot(np.dot(delta,inv_covmat),delta)
	return -chisquare/2

dl=np.vectorize(dl)
inv_covmat=np.linalg.inv(covmat)




print integrando(2,.5)		#prueba1
print dl (3.,.5,.7)		#prueba2
print zmu[1,1]			#prueba3, pide la componente [1,1] de la matriz de mu_vs_z
print loglike([.5,.5])

#II. Seccion MCMC del codigo

pasoinicial=[0.7,0.7]
cadena=[pasoinicial]
cadenalike=[loglike(cadena[0])] #este renglon guarda las entradas de loglike(cadena[0]) en un archivo
ancho_omega=0.01
ancho_h=0.01

n=10
for i in range(n): #range(10) es el numero de pasos (i.e. 10)
	aleatorio=np.random.normal(0,1,2)
	nuevopunto=cadena[i]+[ancho_omega,ancho_h]*aleatorio #aqui nos estamos moviendo a un nuevo punto, la cadena se mueve del 1 al 10 y en cada paso se cambia el ancho
	nuevolike=loglike(nuevopunto)
	if nuevolike > cadenalike[i]:
		aceptar=1
	else:
		aceptar=np.exp(nuevolike-cadenalike[i]) #notemos que nuevolike-cadenalike es un numero negativo
	
	if aceptar >=np.random.uniform(0.,1.):
		cadena.append(nuevopunto)
		cadenalike.append(nuevolike)
	else:
		cadena.append(cadena[i])
		cadenalike.append(cadenalike[i])


cadena=np.array(cadena) #esto convierte la lista en un array

omegaestimada=np.mean(cadena[:,0])
hestimada=np.mean(cadena[:,1])
cadenaomega=cadena[:,0]
sigmaomega=np.sqrt(np.mean(np.square(cadenaomega))-omegaestimada**2)

print 'con',n,'pasos'
print 'omega_lambda=', omegaestimada, 'h=', hestimada
print 'sigma_omega=', sigmaomega

