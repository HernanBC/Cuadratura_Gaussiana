import numpy as np


def gaussxw(N):

	"""
    Calcula los puntos de muestreo y los pesos de la cuadratura de Gauss.

    La cuadratura de Gauss se utiliza para aproximar integrales definidas
    mediante una suma ponderada de valores de la función en puntos específicos.

    Args:
		N (int): Número de puntos de muestreo y pesos a calcular, corresponde al número de subregiones a usar.

    Returns:
		Output (tuple):
			x (numpy.ndarray): Puntos de muestreo de la cuadratura de Gauss.
			w (numpy.ndarray): Pesos correspondientes a cada punto de muestreo.
	
	Example:
		>>> x, w = gaussxw(5)
		>>> print(x)
		>>> print(w)
    """

	# SE OPTA POR MANTENER LA RUTINA TAL Y COMO SE ENCUENTRA EN EL JUPYTER ORIGINAL.
	# Aproximación inicial
	a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
	# Note código vectorial aquí
	x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

	# Ahora calculamos las raíces de los polinomios utilizando el método de Newton
	# Este es un tema que veremos la próxima semana!
	# De momento, puede ignorar el siguiente flujo de control con el "while" y saber que esto 
	# devuelve los puntos de muestreo obtenidos con los polinomios de Legendre
	epsilon = 1e-15
	delta = 1.0
	while delta > epsilon:
		p0 = np.ones(N, dtype = float)
		# Deep copy
		p1 = np.copy(x)
		for k in range(1, N):
			p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
		dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
		dx = p1 / dp
		x -= dx
		delta = np.max(np.abs(dx))

	# Ahora calculamos los pesos
	w = 2 * (N + 1) * (N + 1)/(N * N * (1 - x * x) * dp * dp)

	# Note que la función devuelve un tuple
	return x,w

def gaussxwab(a, b, x, w):
	"""
    Transforma los puntos y pesos de la cuadratura de Gauss
    desde el intervalo [-1, 1] al intervalo [a, b].

    Args:
		a (float): Límite inferior del nuevo intervalo de integración.
		b (float): Límite superior del nuevo intervalo de integración.
		x (numpy.ndarray): Puntos de muestreo en el intervalo [-1, 1].
		w (numpy.ndarray): Pesos asociados a los puntos de muestreo en [-1, 1].

    Returns:
		Output (tuple):
			x_escalado (numpy.ndarray): Puntos de muestreo transformados al intervalo [a, b].
			w_escalado (numpy.ndarray): Pesos transformados correspondientes a los nuevos puntos de muestreo en el intervalo [a, b].

	Example:
		>>> x_escalado, w_escalado = gaussxw(1,3,x,w)
		>>> print(x_escalado)
		>>> print(w_escalado)
    """
	# SE OPTA POR MANTENER LA RUTINA TAL Y COMO SE ENCUENTRA EN EL JUPYTER ORIGINAL.
	# Obtenido de pag 168 Newman)
	# Note que esta función está escrita a-la-Pitón: si funciona bueno; si no, es culpa suya.
	return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def int_cuadraGauss(integrando, x_inf, x_sup, steps, noEscalado = gaussxw, escalado = gaussxwab ):
	"""
    Calcula la integral de una función utilizando la cuadratura de Gauss con polinomio ortogonal (Legendre).

    Esta función evalúa una integral definida en el intervalo [x_inf, x_sup]
    mediante la cuadratura de Gauss, escalando los puntos y pesos de la cuadratura 
    de Gauss mediante polinomios ortogonales al intervalo deseado. Por defecto se utiliza 
	polinomios de Legendre.

    Args:
		integrando (callable): Función a integrar. Debe aceptar un array de puntos como entrada.
		x_inf (float): Límite inferior del intervalo de integración.
		x_sup (float): Límite superior del intervalo de integración.
		steps (int): Número de puntos de muestreo a utilizar en la cuadratura.
		noEscalado (callable): Función que calcula los puntos y pesos en el intervalo [-1, 1]. Por defecto, utiliza `gaussxw`.
		escalado (callable): Función que transforma los puntos y pesos al intervalo [x_inf, x_sup]. Por defecto, utiliza `gaussxwab`.

    Returns:
		output (float): Aproximación de la integral definida de `integrando` en el intervalo [x_inf, x_sup].
	
	Example:
		>>> # Integrar x**6 -(x**2)*sin(2*x) en el intervalo [1,3]
		>>> def integrando(x):
		>>>		return x**6 -(x**2)*np.sin(2*x)

		>>> int_cuadraGauss(integrando,1,3)
    """
	
	x_sin_escalar, w_sin_escalar = noEscalado(steps)
	x , w = escalado(x_inf,x_sup,x_sin_escalar,w_sin_escalar)

	I = np.dot(integrando(x),w)
	return I

# Resolución del problema solicitado
# Integrar x**6 -(x**2)*sin(2*x) en el intervalo [1,3]

# Definiendo la expresión a integrar como una función de una variable independiente x
def integrando(x):
	return x**6 -(x**2)*np.sin(2*x)

N = [2,3,4,5,6,7] # Número de regiones que se van a utilizar para integración

for i in N:
	print(f"Resultado de la integral: {int_cuadraGauss(integrando,1,3,i)} (N = {i})") # Printeo de la solución al usar cuadratura de Gauss y el respectivo N

# El valor exacto se obtiene con N = 7