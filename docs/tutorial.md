# Tutorial

Bajo esta implementación, usar cuadratura de Gauss mediante polinomios de Legendre es bastante sencilla y se resume en definir la función a integrar (`integrando`) y utilizar la función `int_cuadraGauss` para obtener el resultado. A continuación se explica con detalle cada paso:

* ## Paso 1:
### Definir la función a integrar y sus límites de integración.
  Por ejemplo, se desea integrar: $f(x) = x^6 - x^{2} sin(2x)$ de $[1,3]$.

```python
def integrando(x):
	return x**6 -(x**2)*np.sin(2*x)
```


* ## Paso 2:
### Llamar a la función `int_cuadraGauss`.
```python
result = int_cuadraGauss(integrando, 1,3)
```

* ## Paso 3:
### Imprimir el resultado.
```python
print(result)
```
