Clases 2 y 3
===================

Multiplicación de matrices
-----------------------------

* Para "3) EJ DGEMM", instalar también por si acaso:
	~ gfortran
	~ liblapack-dev (LAPACK)

Tarea
---------------------------

* Leer basic C
* Checar cómo usar una sola instalación
* Investigar gather
* Investigar DGEMM (dgemm_(...))
* Comparar resultados de nuestra multiplicación de matrices vs DGEMM de LAPACK
* Terminar código para multiplicar matrices en paralelo usando MPI + DGEMM

SVD
--------------------

Si $A \in \mathbb{C}^{mxn}$ entonces:
$$
	A = U \Sigma V^H,
$$
con $U$ de mxm, $V$ de nxn matrices unitarias y $\Sigma$ diagonal de mxn. Aquí $M^H$ denota la transpuesta conjugada de M; si M tiene entradas reales entonces $M^H = M^T$. Como $a^H = V \Sigma^T U^H$, asumimos s.p.g. que $m \geq n$ y que $r = rango(A) \leq n = min(m,n)$. Denotamos $u_j$, $v_j$ las columnas de $U$ (vectores singulares izquierdos de $A$) y de $V$ (vectores singulares derechos de $A$) respectivamente y $\sigma_i, i = 1, ..., n$ los valores singulares de A (reales y no negativos). Tí\picamente las entradas de $\Sigma$ y las columnas de $U$ y de $V$ se ordenan de acuerdo a $\sigma_1 \geq ... \geq \sigma_r > 0$, con $\sigma_{r+1} = ... = \sigma_n = 0$.

### Algoritmo de Jacobi

El algoritmo de Jacobi aplica al caso de matrices reales y se extienden también al caso complejo. Trabajaremos con matrices reales. Bajo este supuesto, las matrices $U$ y $V$ son reales también (y por lo tanto ortogonales). Si A es simétrica se pueden ajustar los signos de las entradas de $\Sigma$ de modo que  $U = V$ y $A = U^T D U$ corresponde con la descomposición espectral.

Hay dos tipos de métodos basados en Jacobi: one-sided y two-sided.

**Método de Jacobi: one-sided (MJOS)**

Idea: Construimos una matriz $V$ ortogonal tal que $AV = W$ tenga columnas ortogonales. Normalizando las columnas de $W$ distintas de cero (norma euclidiana), tenemos que

$$
	W = \left(U_r 0\right) \left(\begin{array}{cc} \Sigma_r & 0 \\ 0 & 0 \end{array} \right),
$$

con $U_r$ de mxr las columnas ortonormales distintas de cero y $\Sigma_r = diag(\sigma_1, ..., sigma_r)$ las normas de dichas columnas. Una SVD para A está dada por

$$
	A = U_r \Sigma_r V_r^T,
$$
donde $V_r$ de nxr tiene r columnas de V, que es de nxn.

*Algoritmo:*

| $a^0 = A$
| $V^0 = I_n$

for $k > 0$:

* $A^{k+1} = A^k U^k$
* $V^{k+1} = V^k U^k$, donde $U^k$ son matrices de rotación del plano para una $(i,j)$ predefinida para cada $k$, es decir, una matriz identidad pero con elementos:
	+ $U^k_{i,i} = U^k_{j,j} = cos(\theta), U^k_{i,j} = sen(\theta), U^k_{j,i} = -sen(\theta), i < j$
	+ La multiplicación de $a^k U^k$ sólo afecta a 2 columnas de $a^k$: $(a^{k+1}_i, a^{k+1}_j) = (a^k_i, a^k_j) \left(\begin{array}{cc} cos(\theta) & sen(\theta) \\ -sen(\theta) & cos(\theta) \end{array}\right) = (a^k_i, a^k_j)\hat(U_k)$
	+ El ángulo \theta se elige:
		- 0 si $(a^k_i)^T a^k_j = 0$, ie. no hacemos rotación
		- $\theta$ en $(-\pi/4, \pi/4)$ tal que $(a^{k+1}_i)^T a^{k+1}_j = 0$. Para ello calculamos:

			$cos(\theta) = 1/sqrt(1 + t^2)$, $sen(\theta) = t \times cos(\theta)$, donde $t = signo(\xi)/(|\xi| + sqrt(1 +(?) \xi^2))$ y $\xi = ((a^k_j)^T a^k_j - (a^k_i)^T a^k_i)/(2*(a^k_i)^T a^k_j)$

*Observaciones:*

* Para la implementación, usar $((a^k_i)^T a^k_j) / (\|a^k_i\| \|a^k_j\|) < tol$ para la condición de ortogonalidad de $(a^k_i)^T a^k_j = 0$, con $tol \leq 10^{-8}$.
* Sin $1/\xi \leq \epsilon_m$, entonces $\theta = 0$; en otro caso, calculamos $t,\xi y U^k$ normalmente.
* $\left(\begin{array}{cc} a^{k+1}_i & a^{k+1}_j\\ v^{k+1}_i & v^{k+1}_j \end{array}\right) = \left(\begin{array}{cc} a^k_i & a^k_j\\ v^k_i & v^k_j\end{array}\right)  U^k$, es decir, sólo actualizamos lo afectado en cada paso.
* Usar un contador $rot$ que se incrementa cuando una pareja de columnas sean ortogonales. El algoritmo termina si $rot = n(n-1)/2$
* En los algoritmos tradicionales de Jacobi one-sided las rotaciones se realizan en una secuencia **fija** llamada "sweep". Cada sweep consiste en $n(n-1)/2$ rotaciones y en cada sweep se ortogonalizan 2 columnas. El algoritmo termina si en un sweep todas las columnas son ortogonales. 


Clase 3
---------------------------

**Ordenamiento cíclico por renglones**

Un _sweep_ incolucra los pares $(1,2), ..., (1,n),(2,3), ..., (2,n), ..., (n-1,n)$. Hay $\binom nk = n(n-1)/2$ pares. Una desventaja es que el algoritmo de Jacobi con este ordenamiento es secuencial, pero una ventaja es que siempre converge si $-\pi/4 \leq \theta \leq \pi/4$ y además la convergencia es cuadrática. 

**Ordenamiento _round-robin_**

+) Permite ejecición de Jacobi en paralelo.
+) Genera $n(n-1)/2$ pares de columnas en $n-1$ pasos con $n/2$ procesadores. Si $n$ es impar, se agrega una columna de ceros.

_Ejemplo:_ $n = 8, p = 4$

1) $\begin{array}(c|c|c|c) CPU1&CPU2&CPU3&CPU4\\ \hline \\col1&col3&col5&col7\\ col2&col4&col6&col8 \end{array} 
2) $\begin{array}(c|c|c|c) CPU1&CPU2&CPU3&CPU4\\ \hline \\col1&col2&col3&col5\\ col4&col6&col8&col7 \end{array} 
3) etc.

_Ventajas:_

* En paralelo

_Desventajas_:

* Convergencia?
* Se extiende?







































