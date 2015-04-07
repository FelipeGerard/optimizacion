Clase 2
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
con $U$ de mxm, $V$ de nxn matrices unitarias y $\Sigma$ diagonal de mxn. Aquí $M^H$ denota la transpuesta conjugada de M; si M tiene entradas reales entonces $M^H = M^T$. Como $A^H = V \Sigma^T U^H$, asumimos s.p.g. que $m \geq n$ y que $r = rango(A) \leq n = min(m,n)$. Denotamos $u_j$, $v_j$ las columnas de $U$ (vectores singulares izquierdos de $A$) y de $V$ (vectores singulares derechos de $A$) respectivamente y $\sigma_i, i = 1, ..., n$ los valores singulares de A (reales y no negativos). Típicamente las entradas de $\Sigma$ y las columnas de $U$ y de $V$ se ordenan de acuerdo a $\sigma_1 \geq ... \geq \sigma_r > 0$, con $\sigma_{r+1} = ... = \sigma_n = 0$.

### Algoritmo de Jacobi

El algoritmo de Jacobi aplica al caso de matrices reales y se extienden también al caso complejo. Trabajaremos con matrices reales. Bajo este supuesto, las matrices $U$ y $V$ son reales también (y por lo tanto ortogonales). Si A es simétrica se pueden ajustar los signos de las entradas de $\Sigma$ de modo que  $U = V$ y $A = U^T D U$ corresponde con la descomposición espectral.

Hay dos tipos de métodos basados en Jacobi: one-sided y two-sided.

**Método de Jacobi: one-sided (MJOS)**

Idea: Construimos una matriz $V$ ortogonal tal que $AV = W$ tenga columnas ortogonales. Normalizando las columnas de $W$ distintas de cero (norma euclidiana), tenemos que
$$
	W = [U_r 0]\begin{pmat} \Sigma_r & 0 \\ 0 & 0 \end{pmat},
$$
con $U_r$ de mxr las columnas ortonormales distintas de cero y $\Sigma_r = diag(\sigma_1, ..., sigma_r)$ las normas de dichas columnas. Una SVD para A está dada por
$$
	A = U_r \Sigma_r V_r^T,
$$
donde $V_r$ de nxr tiene r columnas de V, que es de nxn.

*Algoritmo:*

A_0 = A
V_0 = I_n
for k > 0:
	A_{k+1} = A_k U_k
	V_{k+1} = V_k U_k,
	donde U_k son matrices de rotación del plano para una (i,j) predefinida para cada k, es decir, una matriz identidad pero con elementos:
	U_k[i,i] = U_k[j,j] = cos(theta), U_k[i,j] = sen(theta), U_k[j,i] = -sen(theta), i < j.
	La multiplicación de A_k U_k sólo afecta a 2 columnas de A_k: (A_{k+1}[,i], A_{k+1}[,j]) = (A_k[,i], A_k[,j]) |cos(theta) sen(theta)|
														      | -sen(theta) cos(theta)|
	El ángulo theta se elige:
		1) 0 si A_k[,i]^T A_k[,j] = 0, ie. no hacemos rotación
		2) theta en (-pi/4, pi/4) tal que A_{k+1}[,i]^T A_{k+1}[,j] = 0. Para ello calculamos
			* cos(theta) = 1/sqrt(1 + t^2), sen(theta) = t*cos(theta), donde t = signo(xi)/(abs(xi) + sqrt(1+xi^2)) y xi = (A_k[,j]^T A_k[,j] - A_k[,i]^T A_k[,i])/(2*A_k[,i]^T A_k[,j])

*Observaciones:*
	* Para la implementación, usar (A_k[,i]^T A_k[,j]) / (norma(A_k[,i])*norma(A_k[,j])) < tol para la condición de ortogonalidad de A_k[,i]^T A_k[,j] = 0, con tol <= 10^-8.
	* Sin 1/xi <= eps_maq, entonces theta = 0; en otro caso, calculamos t, xi y U_k normalmente.
	* | A_{k+1}[,i] A_{k+1}[,j] | = | A_k[,i] A_k[,j] | U_k, es decir, sólo actualizamos lo afectado en cada paso.
	  | V_{k+1}[,i] V_{k+1}[,j] |   | V_k[,i] V_k[,j] |






































