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
