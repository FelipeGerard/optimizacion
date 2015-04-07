
Clase 1
=========================

### Multiplicación de matrices en paralelo.

Sean A, B, C matrices de nxn. Se desea calcular C = AB en un cluster
con *p* procesadores (CPUs). Asumimos que los CPUs están lógicamente
distribuidos en forma de anillo de 0 a p-1.

Cada CPU se comunica exactamente con 2 CPUs de modo que el i-ésimo CPU
intercambia mensajes con el CPU "anterior" $(i-1) mod p$ y con el
"siguiente" $(i+1) mod p$.

Suponemos que p divide a n, es decir $n = mp$. Dividimos entonces a la
matriz A e una forma "block-row wise" como sigue:

A = [A1 \\ ... \\ An], con Ai = [Ai,0 ... Ai,p-1] "bloques renglón" de
mxn y Ai,j de mxm.

Dividimos a B y a C de la misma forma y asumimos que inicialmente el
i-ésimo bloque de A, B y C están en el procesador i. El i-ésimo
procesador calcula el i-ésimo bloque de C, Ci, de acuerdo a:

Ci = Ai * B = sum_{j=0}^p-1 Ai,j*Bj

->> Idea: Ai y Ci residen en el i-ésimo CPU y se transmite B por cada
   CPU (para no tener que almacenar B completa en ningún momento en un
   solo CPU).

* EJ. p = 4, nodo i:
  1) Calcular Ai,i x Bi = Di,i
  2) Calcular Ai,i-1 x Bi-1 = Di,i-1
  ...
  p-1) Calcular Ai,i+1 x Bi+1 = Di,i+1
  END)  Calcular Ci = sum_{j} Di,j
  * Obs: En cada paso, obtiene el siguiente bloque de B del CPU
  anterior Pi-1 y manda el bloque que tiene al CPU siguiente
  * Obs: Usa Ai y Bi primero, luego Ai y Bi-1 que le pasa el CPU
  anterior Pi-1 y le manda Bi a Pi+1, etc
  * Obs: Bj llega a Pi en (i-j) mod p pasos

*Pseudocódigo:*
En el procesador Pi:
Ci <- 0mxn
for j = 0, 1, ..., p-1 do
	Ci <- Ci + A[i, i-j mod p] * B[i-j mod p]
	Enviar B[i-j mod p] al CPU siguiente, P[i+1 mod p]
	Recibir B[(i-j-1) mod p] del CPU anterior, P[i-1 mod p]
end for



