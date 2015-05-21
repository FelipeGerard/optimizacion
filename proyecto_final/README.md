Proyecto final. Algoritmos de gran escala
=============================================

NOTAS DE ESTA VERSIÓN:

En la carpeta `src` se encuentra el código. En la subcarpeta `serie` está el código original de Erick. Aunque el procesamiento estaba en paralelo, le mandaba toda la información a los nodos (ie. las matrices completas). Nuestra versión es la que está en la subcarpeta `paralelo`. Igual que Erick utilizamos `parallel` para hacer el cluster (`snow` tiene una notación muy similar y de hecho `parallel` puede usar clusters creados por `snow`). En esa versión sólo mandamos lo necesario en cada caso. El código aún puede mejorar porque en esta versión mandamos la información para cada proceso (usando `clusterApply`) en lugar de una sola vez al principio. De igual forma, está más legible el código.

FALTA:

* Automatizar el intercepto, las predicciones, etc. Hacer funciones para hacerlo fácilmente sería ideal, porque si no a la hora de hacer las pruebas el proceso será muy desgastante.
* Probar la predicción en una base de datos distinta, contra otras SVM o contra `glmnet` por ejemplo.
* El reporte para Erick.
