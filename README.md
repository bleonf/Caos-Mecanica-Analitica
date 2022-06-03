# Caos-Mecanica-Analitica
Codigo en mathematica y python que se utilizó para generar secciones de Poincare, graficar trayectorias y calcular exponenetes de Lyapunov. Instrucciones de uso en archivo README.
El archivo de python esta escrito para python3 y puede no ser compatible con versiones anteriores. 

Para utilizar ambos archivos se deben introducir valores iniciales para la energía y para el momento inicial en el eje y. 
  Valores recomendados:
  E0=-1
  Pyi=1
 
Despues de esto se debe introducir el winding number (decimal o fraccion) que se quiera utilizar. 

Para el codigo en python: introducir numero de cuerpos que se quieran simular, y el valor de step entre cada condicion inicial de los cuerpos. Cambiar el parametro "Save" entre True y False para que se guarde la figura de trayectorias o que solo se muestre (plt.show). El valor "lim" originalmente es el maximo y minimo valor en el eje Y que se quiera graficar.
 
Para el codigo en Mathematica: se debe poner el limite inferior y superior al rededor del eje X para el cual se calcularan las scciones de poincare. De estos limites y el salto que se propone sobre la variable i se determina cuantos cuerpos seran simulados.


