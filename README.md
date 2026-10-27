# Algoritmo Genético - Función de Rastrigin 2D

## Descripción

Implementación en lenguaje C de un algoritmo genético con codificación binaria para minimizar la función de Rastrigin en dos dimensiones.

La función objetivo es:

\[
f(x,y)=20+x^2-10\cos(2\pi x)+y^2-10\cos(2\pi y)
\]

El dominio de búsqueda es:

- \(x \in [-5.12, 5.12]\)
- \(y \in [-5.12, 5.12]\)

El mínimo global teórico es \(f(0,0)=0\).

## Características

- Codificación binaria de 32 bits por variable.
- Población inicial aleatoria.
- Selección aleatoria de padres.
- Cruza de cromosomas mediante máscaras de bits.
- Mutación mediante inversión de un bit.
- Evaluación y ordenamiento por aptitud.
- Conservación de los individuos más aptos.
- Registro de la convergencia en un archivo CSV.

## Requisitos

- GCC.
- Biblioteca matemática de C.

## Compilación

```bash
gcc -Wall -Wextra -std=c11 -O2 ag_rastrigin.c -o ag -lm
```

## Ejecución

```bash
./ag
```

También se pueden proporcionar los siguientes argumentos:

```bash
./ag [semilla] [n_individuos] [n_generaciones] [p_cruza] [p_mutacion]
```

Ejemplo:

```bash
./ag 7 200 100 0.9 0.1
```

## Resultados

El programa muestra la evolución del mejor individuo encontrado y sus coordenadas. También genera `convergencia.csv`, que contiene información sobre la mejor aptitud y la aptitud promedio por generación.

## Lenguaje

C
