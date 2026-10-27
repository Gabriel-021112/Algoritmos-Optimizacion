/*
 * Algoritmo genetico con codificacion binaria - Funcion de Rastrigin 2D
 *
 *   f(x,y) = 20 + (x^2 - 10cos(2*pi*x)) + (y^2 - 10cos(2*pi*y))
 *   Dominio: x, y en [-5.12, 5.12]      Minimo global: f(0,0) = 0
 *
 * Individuo:  unsigned int xy[2]  +  f(x)
 *   cromo[0] = 32 bits de x      cromo[1] = 32 bits de y
 *
 * Decodificacion:  X = A + (ValorDecimal / (2^L - 1)) * (B - A),  L = 32
 *
 * Pasos del algoritmo (pizarron):
 *   1. Definir individuo
 *   2. Generar poblacion aleatoria (# individuos, % cruza, % mutacion, # generaciones)
 *   3. Cruza: seleccionar 2 individuos de forma aleatoria
 *   4. Mutacion a hijos
 *   5. Evaluar f(x) para toda la poblacion
 *   6. Ordenar poblacion (eliminar menos aptos)      -> se repite # generaciones
 *
 * Compilar:  gcc ag_rastrigin.c -o ag -lm
 * Ejecutar:  ./ag [semilla] [n_individuos] [n_generaciones] [p_cruza] [p_mutacion]
 *   Ej.:     ./ag 7 200 100 0.9 0.1   (valores por defecto salvo la semilla)
 * Genera tambien "convergencia.csv" (generacion, mejor f, f promedio) para graficar.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#define L_BITS      32                 /* bits por variable               */
#define LIM_INF     (-5.12)            /* A                               */
#define LIM_SUP     (5.12)             /* B                               */
#define MAX_VALOR   4294967295.0       /* 2^L - 1                         */
#define PI          3.14159265358979323846

/* ---------------- Paso 1: definir individuo ---------------- */
typedef struct {
    uint32_t cromo[2];   /* cromo[0] = x (32 bits), cromo[1] = y (32 bits) */
    double   x, y;       /* valores reales decodificados                   */
    double   fx;         /* f(x,y): mientras menor, mas apto               */
} Individuo;

typedef struct {
    int    n_ind;        /* numero de individuos   */
    int    n_gen;        /* numero de generaciones */
    double p_cruza;      /* porcentaje de cruza    */
    double p_mut;        /* porcentaje de mutacion */
} Parametros;

/* ---------------- Generador pseudoaleatorio (xorshift64*) ----------------
 * Se usa uno propio para tener 64 bits buenos y resultados iguales en
 * cualquier sistema (rand() cambia de un compilador a otro).              */
static uint64_t estado_rng;

static uint64_t siguiente64(void) {
    estado_rng ^= estado_rng >> 12;
    estado_rng ^= estado_rng << 25;
    estado_rng ^= estado_rng >> 27;
    return estado_rng * 2685821657736338717ULL;
}

static double aleatorio01(void) {              /* uniforme en [0, 1) */
    return (siguiente64() >> 11) * (1.0 / 9007199254740992.0);
}

static int entero_hasta(int n) {               /* entero en [0, n-1] */
    return (int)(aleatorio01() * n);
}

/* ---------------- Decodificacion y funcion objetivo ---------------- */
static double decodificar(uint32_t valor) {
    return LIM_INF + ((double)valor / MAX_VALOR) * (LIM_SUP - LIM_INF);
}

static double rastrigin(double x, double y) {
    return 20.0 + (x * x - 10.0 * cos(2.0 * PI * x))
                + (y * y - 10.0 * cos(2.0 * PI * y));
}

/* Paso 5: evaluar f(x) de un individuo */
static void evaluar(Individuo *ind) {
    ind->x  = decodificar(ind->cromo[0]);
    ind->y  = decodificar(ind->cromo[1]);
    ind->fx = rastrigin(ind->x, ind->y);
}

/* Paso 2: individuo con cadena de 64 bits aleatoria */
static Individuo individuo_aleatorio(void) {
    Individuo ind;
    uint64_t r = siguiente64();
    ind.cromo[0] = (uint32_t)(r >> 32);
    ind.cromo[1] = (uint32_t)(r & 0xFFFFFFFFu);
    evaluar(&ind);
    return ind;
}

/* ---------------- Paso 3: cruza con mascaras ----------------
 * En cada variable (x e y) se elige un punto de corte k al azar (1..31).
 *   mascara_baja = 0...01...1  (k unos)      mascara_alta = ~mascara_baja
 *   hijo1 = (padre1 & mascara_alta) | (padre2 & mascara_baja)
 *   hijo2 = (padre2 & mascara_alta) | (padre1 & mascara_baja)
 * Asi los hijos mezclan bits DENTRO de x y de y (no solo intercambian x con y). */
static void cruzar(const Individuo *p1, const Individuo *p2,
                   Individuo *h1, Individuo *h2) {
    for (int g = 0; g < 2; g++) {
        int k = 1 + entero_hasta(L_BITS - 1);
        uint32_t mascara_baja = (1u << k) - 1u;
        uint32_t mascara_alta = ~mascara_baja;
        h1->cromo[g] = (p1->cromo[g] & mascara_alta) | (p2->cromo[g] & mascara_baja);
        h2->cromo[g] = (p2->cromo[g] & mascara_alta) | (p1->cromo[g] & mascara_baja);
    }
}

/* ---------------- Paso 4: mutacion ----------------
 * Con probabilidad p_mut se invierte (XOR) un bit al azar de los 64. */
static void mutar(Individuo *ind, double p_mut) {
    if (aleatorio01() < p_mut) {
        int bit = entero_hasta(2 * L_BITS);          /* 0..63 */
        ind->cromo[bit / L_BITS] ^= 1u << (bit % L_BITS);
    }
}

/* ---------------- Paso 6: ordenar de menor a mayor f(x) ---------------- */
static int por_aptitud(const void *a, const void *b) {
    double fa = ((const Individuo *)a)->fx;
    double fb = ((const Individuo *)b)->fx;
    return (fa > fb) - (fa < fb);
}

static void imprimir_bits(uint32_t v) {
    for (int i = L_BITS - 1; i >= 0; i--)
        putchar(((v >> i) & 1u) ? '1' : '0');
}

static double promedio_fx(const Individuo *pob, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += pob[i].fx;
    return s / n;
}

int main(int argc, char *argv[]) {
    Parametros p = { 200, 100, 0.9, 0.1 };   /* 0.9 y 0.1 vienen del pizarron */
    uint64_t semilla = (argc > 1) ? (uint64_t)strtoull(argv[1], NULL, 10)
                                  : (uint64_t)time(NULL);
    if (argc > 2) p.n_ind  = atoi(argv[2]);
    if (argc > 3) p.n_gen  = atoi(argv[3]);
    if (argc > 4) p.p_cruza = atof(argv[4]);
    if (argc > 5) p.p_mut   = atof(argv[5]);

    if (p.n_ind < 2 || p.n_gen < 1) {
        fprintf(stderr, "Error: se requieren n_individuos >= 2 y n_generaciones >= 1\n");
        return 1;
    }

    estado_rng = semilla * 0x9E3779B97F4A7C15ULL + 0x1234567ULL;
    if (estado_rng == 0) estado_rng = 88172645463325252ULL;
    for (int i = 0; i < 10; i++) siguiente64();      /* calentar el generador */

    /* padres en [0, n_ind), hijos en [n_ind, 2*n_ind) */
    Individuo *pob = malloc(sizeof(Individuo) * 2 * p.n_ind);
    if (!pob) { fprintf(stderr, "Sin memoria\n"); return 1; }

    FILE *csv = fopen("convergencia.csv", "w");
    if (csv) fprintf(csv, "generacion,mejor_fx,promedio_fx\n");

    printf("AG Rastrigin 2D | individuos=%d generaciones=%d cruza=%.3f mutacion=%.3f semilla=%llu\n\n",
           p.n_ind, p.n_gen, p.p_cruza, p.p_mut, (unsigned long long)semilla);

    /* Paso 2: poblacion inicial aleatoria */
    for (int i = 0; i < p.n_ind; i++) pob[i] = individuo_aleatorio();
    qsort(pob, p.n_ind, sizeof(Individuo), por_aptitud);

    printf("%-6s %-12s %-12s %-12s\n", "Gen", "x", "y", "f(x,y)");
    printf("%-6d %-12.6f %-12.6f %-12.6f\n", 0, pob[0].x, pob[0].y, pob[0].fx);
    if (csv) fprintf(csv, "0,%.10f,%.10f\n", pob[0].fx, promedio_fx(pob, p.n_ind));

    int paso_impresion = (p.n_gen >= 10) ? p.n_gen / 10 : 1;

    for (int gen = 1; gen <= p.n_gen; gen++) {
        /* Pasos 3 y 4: generar n_ind hijos */
        for (int n = 0; n < p.n_ind; n += 2) {
            const Individuo *p1 = &pob[entero_hasta(p.n_ind)];   /* 2 al azar */
            const Individuo *p2 = &pob[entero_hasta(p.n_ind)];
            Individuo h1 = *p1, h2 = *p2;                        /* sin cruza: copias */

            if (aleatorio01() < p.p_cruza) cruzar(p1, p2, &h1, &h2);

            mutar(&h1, p.p_mut);
            mutar(&h2, p.p_mut);

            pob[p.n_ind + n] = h1;
            if (n + 1 < p.n_ind) pob[p.n_ind + n + 1] = h2;
        }

        /* Paso 5: evaluar f(x) de los hijos */
        for (int i = p.n_ind; i < 2 * p.n_ind; i++) evaluar(&pob[i]);

        /* Paso 6: ordenar padres + hijos y conservar los n_ind mejores */
        qsort(pob, 2 * p.n_ind, sizeof(Individuo), por_aptitud);

        if (csv) fprintf(csv, "%d,%.10f,%.10f\n", gen, pob[0].fx, promedio_fx(pob, p.n_ind));
        if (gen % paso_impresion == 0 || gen == p.n_gen)
            printf("%-6d %-12.6f %-12.6f %-12.6f\n", gen, pob[0].x, pob[0].y, pob[0].fx);
    }

    printf("\nMejor individuo encontrado:\n");
    printf("  x (bits) = "); imprimir_bits(pob[0].cromo[0]);
    printf("\n  y (bits) = "); imprimir_bits(pob[0].cromo[1]);
    printf("\n  x = %.8f\n  y = %.8f\n  f(x,y) = %.10f   (optimo teorico: 0 en (0,0))\n",
           pob[0].x, pob[0].y, pob[0].fx);

    if (csv) { fclose(csv); printf("\nConvergencia guardada en convergencia.csv\n"); }
    free(pob);
    return 0;
}
