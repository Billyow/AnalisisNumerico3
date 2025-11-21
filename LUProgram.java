import java.util.Scanner;


public class LUProgram {

    // ----------------------------------------------------------
    // SUB LuDecomp (a, b, n, tol, x, er)
    //    Versión del pseudocódigo que, dado A y b, resuelve A x = b
    //    usando descomposición LU + sustitución hacia adelante y atrás.
    // ----------------------------------------------------------
    public static void LUdecomp(double[][] a, double[] b, int n, double tol, double[] x, int[] er) {
        double[] s = new double[n]; // vector de escala
        er[0] = 0;                  // código de error (0 = OK, -1 = mal condicionado)

        // 1. Factoriza A en sus factores LU
        Decompose(a, n, tol, s, er);

        // 2. Si no hubo error, resuelve A x = b con la misma matriz a (ya transformada en LU)
        if (er[0] != -1) {
            Substitute(a, 0, n, b, x);
        }
    }

    // ----------------------------------------------------------
    // SUB Decompose(a, n, tol, s, er)
    //
    // Realiza la factorización LU de A con pivoteo parcial.
    //  - a se sobreescribe con L y U mezcladas:
    //      * parte inferior (debajo de la diagonal) = L (sin la diagonal, que es 1)
    //      * parte superior (incluyendo diagonal) = U
    //  - s: vector de escala para cada fila
    //  - tol: tolerancia para decidir si el pivote es demasiado pequeño
    //  - er: código de error (0 = OK, -1 = sistema mal condicionado)
    // ----------------------------------------------------------
    public static void Decompose(double[][] a, int n, double tol, double[] s, int[] er) {

        // 1. Cálculo del vector de escala s(i) = máximo valor absoluto de la fila i
        for (int i = 0; i < n; i++) {
            double o = 0.0;
            for (int j = 0; j < n; j++) {
                double abs = Math.abs(a[i][j]);
                if (abs > o) {
                    o = abs;
                }
            }
            s[i] = o;
        }

        // 2. Proceso de eliminación con pivoteo parcial
        for (int k = 0; k < n - 1; k++) {

            // Seleccionar fila pivote óptima a partir de k usando s
            Pivot(a, s, n, k);

            // Verificar si el pivote normalizado es demasiado pequeño
            if (Math.abs(a[k][k] / s[k]) < tol) {
                System.out.println("Pivot demasiado pequeño en la fila " + (k + 1));
                er[0] = -1; // matriz mal condicionada
                return;     // salimos sin completar la descomposición
            }

            // Eliminación hacia adelante (construcción de L y U)
            for (int i = k + 1; i < n; i++) {
                double factor = a[i][k] / a[k][k]; // L(i,k)
                a[i][k] = factor;                  // guardamos L en la parte inferior de a

                for (int j = k + 1; j < n; j++) {
                    // U(i,j) = A(i,j) - L(i,k)*U(k,j)
                    a[i][j] = a[i][j] - factor * a[k][j];
                }
            }
        }

        // Revisión final del último pivote
        if (Math.abs(a[n - 1][n - 1] / s[n - 1]) < tol) {
            er[0] = -1;
        }

        System.out.println("\nMatriz A después de la descomposición LU (L y U mezcladas):");
        printMatrix(a);
    }

    // ----------------------------------------------------------
    // SUB Pivot(a, s, n, k)
    //
    // Selecciona la mejor fila pivote entre k y n-1 usando el
    // criterio |a(i,k)/s(i)| máximo y, si es necesario, intercambia filas.
    // ----------------------------------------------------------
    public static void Pivot(double[][] a, double[] s, int n, int k) {
        int p = k; // fila candidata a pivote
        double big = Math.abs(a[k][k] / s[k]);

        // Buscar la mejor fila pivote por debajo
        for (int i = k + 1; i < n; i++) {
            double dummy = Math.abs(a[i][k] / s[i]);
            if (dummy > big) {
                big = dummy;
                p = i;
            }
        }

        // Si la mejor fila no es la actual, hacemos el intercambio de filas en A y en s
        if (p != k) {
            double[] tempRow = a[k];
            a[k] = a[p];
            a[p] = tempRow;

            double tempS = s[k];
            s[k] = s[p];
            s[p] = tempS;
        }
    }

    // ----------------------------------------------------------
    // SUB Substitute(a, o, n, b, x)
    //
    // Resuelve el sistema A x = b sabiendo que A ya está factorizada en LU.
    // Se realiza en dos pasos:
    //   1) Resolución de L y = b  (sustitución hacia adelante, guardando y en b)
    //   2) Resolución de U x = y  (sustitución hacia atrás, resultado en x)
    //
    // El parámetro 'o' existe en el pseudocódigo, pero aquí no lo usamos.
    // ----------------------------------------------------------
    public static void Substitute(double[][] a, int o, int n, double[] b, double[] x) {

        // 1) Sustitución hacia adelante: L * y = b
        //    Recordar que la diagonal de L es 1, y L está almacenada debajo de la diagonal de 'a'.
        for (int i = 1; i < n; i++) {
            double sum = b[i];
            for (int j = 0; j <= i - 1; j++) {
                // L(i,j) = a[i][j]
                sum -= a[i][j] * b[j];
            }
            b[i] = sum; // aquí b[i] se convierte en y(i)
        }

        // 2) Sustitución hacia atrás: U * x = y
        //    La parte superior de 'a' (incluida diagonal) es U.
        x[n - 1] = b[n - 1] / a[n - 1][n - 1];

        for (int i = n - 2; i >= 0; i--) {
            double sum = b[i];
            for (int j = i + 1; j < n; j++) {
                sum -= a[i][j] * x[j];
            }
            x[i] = sum / a[i][i];
        }
    }

    // ----------------------------------------------------------
    // Función auxiliar: imprime una matriz en formato legible
    // ----------------------------------------------------------
    public static void printMatrix(double[][] m) {
        for (double[] row : m) {
            for (double val : row) {
                System.out.printf("%12.6f ", val);
            }
            System.out.println();
        }
        System.out.println();
    }

    // ----------------------------------------------------------
    // PROGRAMA PRINCIPAL (el que genera la matriz inversa)
    //
    // Es la traducción del último pseudocódigo:
    //   - Llama a Decompose
    //   - Si no hay error, va formando la inversa columna por columna
    // ----------------------------------------------------------
    public static void main(String[] args) {

        Scanner sc = new Scanner(System.in);

        System.out.println("==============================================");
        System.out.println("  DESCOMPOSICIÓN LU + CÁLCULO DE LA INVERSA  ");
        System.out.println("==============================================\n");

        System.out.println("A * x = e_i para cada columna.\n");

        // Ejemplos de entrada para que el usuario se ubique
        System.out.println("Ejemplo 1 (matriz 2x2):");
        System.out.println("  n = 2");
        System.out.println("  A = [ 4   3 ]");
        System.out.println("      [ 6   3 ]\n");

        System.out.println("Ejemplo 2 (matriz 3x3):");
        System.out.println("  n = 3");
        System.out.println("  A = [ 4  -1   2 ]");
        System.out.println("      [ -2  5   1 ]");
        System.out.println("      [  1  2   3 ]\n");

        System.out.print("Ingrese el tamaño n de la matriz cuadrada A (n x n): ");
        int n = sc.nextInt();

        double[][] a = new double[n][n];   // matriz original A (se modificará a LU)
        double[][] ai = new double[n][n];  // matriz inversa A^{-1}
        double[] b = new double[n];        // vector del lado derecho (columna de la identidad)
        double[] x = new double[n];        // solución temporal x para cada columna de la inversa
        double tol = 1e-6;                 // tolerancia para pivotes pequeños
        int[] er = new int[1];             // código de error

        System.out.println("\nAhora ingrese los elementos de la matriz A fila por fila.");
        System.out.println("Por ejemplo, para el Ejemplo 2 (3x3) escribirías:");
        System.out.println("Fila 1: 4  -1   2");
        System.out.println("Fila 2: -2  5   1");
        System.out.println("Fila 3: 1   2   3\n");

        // Lectura de la matriz A
        for (int i = 0; i < n; i++) {
            System.out.println("Fila " + (i + 1) + ":");
            for (int j = 0; j < n; j++) {
                System.out.print("A[" + (i + 1) + "][" + (j + 1) + "] = ");
                a[i][j] = sc.nextDouble();
            }
        }

        System.out.println("\nMatriz A ingresada:");
        printMatrix(a);

        // Vector de escala para Decompose
        double[] s = new double[n];

        // Llamamos a Decompose (factorización LU)
        Decompose(a, n, tol, s, er);

        if (er[0] == 0) {
            // No hubo error: calculamos la inversa

            System.out.println("Calculando la matriz inversa A^{-1}...\n");

            // Recorremos cada columna r de la matriz identidad
            for (int r = 0; r < n; r++) {

                // Construimos b = e_r (vector canónico)
                for (int j = 0; j < n; j++) {
                    b[j] = (j == r) ? 1.0 : 0.0;
                }

                // Resolvemos A x = e_r usando la LU almacenada en 'a'
                Substitute(a, 0, n, b, x);

                // Guardamos el resultado x como la columna r de la inversa
                for (int j = 0; j < n; j++) {
                    ai[j][r] = x[j];
                }
            }

            System.out.println("Matriz inversa A^{-1}:");
            printMatrix(ai);


        } else {
            // Hubo error en la descomposición (pivote demasiado pequeño)
            System.out.println("El sistema está mal condicionado. No es seguro calcular la inversa.");
        }

        sc.close();
    }
}
